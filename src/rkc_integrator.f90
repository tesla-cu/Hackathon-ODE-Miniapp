module rkc_integrator
!< The error-controlled adaptive Runge-Kutta-Chebyshev integrator, a
!< stabilized explicit RK scheme of 2nd order accuracy and variable stages.
!<
!< This particular implementation is a refactored and significantly modified
!< version of RKC compared to what was added to NCAR-LES by Kat Smith and
!< Kyle Niemeyer.
    use integrand, only: integrand_type, IK, RK
    implicit none
    private
    public :: rkc_type

    real, parameter :: UROUND = epsilon(1.0)

    type rkc_type
        private
        ! public variables
        real(RK), protected, public :: t, dt !< current time, timestep size
        class(integrand_type), pointer, public :: y !< user-supplied Integrand

        ! private variables with persistent values between calls to integrate() or step()
        real(RK) :: rtol, atol !< error tolerances
        real(RK) :: rho !< spectral radius
        integer(IK) :: s, s_max, nstep !< current and max stage count, num steps since last computing rho
        real(RK), allocatable(:) :: dydt !< current time derivative
        class(integrand_type) :: eigenv !< current eigenvector

        ! private temporary storage
        real(RK), allocatable(:) :: work1, work2 !< internal work arrays
        class(integrand_type) :: y_work !< internal work object

    contains
        procedure, pass(self), public :: initialize
        procedure, pass(self), public :: destroy
        procedure, pass(self), public :: integrate
    end type rkc_type

contains

    ! Helper function not a part of rkc_type -----------------------------------
    pure function amaxabs(a1, a2)
        !! Element-wise maximum of the absolute values of two arrays
        real(RK), intent(in) :: a1(:)
        real(RK), intent(in) :: a2(size(a1))
        real(RK) :: amaxabs(size(a1))
        integer :: i
        do i = 1, size(a1)
            amaxabs(i) = max(abs(a1(i)), abs(a2(i)))
        end do
    end function amaxabs

    !---------------------------------------------------------------------------
    ! RKC_TYPE METHODS
    !---------------------------------------------------------------------------
    subroutine destroy(self)
        !< Destroy integrator.
        class(rkc_type), intent(inout) :: self
        type(rkc_type) :: fresh

        self = fresh ! derived-type dynamic lhs reallocation
    end subroutine destroy

    subroutine initialize(self, y, rtol, atol)
        !! Initialize the RKC integrator working memory and RHS integrand
        class(rkc_type), intent(inout) :: self
            !! this integrator object
        class(integrand_type), intent(in) :: y
            !! the system of ODEs to be solved
        real(RK), intent(in), optional :: rtol
            !! relative tolerance value
        real(RK), intent(in), optional :: atol
            !! absolute tolerance value

        real(RK) :: hmin, hmax

        self%y => y
        self%eigenv = y ! lhs allocation by overloaded assignment
        self%y_work = y ! lhs allocation by overloaded assignment

        self%rtol = 1.0e-6; if (present(rtol)) self%rtol = rtol
        self%atol = 1.0e-10; if (present(atol)) self%atol = atol
        ! maximum number of RKC stages based on rtol (minimum is 2)
        self%s_max = max(2, nint(sqrt(0.1 * self%rtol / UROUND)))
        self%t = 0.0

        !> Get the spectral radius and eigenvector of the Jacobian
        hmin = self%atol        ! dummy tiny timestep
        hmax = 1.0 / self%atol  ! dummy huge timestep

        associate( &
            y_new => self%y_work, &
            dydt => self%dydt, &
            work => self%work1 &
        )
        dydt(:) = y%d_dt(self%t)
        self%eigenv = dydt
        rho = self%spec_rad(hmax) ! also updates eigenv!

        !> Estimate the timestep size
        self%dt = min(max(1.0 / self%rho, hmin), hmax)
        y_new = y + self%dt * dydt
        work(:) = y_new%d_dt(0.0)
        work(:) = (work - dydt) / (self%atol + self%rtol * abs(y%state()))
        err = self%dt * sqrt(sum(work**2) / real(y%size())) ! err is normalized RMS
        self%dt = min(max(0.1 * self%dt / sqrt(err), hmin), hmax) ! ensure hmin < h < hmax

        end associate

        !> Compute number of RKC stages for first step based on dt estimate
        self%s = 1 + nint(sqrt(1.54 * self%dt * self%rho + 1.0))
        if (self%s > self%s_max) then ! correct dt
            self%s = self%s_max
            self%dt = real(self%s**2 - 1) / (1.54 * self%rho)
        end if

    end subroutine initialize

    subroutine integrate(self, t_i, t_f)
        !! Integrate forward in time using adaptive Runge-Kutta-Chebyshev as an
        !! inner timestepper for stiff chemistry
        class(rkc_type), intent(inout) :: self
            !! this integrator object
        real(RK), intent(in) :: t_i
            !! The initial time
        real(RK), intent(in) :: t_f
            !! The final time

        real(RK) :: err_old, h_old, inv_size
        integer(IK) :: nstep, i
        real(RK) :: hmin, hmax, err, est, adapt, temp1, temp2

        self%t = t_i
        hmax = abs(t_f - t_i)                       ! maximum timestep size
        hmin = 10.0 * UROUND * max(abs(t_i), hmax)  ! minimum timestep size
        self%dt = min(max(self%dt, hmin), hmax)     ! ensure hmin < dt < hmax
        err_old = 0.0
        h_old = 0.0
        inv_size = 1.0 / real(self%y%size())

        associate( &
            dt => self%dt, &
            rho => self%rho, &
            y => self%y, &
            y_new => self%y_work, &
            dydt => self%dydt, &
            dydt_new => self%work1, &
            work => self%work2 &
        )
        do
            ! perform tentative time step
            y_new = self%internal_step()

            ! calculate F_np1 with tentative y_np1
            dydt_new(:) = y_new%d_dt(self%t)

            ! estimate error
            work(:) = 0.8 * (y - y_new) + 0.4 * dt * (dydt + dydt_new)
            work(:) = work / (atol + rtol * amaxabs(y_new%state(), y%state()))
            err = inv_size * sqrt(sum(work**2))

            ! If error too large, reject step and do not update self%t, etc.
            if (err >= 1.0) then
                dt = 0.8 * dt / (err**(1.0 / 3.0))
                rho = self%spec_rad(hmax)
                cycle
            end if

            ! Otherwise, accept the step and increment forward in time
            self%t = self%t + dt
            y = y_new
            dydt(:) = dydt_new
            nstep = nstep + 1

            ! if at t_f, exit loop (only loop termination point)
            if (self%t >= t_f) exit

            ! compute adaptive time step, stage count, based on error
            adapt = 10.0
            if (nstep > 1) then
                temp1 = 0.8 * dt * (err_old**(1.0 / 3.0))
                temp2 = h_old * (err**(2.0 / 3.0))
                if (temp1 < (adapt * temp2)) then
                    adapt = temp1 / temp2
                end if
            else ! at end of first time step, err_old and h_old do not exist
                temp2 = err**(1.0 / 3.0)
                if (0.8 < (adapt * temp2)) then
                    adapt = 0.8 / temp2
                end if
            end if
            err_old = err
            h_old = dt
            hmax = abs(t_f - self%t)
            hmin = 10.0 * UROUND * max(abs(self%t), hmax)
            dt = max(0.1, adapt) * dt ! adapt dt with a limit on shrink rate but not growth rate!?!?
            self%s = 1 + nint(sqrt(1.54 * dt * rho + 1.0))
            if (self%s > self%s_max) then
                self%s = self%s_max
                dt = real((self%s**2 - 1) / (1.54 * rho))
            end if
            dt = max(hmin, min(hmax, dt))

            ! re-estimate Jacobian spectral radius every 5 steps
            if (mod(nstep, 5) == 0) then
                rho = self%spec_rad(hmax)
            end if
        end do ! while loop
        end associate
    end subroutine integrate

    function spec_rad(hmax)
        !! Function to estimate upper bound of the spectral radius of stability
        real(RK), intent(in) :: hmax
            !! Maximum time step size
        real(RK) :: spec_rad
            !! function result

        integer(IK), parameter :: itmax = 50
        integer(IK) :: iter, ind
        real(RK) :: small, ny, y_rms, v_rms, dF_rms, dy_norm, sigma0, sigma1, tol, tol2

        small = 1.0 / hmax
        ny = real(self%y%size())

        associate( &
            y1d => self%y%state(), & ! state() calls return a pointer
            v1d => self%eigenv%state(), & ! state() calls return a pointer
            v => self%eigenv, & ! "inout" integrand object association
            F => self%dydt, &
            Fv => self%work1 &
        )

        y_rms = sqrt(sum(y1d**2))
        tol = sqrt(ny * UROUND * maxval(y1d**2)) ! floating-point roundoff near 0.0 for y_rms
        v_rms = sqrt(sum(v1d**2))
        tol2 = sqrt(ny * UROUND * maxval(v1d**2)) ! floating-point roundoff near 0.0 for v_rms

        ! dy_norm is the normalization for (v - y) and (Fv - F)
        if ((y_rms > tol) .and. (v_rms > tol2)) then
            dy_norm = y_rms * sqrt(UROUND)
            v = y1d + v1d * (dy_norm / v_rms)
        elseif (y_rms > tol) then
            print *, 'RKC TESTING COMMENT: hit v_rms < tol branch inside spec_rad'
            dy_norm = y_rms * sqrt(UROUND)
            v = y1d * (1.0 + sqrt(UROUND))
        elseif (v_rms > tol2) then
            print *, 'RKC TESTING COMMENT: hit y_rms < tol branch inside spec_rad'
            dy_norm = UROUND
            v = v1d * (dy_norm / v_rms)
        else
            print *, 'RKC TESTING COMMENT: hit v_rms, y_rms < tol branch inside spec_rad'
            dy_norm = UROUND
            Fv(:) = UROUND ! just using as temp array for integrand assignment
            v = Fv
        end if

        ! Now iterate using nonlinear power method
        Fv(:) = v%d_dt(self%t)
        dF_rms = sqrt(sum((Fv - F)**2))
        sigma1 = dF_rms / dy_norm

        do iter = 0, itmax - 1
            Fv(:) = v%d_dt(self%t)
            dF_rms = sqrt(sum((Fv - F)**2))
            sigma0 = sigma1
            sigma1 = dF_rms / dy_norm
            spec_rad = 1.2 * sigma1
            if (abs(sigma1 - sigma0) < (0.01 * max(sigma1, small))) then
                v = v1d - y1d
                self%nstep = 0 ! reset nstep
                return ! immediately terminate the subroutine
            end if

            tol = sqrt(ny * UROUND * maxval((Fv - F)**2)) ! FP tolerance for dF_rms
            if (dF_rms > tol) then
                v = y1d + ((Fv - F) * (dy_norm / dF_rms))
            else
                ind = mod(iter, v%size()) + 1
                v1d(ind) = -v1d(ind)
                print *, 'RKC TESTING COMMENT: hit dF_rms < tol branch inside spec_rad'
            end if
        end do
        end associate

        ! if you get to the end of the loop without hitting the alternate return ...
        print *, 'RKC NON-BLOCKING ERROR: spec_rad failed to converge!'

    end function spec_rad

    pure function internal_step() result(y_j)
        !! Function to take a single RKC integration step of variable stage count.
        real(RK), allocatable :: y_j(:)
            !! Final state at end of step


        ! variable RK stage coefficients, and related variables
        real(RK) :: w0, temp1, temp2, arg, w1, b_jm1, b_jm2, mu_t
        real(RK) :: c_jm2, c_jm1, zjm1, zjm2, dzjm1, dzjm2, d2zjm1, d2zjm2
        real(RK) :: zj, dzj, d2zj, b_j, gamma_t, nu, mu, c_j

        ! loop index
        integer :: j

        allocate(y_j(self%y%size()))
        associate( &
            dt => self%dt, &
            y_0 => self%y%state(), &
            F_0 => self%dydt, &
            y_jm1 => self%y_work, &
            y_jm2 => self%work1, &
        )

        w0 = 1.0 + 2.0 / (13.0 * real(self%s**2))
        temp1 = w0**2 - 1.0
        temp2 = sqrt(temp1)
        arg = real(self%s) * log(w0 + temp2)
        w1 = sinh(arg) * temp1 / (cosh(arg) * real(self%s) * temp2 - w0 * sinh(arg))

        b_jm1 = 0.25 / w0**2
        b_jm2 = b_jm1

        ! calculate y_1
        mu_t = w1 * b_jm1
        y_jm2(:) = y_0
        y_jm1 = y_0 + (mu_t * dt * F_0)

        c_jm2 = 0.0
        c_jm1 = mu_t
        zjm1 = w0
        zjm2 = 1.0
        dzjm1 = 1.0
        dzjm2 = 0.0
        d2zjm1 = 0.0
        d2zjm2 = 0.0

        do j = 2, self%s
            zj = 2.0 * w0 * zjm1 - zjm2
            dzj = 2.0 * w0 * dzjm1 - dzjm2 + 2.0 * zjm1
            d2zj = 2.0 * w0 * d2zjm1 - d2zjm2 + 4.0 * dzjm1
            b_j = d2zj / dzj**2
            gamma_t = 1.0 - (zjm1 * b_jm1)

            nu = -b_j / b_jm2
            mu = 2.0 * b_j * w0 / b_jm1
            mu_t = mu * w1 / w0

            y_j(:) = y_jm1%d_dt(self%t) ! calculate derivative

            y_j(:) = (1.0 - mu - nu) * y_0 + (mu * y_jm1) + (nu * y_jm2) &
                     + dt * mu_t * (y_j - (gamma_t * F_0))
            c_j = (mu * c_jm1) + (nu * c_jm2) + mu_t * (1.0 - gamma_t)

            y_jm2(:) = y_jm1
            y_jm1 = y_j

            c_jm2 = c_jm1
            c_jm1 = c_j
            b_jm2 = b_jm1
            b_jm1 = b_j
            zjm2 = zjm1
            zjm1 = zj
            dzjm2 = dzjm1
            dzjm1 = dzj
            d2zjm2 = d2zjm1
            d2zjm1 = d2zj
        end do
        end associate

    end function internal_step

end module rkc_integrator
