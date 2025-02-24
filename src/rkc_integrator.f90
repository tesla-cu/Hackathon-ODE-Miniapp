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
        real(RK), public :: t, dt !< current time, timestep size
        class(integrand_type), pointer, public :: y => null() !< user-supplied Integrand

        ! private variables with persistent values between calls to integrate() or step()
        real(RK) :: rtol, atol !< error tolerances
        real(RK) :: rho !< spectral radius
        integer(IK) :: s, s_max !< current and max stage count, num steps since last computing rho
        real(RK), allocatable :: dydt(:) !< current time derivative
        class(integrand_type), allocatable :: eigenv !< current eigenvector

        ! private temporary storage
        real(RK), allocatable :: work1(:), work2(:) !< internal work arrays
        class(integrand_type), allocatable :: y_work !< internal work object

    contains
        procedure, pass(self), public :: initialize
        procedure, pass(self), public :: destroy
        procedure, pass(self), public :: integrate
        procedure, pass(self), public :: spec_rad
        procedure, pass(self), public :: internal_step
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
        !! Destroy integrator.
        !!
        !! Using intrinsic assignment here overwrites self with a new object
        !! that has not allocated any of its data, thus implicitly deallocating
        !! all data associated with the old object. Hopefully. That's how FOODIE
        !! does it, but I'm definitely going to look at this in a profiler.
        class(rkc_type), intent(inout) :: self
        deallocate(self%eigenv, self%y_work, self%dydt, self%work1, self%work2)
        nullify(self%y)
    end subroutine destroy

    subroutine initialize(self, y, rtol, atol)
        !! Initialize the RKC integrator working memory and RHS integrand
        class(rkc_type), intent(inout) :: self
            !! this integrator object
        class(integrand_type), target, intent(in) :: y
            !! the system of ODEs to be solved
        real(RK), intent(in), optional :: rtol
            !! relative tolerance value
        real(RK), intent(in), optional :: atol
            !! absolute tolerance value

        real(RK) :: hmin, hmax, err

        self%y => y ! point to this, don't make a copy
        self%rtol = 1.0e-6; if (present(rtol)) self%rtol = rtol
        self%atol = 1.0e-10; if (present(atol)) self%atol = atol

        allocate (self%dydt(y%size()), self%work1(y%size()), self%work2(y%size()))
        allocate (self%eigenv, self%y_work, mold=y)
        self%eigenv = y
        self%y_work = y

        !> maximum number of RKC stages based on rtol (minimum is 2)
        self%s_max = max(2, nint(sqrt(0.1 * self%rtol / UROUND)))

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
        real(RK) :: hmin, hmax, err, adapt, temp1, temp2
        integer :: nstep, ncycles

        self%t = t_i
        hmax = abs(t_f - t_i)                       ! maximum timestep size
        hmin = 10.0 * UROUND * max(abs(t_i), hmax)  ! minimum timestep size

        associate( &
            dt => self%dt, &
            rho => self%rho, &
            y => self%y, &
            y_new => self%y_work, &
            dydt => self%dydt, &
            dydt_new => self%work1, &
            work => self%work2 &
        )

        !> initialize rho and eigenv
        dydt(:) = y%d_dt(self%t)
        self%eigenv = dydt
        rho = self%spec_rad(hmax) ! also updates eigenv!

        !> Estimate the timestep size
        dt = min(max(1.0 / rho, hmin), hmax)
        y_new = y + dt * dydt
        work(:) = y_new%d_dt(self%t)
        work(:) = (work - dydt) / (self%atol + self%rtol * abs(y%state()))
        err = dt * sqrt(sum(work**2) / real(y%size())) ! err is normalized RMS
        dt = min(max(0.1 * dt / sqrt(err), hmin), hmax) ! ensure hmin < h < hmax

        !> Compute number of RKC stages for first step based on dt estimate
        self%s = 1 + nint(sqrt(1.54 * dt * rho + 1.0))
        if (self%s > self%s_max) then ! correct dt
            self%s = self%s_max
            dt = real(self%s**2 - 1) / (1.54 * rho)
        end if

        ! print *, 'RKC initial dt = ', dt
        err_old = 0.0
        h_old = 0.0
        nstep = 0
        ! ncycles = 0
        do
            !> perform tentative time step
            y_new = self%internal_step()
            ! ncycles = ncycles + 1

            !> calculate F_np1 with tentative y_np1
            dydt_new(:) = y_new%d_dt(self%t)

            !> estimate error
            work(:) = 0.8 * (y - y_new) + 0.4 * dt * (dydt + dydt_new)
            work(:) = work / (self%atol + self%rtol * amaxabs(y_new%state(), y%state()))
            err = sqrt(sum(work**2) / real(y%size()))

            !> If error too large, reject step and do not update self%t, etc.
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
            else ! at end of first time step, err_old and h_old == 0.0
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
            if (mod(nstep, 25) == 0) then
                rho = self%spec_rad(hmax)
            end if
        end do ! while loop
        ! print *, 'RKC final dt, steps = ', dt, nstep, ncycles
        end associate
    end subroutine integrate

    function spec_rad(self, hmax)
        !! Function to estimate upper bound of the spectral radius of stability
        class(rkc_type), intent(inout) :: self
            !! this integrator object
        real(RK), intent(in) :: hmax
            !! Maximum time step size
        real(RK) :: spec_rad
            !! function result

        integer(IK), parameter :: itmax = 50
        integer(IK) :: iter, ind
        real(RK) :: small, ny, y_rms, v_rms, dF_rms, dy_norm, sigma0, sigma1, tol, tol2
        real(RK), contiguous, pointer :: vstate(:), ystate(:)

        small = 1.0 / hmax
        ny = real(self%y%size())

        associate( &
            y => self%y, &
            v => self%eigenv, &
            F => self%dydt, &
            Fv => self%work1 &
        )
        ystate => y%state()
        y_rms = sqrt(sum(ystate**2))
        tol = sqrt(ny * UROUND * maxval(ystate**2)) ! floating-point roundoff near 0.0 for y_rms

        vstate => v%state()
        v_rms = sqrt(sum(vstate**2))
        tol2 = sqrt(ny * UROUND * maxval(vstate**2)) ! floating-point roundoff near 0.0 for v_rms

        ! dy_norm is the normalization for (v - y) and (Fv - F)
        if ((y_rms > tol) .and. (v_rms > tol2)) then
            dy_norm = y_rms * sqrt(UROUND)
            v = y + v * (dy_norm / v_rms)
        elseif (y_rms > tol) then
            print *, 'RKC TESTING COMMENT: hit v_rms < tol branch inside spec_rad'
            dy_norm = y_rms * sqrt(UROUND)
            v = y * (1.0 + sqrt(UROUND))
        elseif (v_rms > tol2) then
            print *, 'y_rms = ', y_rms
            print *, 'tol = ', tol
            print *, 'ystate = ', ystate
            error stop 'RKC TESTING: hit y_rms < tol branch inside spec_rad'
            dy_norm = UROUND
            v = v * (dy_norm / v_rms)
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
        v = y + ((Fv - F) * (dy_norm / dF_rms))

        do iter = 0, itmax - 1
            Fv(:) = v%d_dt(self%t)
            dF_rms = sqrt(sum((Fv - F)**2))
            sigma0 = sigma1
            sigma1 = dF_rms / dy_norm
            spec_rad = 1.2 * sigma1
            if (abs(sigma1 - sigma0) < (0.01 * max(sigma1, small))) then
                v = v - y
                return ! immediately terminate the subroutine
            end if

            tol = sqrt(ny * UROUND * maxval((Fv - F)**2)) ! FP tolerance for dF_rms
            if (dF_rms > tol) then
                v = y + ((Fv - F) * (dy_norm / dF_rms))
            else
                ind = mod(iter, v%size()) + 1
                vstate(ind) = -vstate(ind)
                print *, 'dF_rms = ', dF_rms
                print *, 'tol = ', tol
                print *, 'vstate = ', vstate
                print *, 'ystate = ', ystate
                error stop 'RKC TESTING: hit dF_rms < tol branch inside spec_rad'
            end if
        end do
        end associate

        ! if you get to the end of the loop without hitting the alternate return ...
        print *, 'RKC NON-BLOCKING ERROR: spec_rad failed to converge!'

    end function spec_rad

    function internal_step(self) result(y_j)
        !! Function to take a single RKC integration step of variable stage count.
        class(rkc_type), intent(inout) :: self
            !! this integrator object
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
            y_jm2 => self%work1 &
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

            y_jm2(:) = y_jm1%state()
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
