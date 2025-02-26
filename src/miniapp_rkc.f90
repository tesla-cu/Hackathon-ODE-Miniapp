module miniapp_rkc
    !! The error-controlled adaptive Runge-Kutta-Chebyshev integrator, a
    !! stabilized explicit RK scheme of 2nd order accuracy and variable stages.
    !!
    !! This particular implementation is a refactored and somewhat modified
    !! version of RKC compared to what was added to NCAR-LES by Kat Smith and
    !! Kyle Niemeyer. For a version of RKC faithful to the NCAR-LES implementation,
    !! see the `ncarles_rkc` module.
    implicit none
    private
    public :: initialize_rkc, rkc_integrate, rkc_inplace_step

    real, parameter :: UROUND = epsilon(1.0)
    real :: rel_tol, abs_tol
    integer :: s_max

    interface
        pure subroutine time_derivative(t, y, ydot, p)
            real, intent(in) :: t
            real, intent(in) :: y(0:)
            real, intent(inout) :: ydot(0:)
            real, intent(in), optional :: p(0:)
        end subroutine time_derivative
    end interface

contains

    subroutine initialize_rkc(rtol, atol)
        !! Initialize the RKC integrator's working memory and RHS function pointer.
        real, intent(in), optional :: rtol
            !! relative tolerance value
        real, intent(in), optional :: atol
            !! absolute tolerance value

        rel_tol = 1.0e-6; if (present(rtol)) rel_tol = rtol
        abs_tol = 1.0e-10; if (present(atol)) abs_tol = atol
        ! maximum number of RKC stages based on rel_tol
        s_max = max(2, nint(sqrt(0.1 * rel_tol / UROUND)))

    end subroutine initialize_rkc

    subroutine rkc_inplace_step(rhs, t, dt, y, p)
        procedure(time_derivative) :: rhs
        real, intent(in) :: t, dt
        real, intent(inout) :: y(:)
        real, intent(in), optional :: p(:)
        real :: ydot(size(y))  ! initial derivative at time `t`
        call rhs(t, y, ydot, p)
        y(:) = rkc_step(rhs, t, dt, s_max, y, ydot, p)
    end subroutine rkc_inplace_step

    subroutine rkc_integrate(rhs, t_i, t_f, y, p)
        !! Integrate forward in time using adaptive Runge-Kutta-Chebyshev as an
        !! inner timestepper for stiff chemistry
        procedure(time_derivative) :: rhs
        real, intent(in) :: t_i
            !! The initial time.
        real, intent(in) :: t_f
            !! The desired final time.
        real, intent(inout) :: y(1:)
            !! Solution vector at `t_i` on input, at `t_f` on output
        real, intent(in), optional :: p(:)
            !! optional parameters to pass on to rhs subroutine

        ! work arrays
        real, dimension(size(y)) :: y_end, ydot, vtemp1, vtemp2, eigenv
        ! internal work arrays: not all may be necessary, haven't figured out minimum
        ! needed working memory. `eigenv` used to be stored in work(4:)

        real :: err_old, h_old, h_n, rho
        ! these variables used to be stored in work(0:3)
        integer :: ny, nstep, s, i
        real :: t_rkc, hmax, hmin, err, est, adapt, temp1, temp2

        ny = size(y)
        t_rkc = t_i
        nstep = 0

        ! Initialize integration limiters
        hmax = abs(t_f - t_i) ! maximum timestep size
        hmin = 10.0 * UROUND * max(abs(t_i), hmax) ! minimum timestep size
        call rhs(t_rkc, y, ydot, p) ! calculate RHS for initial y input
        eigenv = ydot ! initial estimate of eigenvector
        rho = rkc_spec_rad(rhs, t_rkc, hmax, y, ydot, eigenv, vtemp2, p)
        err_old = 0.0
        h_old = 0.0

        ! --> Estimate the timestep size, h_n
        h_n = hmax
        if (rho * h_n > 1.0) h_n = 1.0 / rho
        h_n = max(h_n, hmin)

        vtemp1 = y + h_n * ydot
        call rhs(t_rkc, vtemp1, vtemp2, p)

        vtemp1 = (vtemp2 - ydot) / (abs_tol + rel_tol * abs(y))
        err = h_n * sqrt(sum(vtemp1**2) / real(ny))
        h_n = min(max(0.1 * h_n / sqrt(err), hmin), hmax)
        ! <-- end of estimate h_n

        ! compute number of RKC stages for first step based on h_n estimate
        s = 1 + nint(sqrt(1.54 * h_n * rho + 1.0))
        if (s > s_max) then ! correct h_n
            s = s_max
            h_n = real((s**2 - 1) / (1.54 * rho))
        end if

        ! INTEGRATE TO END TIME
        do while (t_rkc < t_f)
            ! perform tentative time step
            y_end = rkc_step(rhs, t_rkc, h_n, s, y, ydot, p)

            ! calculate F_np1 with tenative y_np1
            call rhs(t_rkc, y_end, vtemp1, p)

            ! estimate error
            err = 0.0
            do i = 1, ny
                est = 0.8 * (y(i) - y_end(i)) + 0.4 * h_n * (ydot(i) + vtemp1(i))
                est = est / (abs_tol + rel_tol * max(abs(y_end(i)), abs(y(i))))
                err = err + est**2
            end do
            err = sqrt(err / real(ny))

            ! If error too large, reject step and do not update t_rkc, etc.
            if (err >= 1.0) then
                h_n = 0.8 * h_n / (err**(1.0 / 3.0))
                rho = rkc_spec_rad(rhs, t_rkc, hmax, y, ydot, eigenv, vtemp2, p)
                cycle
            end if

            ! Otherwise, accept the step and increment forward in time
            t_rkc = t_rkc + h_n
            y = y_end
            ydot = vtemp1
            nstep = nstep + 1

            ! if at t_f, exit loop immediately
            ! if (t_rkc >= t_f) exit

            ! compute adaptive time step, stage count, based on error
            adapt = 10.0
            if (nstep > 1) then
                temp1 = 0.8 * h_n * (err_old**(1.0 / 3.0))
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
            h_old = h_n
            hmax = abs(t_f - t_rkc)
            hmin = 10.0 * UROUND * max(abs(t_rkc), hmax)
            h_n = max(0.1, adapt) * h_n
            ! grow/shrink h_n with a limit on shrink rate but not growth rate!?!?
            ! Normally Colin would allow any shrink rate, but limit growth rate
            s = 1 + nint(sqrt(1.54 * h_n * rho + 1.0))
            if (s > s_max) then
                s = s_max
                h_n = real((s**2 - 1) / (1.54 * rho))
            end if
            h_n = max(hmin, min(hmax, h_n)) ! bound h_n by min/max values

            ! re-estimate Jacobian spectral radius every 25 steps
            if (mod(nstep, 25) == 0) then
                rho = rkc_spec_rad(rhs, t_rkc, hmax, y, ydot, eigenv, vtemp2, p)
            end if

        end do ! while loop

    end subroutine rkc_integrate

    function rkc_spec_rad(rhs, t_rkc, hmax, y, F, v, Fv, p)
        !! Function to estimate upper bound of the spectral radius of stability
        procedure(time_derivative) :: rhs
        real, intent(in) :: t_rkc
            !! The current time.
            !! NOTE: unused dummy argument! Kept for backwards compatibility.
            !! (Generally required as an argument to ODE function integrators)
        real, intent(in) :: hmax
            !! Maximum time step size
        real, intent(in) :: y(1:)
            !! Current solution
        real, intent(in) :: F(1:size(y))
            !! Time derivative of solution, dy/dt = F(y)
        real, intent(inout) :: v(1:size(y))
            !! Estimate of ODE system eigenvalues
        real, intent(out) :: Fv(1:size(y))
            !! Time derivative of eigenalues, dv/dt = F(v)
        real, intent(in), optional :: p(:)
            !! optional parameters to pass on to rhs subroutine
        real :: rkc_spec_rad
            !! function result

        integer, parameter :: itmax = 50
        integer :: iter, ind
        real :: small, y_rms, v_rms, dF_rms, dy_norm, sigma0, sigma1

        small = 1.0 / hmax

        y_rms = sqrt(sum(y**2))
        v_rms = sqrt(sum(v**2))
        if ((y_rms > 0.0) .and. (v_rms > 0.0)) then
            dy_norm = y_rms * sqrt(UROUND)
            v = y + v * (dy_norm / v_rms)
        elseif (y_rms > 0.0) then
            dy_norm = y_rms * sqrt(UROUND)
            v = y * (1.0 + sqrt(UROUND))
        elseif (v_rms > 0.0) then
            dy_norm = UROUND
            v = v * (dy_norm / v_rms)
        else
            dy_norm = UROUND
            v = UROUND
        end if

        ! now iterate using nonlinear power method
        sigma1 = 0.0
        do iter = 1, itmax
            call rhs(t_rkc, v, Fv, p)
            dF_rms = sqrt(sum((Fv - F)**2))
            sigma0 = sigma1
            sigma1 = dF_rms / dy_norm
            rkc_spec_rad = 1.2 * sigma1
            if ((iter >= 2) .and. (abs(sigma1 - sigma0) <= (max(sigma1, small) * 0.01))) then
                v = v - y
                return ! immediately terminate the subroutine
            end if

            if (dF_rms > 0.0) then
                v = y + ((Fv - F) * (dy_norm / dF_rms))
            else
                ind = mod(iter, size(y))
                v(ind) = y(ind) - (v(ind) - y(ind))
                ! this part is from original RKC code, but doesn't make sense to Colin.
                ! Based on RKC paper, what is wanted is v(ind) = -v(ind). If we
                ! had infinite precision, then this would set v(ind) = 2*y(ind) - v(ind).
                print *, 'RKC TESTING: reached dF_rms == 0 branch inside rkc_spec_rad'
            end if
        end do

        ! if you get to the end of the loop without hitting the alternate return ...
        print *, 'RKC ERROR: rkc_spec_rad failed to converge!'
        print *, 't_rkc, 1/hmax, rho = ', t_rkc, 1/hmax, rkc_spec_rad
        print *, 'lbdound = ', lbound(y), 'y = ', y
        print *, 'lbdound = ', lbound(v), 'v = ', v
        print *, 'lbdound = ', lbound(F), 'F = ', F
        print *, 'lbdound = ', lbound(Fv), 'Fv = ', Fv
        error stop 'ERROR STOP'

    end function rkc_spec_rad

    function rkc_step(rhs, t_rkc, h, s, y_0, F_0, p) result(y_j)
        !! Function to take a single RKC integration step of variable stage count.
        procedure(time_derivative) :: rhs
        real, intent(in) :: t_rkc
            !! The current time. Here for generic integrator compatibility.
        real, intent(in) :: h
            !! The time step size
        integer, intent(in) :: s
            !! number of stages to compute
        real, intent(in) :: y_0(1:)
            !! The current solution
        real, intent(in) :: F_0(size(y_0))
            !! The time derivative of current solution, dy/dt = F(y)
        real, intent(in), optional :: p(:)
            !! optional parameters to pass on to rhs subroutine
        real :: y_j(size(y_0))
            !! The solution at the end of the step

        ! internal work memory
        real, dimension(size(y_0)) :: y_jm1, y_jm2

        ! variable RK stage coefficients, and related variables
        real :: w0, temp1, temp2, arg, w1, b_jm1, b_jm2, mu_t
        real :: c_jm2, c_jm1, zjm1, zjm2, dzjm1, dzjm2, d2zjm1, d2zjm2
        real :: zj, dzj, d2zj, b_j, gamma_t, nu, mu, c_j

        ! loop index
        integer :: j

        w0 = 1.0 + 2.0 / (13.0 * real(s**2))
        temp1 = w0**2 - 1.0
        temp2 = sqrt(temp1)
        arg = real(s) * log(w0 + temp2)
        w1 = sinh(arg) * temp1 / (cosh(arg) * real(s) * temp2 - w0 * sinh(arg))

        b_jm1 = 0.25 / w0**2
        b_jm2 = b_jm1

        ! calculate y_1
        mu_t = w1 * b_jm1
        y_jm2(:) = y_0
        y_jm1(:) = y_0 + (mu_t * h * F_0)

        c_jm2 = 0.0
        c_jm1 = mu_t
        zjm1 = w0
        zjm2 = 1.0
        dzjm1 = 1.0
        dzjm2 = 0.0
        d2zjm1 = 0.0
        d2zjm2 = 0.0

        do j = 2, s
            zj = 2.0 * w0 * zjm1 - zjm2
            dzj = 2.0 * w0 * dzjm1 - dzjm2 + 2.0 * zjm1
            d2zj = 2.0 * w0 * d2zjm1 - d2zjm2 + 4.0 * dzjm1
            b_j = d2zj / dzj**2
            gamma_t = 1.0 - (zjm1 * b_jm1)

            nu = -b_j / b_jm2
            mu = 2.0 * b_j * w0 / b_jm1
            mu_t = mu * w1 / w0

            ! calculate derivative
            call rhs(t_rkc, y_jm1, y_j, p)

            y_j(:) = (1.0 - mu - nu) * y_0 + (mu * y_jm1) + (nu * y_jm2) &
                     + h * mu_t * (y_j - (gamma_t * F_0))
            c_j = (mu * c_jm1) + (nu * c_jm2) + mu_t * (1.0 - gamma_t)

            y_jm2(:) = y_jm1
            y_jm1(:) = y_j

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

    end function rkc_step

end module miniapp_rkc
