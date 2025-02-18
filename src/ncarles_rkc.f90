module ncarles_rkc
    !! The error-controlled adaptive Runge-Kutta-Chebyshev integrator, a
    !! stabilized explicit RK scheme of 2nd order accuracy and variable stages.
    !!
    !! This particular implementation is what was added to NCAR-LES by Kat Smith and
    !! Kyle Niemeyer, with some very minor alterations to accomodate being placed
    !! into this miniapp.
    !! Niemeyer, in turn, derived the code from the original RKC, distributed by
    !! Netlib, https://www.netlib.org/ode/rkc.f
    implicit none
    private
    public :: initialize_rkc, rkc_integrate, rkc_inplace_step

    real, parameter :: UROUND = 2.22e-16
    real :: rel_tol, abs_tol
    integer :: m_max

    abstract interface
        subroutine time_derivative(t, y, ydot, p)
            real, intent(in) :: t, y(:)
            real, intent(inout) :: ydot(:)
            real, intent(in), optional :: p(:)
        end subroutine time_derivative
    end interface
    procedure(time_derivative), pointer :: rhs

contains

    subroutine initialize_rkc(dydt, rtol, atol)
        !! Initialize the RKC integrator's working memory and RHS function pointer.
        procedure(time_derivative) :: dydt
            !! User-supplied RHS term of ODE, dy/dt = RHS, with signature
            !! `subroutine dydt(t, y, ydot, p)`
        real, intent(in), optional :: rtol
            !! relative tolerance value
        real, intent(in), optional :: atol
            !! absolute tolerance value

        rhs => dydt

        rel_tol = 1.0e-6; if (present(rtol)) rel_tol = rtol
        abs_tol = 1.0e-10; if (present(atol)) abs_tol = atol
        ! maximum number of RKC stages based on rel_tol
        m_max = max(2, nint(sqrt(0.1 * rel_tol / UROUND)))
    end subroutine initialize_rkc

    subroutine rkc_integrate(ti, tf, y, p)
        real, intent(in) :: ti, tf
        real, intent(inout) :: y(:)
        real, intent(in), optional :: p(:)
        y(:) = rkc_driver(ti, tf, size(y), y, p)
    end subroutine rkc_integrate

    subroutine rkc_inplace_step(t, dt, y, p)
        real, intent(in) :: t, dt
        real, intent(inout) :: y(:)
        real, intent(in), optional :: p(:)
        real :: ydot(size(y))  ! initial derivative at time `t`
        call rhs(t, y, ydot, p)
        y(:) = rkc_step(t, dt, size(y), y, ydot, m_max, p)
    end subroutine rkc_inplace_step

    function rkc_driver(t_0, t_end, nscl, y_0, p) result(y_end)
        ! Driver function for RKC integrator.
        !
        ! t_0    the starting time.
        ! t_end  the desired end time.
        ! y_0    dependent variable array
        implicit none
        real, intent(in) :: t_0
        real, intent(in) :: t_end
        integer, intent(in) :: nscl
        real, intent(in), dimension(0:nscl-1) :: y_0
        real, intent(in), optional :: p(:)
        real, dimension(0:nscl-1) :: y_end

        real, dimension(0:nscl+4) :: work
        real, dimension(0:nscl-1) :: y_n, F_n, temp_arr, temp_arr2
        integer nstep, i, m
        real hmax, hmin, err, est
        real fac, temp1, temp2, t_rkc

        hmax = abs(t_end - t_0)
        hmin = 10.0 * UROUND * max(abs(t_0), hmax)

        ! initialize variables for start of loop
        t_rkc = t_0
        y_n = y_0
        work(:) = 0.0
        nstep = 0

        ! calculate F_n for y_0, and store as eigenvector estimate
        call rhs(t_rkc, y_n, F_n, p)
        work(4:) = F_n

        do while (t_rkc < t_end)
            ! estimate Jacobian spectral radius only if 25 steps passed
            if (mod(nstep, 25) == 0) then
                work(3) = rkc_spec_rad(t_rkc, hmax, y_n, F_n, work(4:), temp_arr2, p)
            end if

            ! first step, estimate step size
            if (work(2) < UROUND) then
                work(2) = hmax
                if ((work(3) * work(2)) > 1.0) then
                    work(2) = 1.0 / work(3)
                end if
                work(2) = max(work(2), hmin)

                temp_arr = y_n + (work(2) * F_n)
                call rhs(t_rkc, temp_arr, temp_arr2, p)

                temp_arr = (temp_arr2 - F_n) / (abs_tol + rel_tol * abs(y_n))
                err = work(2) * sqrt(sum(temp_arr**2) / real(nscl))

                if ((0.1 * work(2)) < (hmax * sqrt(err))) then
                    work(2) = max((0.1 * work(2)) / sqrt(err), hmin)
                else
                    work(2) = hmax
                end if
            end if

            ! check if last step
            if ((1.1 * work(2)) >= abs(t_end - t_rkc)) then
                work(2) = abs(t_end - t_rkc)
            end if

            ! calculate number of steps
            m = 1 + nint(sqrt(1.54 * work(2) * work(3) + 1.0))

            if (m > m_max) then
                m = m_max
                work(2) = real((m * m - 1) / (1.54 * work(3)))
            end if

            hmin = 10.0 * UROUND * max(abs(t_rkc), abs(t_rkc + work(2)))

            ! perform tentative time step
            y_end = rkc_step(t_rkc, work(2), nscl, y_n, F_n, m, p)

            ! calculate F_np1 with tenative y_np1
            call rhs(t_rkc, y_end, temp_arr, p)

            ! estimate error
            err = 0.0
            do i = 0, nscl - 1
                est = 0.0
                est = 0.8 * (y_n(i) - y_end(i)) + 0.4 * work(2) * (F_n(i) + temp_arr(i))
                est = est / (abs_tol + rel_tol * max(abs(y_end(i)), abs(y_n(i))))
                err = err + est * est
            end do
            err = sqrt(err / real(nscl))

            if (err > 1.0) then
                ! error too large, step is rejected.
                !select smaller step size
                work(2) = 0.8 * work(2) / (err**(1.0 / 3.0))
                ! reevaluate spectral radius
                work(3) = rkc_spec_rad(t_rkc, hmax, y_n, F_n, work(4:), temp_arr2, p)

            else
                ! step accepted
                t_rkc = t_rkc + work(2)
                nstep = nstep + 1

                fac = 10.0
                temp1 = 0.0
                temp2 = 0.0
                if (work(1) < UROUND) then
                    temp2 = err**(1.0 / 3.0)
                    if (0.8 < (fac * temp2)) then
                        fac = 0.8 / temp2
                    end if
                else
                    temp1 = 0.8 * work(2) * (work(0)**(1.0 / 3.0))
                    temp2 = work(1) * (err**(2.0 / 3.0))
                    if (temp1 < (fac * temp2)) then
                        fac = temp1 / temp2
                    end if
                end if

                ! set "old" values to those for current time step
                work(0) = err
                work(1) = work(2)

                do i = 0, nscl - 1
                    y_n(i) = y_end(i)
                    F_n(i) = temp_arr(i)
                end do

                ! store next time step
                work(2) = work(2) * max(0.1, fac)
                work(2) = max(hmin, min(hmax, work(2)))

            end if
        end do

    end function rkc_driver

    real function rkc_spec_rad(t_rkc, hmax, yLocal, F, v, Fv, p)
        ! Function to estimate spectral radius.
        !
        ! t_rkc    the time.
        ! hmax     Max time step size.
        ! yLocal   Array of dependent variable.
        ! F        Derivative evaluated at current state
        ! v
        ! Fv
        implicit none
        real, intent(in) :: t_rkc
        real, intent(in) :: hmax
        real, intent(in), dimension(0:) :: yLocal
        real, intent(inout), dimension(0:) :: v, Fv, F
        real, intent(in), optional :: p(:)

        integer nscl, itmax, i, iter, ind
        real small, nrm1, nrm2, dynrm, sigma

        nscl = size(yLocal)
        itmax = 50
        small = 1.0 / hmax
        nrm1 = 0.0
        nrm2 = 0.0
        sigma = 0.0

        do i = 0, nscl - 2
            nrm1 = nrm1 + yLocal(i) * yLocal(i)
            nrm2 = nrm2 + v(i) * v(i)
        end do
        nrm1 = sqrt(nrm1)
        nrm2 = sqrt(nrm2)

        if ((nrm1 > 0.0) .and. (nrm2 > 0.0)) then
            dynrm = nrm1 * sqrt(UROUND)
            do i = 0, nscl - 1
                v(i) = yLocal(i) + v(i) * (dynrm / nrm2)
            end do
        elseif (nrm1 > 0.0) then
            dynrm = nrm1 * sqrt(UROUND)
            do i = 0, nscl - 1
                v(i) = yLocal(i) * (1.0 + sqrt(UROUND))
            end do
        elseif (nrm2 > 0.0) then
            dynrm = UROUND
            do i = 0, nscl - 1
                v(i) = v(i) * (dynrm / nrm2)
            end do
        else
            dynrm = UROUND
            do i = 0, nscl - 1
                v(i) = UROUND
            end do
        end if

        ! now iterate using nonlinear power method
        sigma = 0.0
        do iter = 1, itmax
            call rhs(t_rkc, v, Fv, p)

            nrm1 = 0.0
            do i = 0, nscl - 1
                nrm1 = nrm1 + ((Fv(i) - F(i)) * (Fv(i) - F(i)))
            end do
            nrm1 = sqrt(nrm1)
            nrm2 = sigma
            sigma = nrm1 / dynrm
            if ((iter >= 2) .and. (abs(sigma - nrm2) <= (max(sigma, small) * 0.01))) then
                do i = 0, nscl - 1
                    v(i) = v(i) - yLocal(i)
                end do
                rkc_spec_rad = 1.2 * sigma
            end if

            if (nrm1 > 0.0) then
                do i = 0, nscl - 1
                    v(i) = yLocal(i) + ((Fv(i) - F(i)) * (dynrm / nrm1))
                end do
            else
                ind = mod(iter, int(nscl - 1))
                v(ind) = yLocal(ind) - (v(ind) - yLocal(ind))
            end if
        end do

        rkc_spec_rad = 1.2 * sigma

    end function rkc_spec_rad

    function rkc_step(t_rkc, h, nscl, y_0, F_0, s, p)
        ! Function to take a single RKC integration step
        !
        ! t_rkc    the starting time.
        ! h        Time-step size.
        ! y_0      Initial conditions.
        ! F_0      Derivative function at initial conditions.
        ! s        number of steps.
        ! rkc_step Integrated variables
        implicit none
        real, intent(in) :: t_rkc
        real, intent(in) :: h
        integer, intent(in) :: nscl
        real, intent(in), dimension(0:nscl-1) :: y_0, F_0
        integer, intent(in) :: s
        real, intent(in), optional :: p(:)
        real, dimension(0:nscl-1) :: rkc_step

        real, dimension(0:nscl-1) :: y_j
        real, dimension(0:nscl-1) :: y_jm1, y_jm2

        real w0, temp1, temp2, arg, w1, b_jm1, b_jm2, mu_t
        real c_jm2, c_jm1, zjm1, zjm2, dzjm1, dzjm2, d2zjm1, d2zjm2
        real zj, dzj, d2zj, b_j, gamma_t, nu, mu, c_j
        integer i, j

        w0 = 1.0 + 2.0 / (13.0 * real(s * s))
        temp1 = (w0 * w0) - 1.0
        temp2 = sqrt(temp1)
        arg = real(s) * log(w0 + temp2)
        w1 = sinh(arg) * temp1 / (cosh(arg) * real(s) * temp2 - w0 * sinh(arg))

        b_jm1 = 1.0 / (4.0 * (w0 * w0))
        b_jm2 = b_jm1

        ! calculate y_1
        mu_t = w1 * b_jm1
        do i = 0, nscl - 1
            y_jm2(i) = y_0(i)
            y_jm1(i) = y_0(i) + (mu_t * h * F_0(i))
        end do

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
            b_j = d2zj / (dzj * dzj)
            gamma_t = 1.0 - (zjm1 * b_jm1)

            nu = -b_j / b_jm2
            mu = 2.0 * b_j * w0 / b_jm1
            mu_t = mu * w1 / w0

            ! calculate derivative, use y array for temporary storage
            call rhs(t_rkc, y_jm1, y_j, p)

            do i = 0, nscl - 1
                y_j(i) = (1.0 - mu - nu) * y_0(i) + (mu * y_jm1(i)) + (nu * y_jm2(i)) &
                        + h * mu_t * (y_j(i) - (gamma_t * F_0(i)))
            end do
            c_j = (mu * c_jm1) + (nu * c_jm2) + mu_t * (1.0 - gamma_t)

            if (j < s) then
                do i = 0, nscl - 1
                    y_jm2(i) = y_jm1(i)
                    y_jm1(i) = y_j(i)
                end do
            end if

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

        do i = 0, nscl - 1
            rkc_step(i) = y_j(i)
        end do

    end function rkc_step

end module ncarles_rkc
