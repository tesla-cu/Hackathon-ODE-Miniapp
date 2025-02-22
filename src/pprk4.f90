module pprk4
    !! A positivity-preserving RK4 integrator for systems of ODEs with
    !! strictly-positive solution spaces, such as chemical reactions.

    implicit none
    private
    public :: initialize_pprk4, pprk4_integrate, pprk4_step, rk4_step, finalize_pprk4

    real, parameter :: a(4) = [1./6., 1./3., 1./3., 1./6.]
    real, parameter :: b(4) = [0.5, 0.5, 1.0, 0.0]
    real, parameter :: rtol = 10.0 * epsilon(1.0)
    real, allocatable :: temp(:, :, :, :), y_new(:, :, :, :), dydt(:, :, :, :)
    real :: dt_old

    procedure(time_derivative), pointer :: rhs
    abstract interface
        subroutine time_derivative(t, y, ydot, p)
            real, allocatable, intent(in) :: t, y(:, :, :, :)
            real, allocatable, intent(inout) :: ydot(:, :, :, :)
            real, allocatable, intent(in), optional :: p(:, :, :, :)
        end subroutine time_derivative
    end interface

contains

    subroutine initialize_pprk4(user_y, user_rhs, init_dt)
        !! Initialize the RK4 integrator's working memory and minimum timestep
        !! size (optional).

        real, intent(in) :: user_y(:, :, :, :)
            !! User-supplied solution vector
        procedure(time_derivative) :: user_rhs
            !! User-supplied RHS term of ODE, dy/dt = RHS
        real, intent(in), optional :: init_dt
            !! Optional initial timestep size

        allocate (temp, y_new, dydt, mold=user_y)
        rhs => user_rhs
        dt_old = rtol; if (present(init_dt)) dt_old = init_dt

    end subroutine initialize_pprk4

    subroutine pprk4_integrate(ti, tf, y, p)
        !! Integrate from ti to tf, using adaptive RK4.

        real, intent(in) :: ti, tf
            !! the initial and final times over which to evolve the solution.
        real, intent(inout) :: y(:, :, :, :)
            !! Solution vector at initial time `ti` on input, and at final time
            !! `tf` on output.
        real, intent(in), optional :: p(:, :, :, :)
            !! user-supplied extra arguments (aka [p]arameters) to the RHS function

        real :: t, dt

        t = ti
        dt = dt_old
        do
            dt = max(tf * rtol, min(dt, tf - t)) ! same as: `if (tf - t < dt) dt = tf - t`
            call pprk4_step(t, dt, y, p) ! t, dt, and y all updated here
            if (tf - t <= tf * rtol) exit
        end do
        dt_old = dt

    end subroutine pprk4_integrate

    subroutine pprk4_step(t, dt, y, p)
        !! Integrate forward one RK4 step, using a step-size of `dt` or smaller,
        !! as needed, to preserve the positivity of the solution.

        real, intent(inout) :: t, dt
            !! Current time of the solution on input, final time of RK4 step on output
        real, intent(inout) :: y(:, :, :, :)
            !! Solution vector at current time on input, at end of step on output
        real, intent(in), optional :: p(:, :, :, :)
            !! user-supplied extra arguments (aka [p]arameters) to the RHS function

        real :: dt_new
        integer :: irk

        ! take a step at largest-possible dt, trying multiple times if necessary
        do
            ! Integrate one full RK4 step
            temp(...) = y ! stores intermediate stage states
            y_new(...) = y ! accumulates updated state
            do irk = 1, 4
                call rhs(t, temp, dydt, p)
                temp(...) = y + b(irk) * dt * dydt
                y_new(...) = y_new + a(irk) * dt * dydt
            end do

            ! check positivity of solution
            if (minval(y_new) > 0.0) then
                ! end the loop
                y(...) = y_new
                exit
            else
                ! redo the step with dt = 0.5 * dt
                if (dt <= t * rtol) stop "RUNTIME ERROR: RK4 integration unstable at minimum dt value!"
                dt = max(0.5 * dt, t * rtol)
            end if
        end do

        ! update time
        t = t + dt
        ! compute new dt using a rate-based method whereby no member of y will
        ! grow or shrink more than 10% in a single timestep
        dt_new = max(0.1 * minval(abs(y / dydt)), t * rtol)
        dt = min(dt_new, 2.0 * dt) ! don't let dt grow too fast

    end subroutine pprk4_step

    subroutine rk4_step(t, dt, y, p)
        !! Take a single plain RK4 integration step, of prescribed size dt.
        real, intent(in) :: t, dt
            !! current time and step size
        real, intent(inout) :: y(:)
            !! Solution vector at current time on input, at time+dt on output
        real, intent(in), optional :: p(:)
            !! user-supplied extra arguments (aka [p]arameters) to the RHS function

        integer :: irk

        ! copy current state into RK4 working memory
        temp = y ! stores intermediate stage states
        y_new = y ! accumulates stages into updated state

        ! Integrate one full RK4 step
        do irk = 1, 4
            call rhs(t, temp, dydt, p)
            temp = y + b(irk) * dt * dydt
            y_new = y_new + a(irk) * dt * dydt
        end do
        y = y_new

    end subroutine rk4_step

    subroutine finalize_pprk4()
        !! Deallocate RK4 working memory
        deallocate (temp, y_new, dydt)
    end subroutine finalize_pprk4

end module pprk4
