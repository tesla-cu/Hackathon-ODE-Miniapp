module integrators
    !! Interface connecting Miniapp to ODE solver packages ("integrators")
    use miniapp_rkc, only: initialize_rkc, rkc_integrate, &
                           rkc_step => rkc_inplace_step
    use ncarles_rkc, only: initialize_rkc77 => initialize_rkc, &
                           rkc77_integrate => rkc_integrate, &
                           rkc77_step => rkc_inplace_step
    use pprk4, only: initialize_pprk4, pprk4_integrate, rk4_step, finalize_pprk4
    implicit none
    private
    public :: initialize_integrator, finalize_integrator, solve_interval, solve_step

    character(len=:), allocatable :: solver

    abstract interface
        subroutine time_derivative(t, y, ydot, p)
            real, intent(in) :: t, y(:)
            real, intent(inout) :: ydot(:)
            real, intent(in), optional :: p(:)
        end subroutine time_derivative

        subroutine interval(ti, tf, y, p)
            real, intent(in) :: ti, tf
            real, intent(inout) :: y(:)
            real, intent(in), optional :: p(:)
        end subroutine interval

        subroutine step(t, dt, y, p)
            real, intent(in) :: t, dt
            real, intent(inout) :: y(:)
            real, intent(in), optional :: p(:)
        end subroutine step
    end interface
    procedure(interval), pointer, protected :: solve_interval => null()
    procedure(step), pointer, protected :: solve_step => null()

contains ! -----------------------------------------------------------------

    subroutine initialize_integrator(choice, dydt, y, init_dt, rtol, atol)
        !! Initialize the time integration

        implicit none
        character(len=*), intent(in) :: choice
            !! user-selected solver type
        procedure(time_derivative) :: dydt
            !! right-hand side of the program-defined system of differential equation
        real, intent(in) :: y(:)
            !! pre-allocated solution array
        real, intent(in), optional :: init_dt, rtol, atol
            !! initial timestep size, relative error tolerance, absolute error tolerance

        solver = trim(adjustl(choice))

        if (solver == 'pprk4') then
            call initialize_pprk4(y, dydt, init_dt)
            solve_interval => pprk4_integrate
            solve_step => rk4_step

        else if (solver == 'rkc') then
            call initialize_rkc(dydt, rtol, atol)
            solve_interval => rkc_integrate
            solve_step => rkc_step

        else if (solver == 'rkc77') then
            call initialize_rkc77(dydt, rtol, atol)
            solve_interval => rkc77_integrate
            solve_step => rkc77_step

        else
            stop 'ERROR: Integrator choice not recognized!'

        end if

    end subroutine initialize_integrator

    subroutine finalize_integrator()
        !! Finalize the time integration (only necessary for rk4)
        implicit none

        if (solver == 'pprk4') then
            call finalize_pprk4()

        else if (solver(1:3) /= 'rkc') then
            stop 'ERROR: Integrator choice not recognized!'

        end if

    end subroutine finalize_integrator

end module integrators
