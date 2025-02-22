program chem_ode_miniapp
    !! ADD PROGRAM DOCSTRING(S)
    !!

    use chemistry, only: initialize_chemistry, compute_chemistry
    ! use integrators, only: initialize_integrator, finalize_integrator, &
    !                        solve_interval
    use pprk4, only: initialize_pprk4, pprk4_integrate, finalize_pprk4

    implicit none ! ------------------------------------------------------------

    character(len=*), parameter :: input_file = "user_inputs.nml"

    character(len=20) :: model, integrator
        !! configuration choices
    character(len=128) :: save_name
        !! output filenames
    real :: dt_save = 1e99
        !! data output intervals, [s]
    real :: start_time = 0.0, end_time = 1e-5, time = 0.0
        !! time integration variables, [s]
    real :: temperature = 25.0, salinity = 35.0
        !! temperature [deg C], and salinity [units]
    integer :: nx(3)
        !! 3D size of domain
    integer :: nscl, nargs
    integer :: nt, save_unit, nml_unit
    integer :: ix, jy, kz

    real, allocatable :: tracers(:, :, :, :), y_0(:)
        !! 3D reacting scalars state vector and 0D initial condition
    real, allocatable :: args(:, :, :, :), p_0(:)
        !! 3D non-reacting scalars vector (e.g., temperature, salinity, etc.)
        !! and it's 0D initial condition

    namelist /params/ integrator, start_time, end_time, save_name, dt_save, &
        nx, model, temperature, salinity
    namelist /carbonate_ic/ y_0
    namelist /npzd_ic/ y_0

    ! Configuration and Setup --------------------------------------------------
    ! Read namelists from input file
    open (newunit=nml_unit, file=input_file, status="old")
    read (nml_unit, nml=params)
    rewind (nml_unit)

    ! Initialize chemistry, which associates the `compute_chemistry` pointer
    print *, 'chem model = ', model
    call initialize_chemistry(trim(model), nscl, nargs)

    ! NOTE: allocation of nscl/nargs as intermediate dimension before
    ! z-direction is how NCAR-LES does it currently. This is sure
    ! to be inefficient and should be changed as part of testing.
    ! DON'T FORGET TO CHANGE SAVE_TRACERS AS WELL!
    allocate (tracers(nx(1), nx(2), nx(3), nscl), y_0(nscl))
    allocate (args(nx(1), nx(2), nx(3), nargs), p_0(nargs))

    ! Read in the chemical initial condition from the input file
    if (model == 'carbonate') then
        read (nml_unit, nml=carbonate_ic)
        p_0(1) = temperature
        p_0(2) = salinity
    else if (model == 'npzd') then
        read (nml_unit, nml=npzd_ic)
        p_0(1) = temperature
    end if

    !TODO: Add perturbations to the ICs, like sinusoids or random noise, so that
    !      each spatial point solves a slightly different trajectory in state space
    do kz = 1, nx(3)
        do jy = 1, nx(2)
            do ix = 1, nx(1)
                tracers(ix, jy, kz, :) = y_0
                args(ix, jy, kz, :) = p_0
            end do
        end do
    end do

    close (nml_unit)

    ! Initialize the ODE solver, which associates the `solve_interval` pointer
    ! print *, 'integrator = ', integrator
    call initialize_pprk4(tracers, rhs_wrapped, init_dt=1e-8)

    ! Open file for saving tracer history
    open (newunit=save_unit, file=trim(adjustl(save_name)), action="write", status="replace")

    ! Save initial conditions
    print *, 'saving initial condition'
    call save_tracers(time_in_days=.true.)

    ! Time integration loop ----------------------------------------------------
    nt = 1
    do while (time < end_time)

        call pprk4_integrate(time, time + dt_save, tracers, args)

        time = time + dt_save

        print *, 'saving output', nt
        call save_tracers(time_in_days=.true.)
        nt = nt + 1

    end do

    ! Finalization -------------------------------------------------------------
    close (save_unit)
    call finalize_integrator()
    deallocate (tracers, args, y_0, p_0)

contains ! ---------------------------------------------------------------------

    subroutine rhs_wrapped(t, y, ydot, p)
        real, intent(in) :: t, y(:, :, :, :)
        real, intent(inout) :: ydot(:, :, :, :)
        real, intent(in), optional :: p(:, :, :, :)
        associate (t => t); end associate ! suppress unused dummy argument warning
        call compute_chemistry(y, ydot, p)
    end subroutine rhs_wrapped

    subroutine save_tracers(time_in_days)
        !! DOCSTRING
        implicit none
        logical, intent(in), optional :: time_in_days
            !! Optional T/F whether to output time in days instead of seconds

        real, parameter :: SEC_PER_DAY = 86400.0

        ! Local variables
        real :: io_time ! the time at output, optionally converted to days
        character(len=2) :: str_nscl ! string representation of the integer `nscl`, up to 99
        character(len=:), allocatable :: fmt ! text I/O format string

        io_time = time ! io_time currently in seconds
        if (present(time_in_days) .and. time_in_days) io_time = io_time / SEC_PER_DAY

        write (str_nscl, '(I0)') nscl + 1 ! +1 for time
        fmt = '('//trim(str_nscl)//'ES15.5)'

        write (save_unit, fmt) io_time, tracers(1, 1, 1, :)
    end subroutine save_tracers

end program chem_ode_miniapp
