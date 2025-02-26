program chem_ode_miniapp
    !! ADD PROGRAM DOCSTRING(S)
    !!

    use chemistry, only: time_derivative, nscl, nargs
    use miniapp_rkc, only: initialize_rkc, rkc_integrate

    implicit none ! ------------------------------------------------------------

    character(len=*), parameter :: input_file = "user_inputs.nml"

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
    integer :: nt, save_unit, nml_unit
    integer :: ix, jy, kz

    real, allocatable :: tracers(:, :, :, :), y_0(nscl), y(nscl)
        !! 3D reacting scalars state vector and 0D initial condition
    real, allocatable :: args(:, :, :, :), p_0(nargs), p(nargs)
        !! 3D non-reacting scalars vector (e.g., temperature, salinity, etc.)
        !! and it's 0D initial condition

    namelist /params/ start_time, end_time, save_name, dt_save, &
        nx, temperature, salinity, y_0

    ! Configuration and Setup --------------------------------------------------
    ! Read namelists from input file
    open (newunit=nml_unit, file=input_file, status="old")
    read (nml_unit, nml=params)
    rewind (nml_unit)

    p_0(1) = temperature
    p_0(2) = salinity

    ! NOTE: allocation of nscl/nargs as intermediate dimension before
    ! z-direction is how NCAR-LES does it currently. This is sure
    ! to be inefficient and should be changed as part of testing.
    ! DON'T FORGET TO CHANGE SAVE_TRACERS AS WELL!
    allocate (tracers(nx(1), nx(2), nscl, nx(3)), y_0(nscl), y(nscl))
    allocate (args(nx(1), nx(2), nargs, nx(3)), p_0(nargs), p(nargs))

    !TODO: Add perturbations to the ICs, like sinusoids or random noise, so that
    !      each spatial point solves a slightly different trajectory in state space
    do kz = 1, nx(3)
        do jy = 1, nx(2)
            do ix = 1, nx(1)
                tracers(ix, jy, :, kz) = y_0
                args(ix, jy, :, kz) = p_0
            end do
        end do
    end do

    close (nml_unit)

    call initialize_rkc(1e-6, 1e-10)

    ! Save initial conditions
    open (newunit=save_unit, file=trim(adjustl(save_name)), action="write", status="replace")
    call save_tracers(time_in_days=.true.)

    ! Time integration loop ----------------------------------------------------
    nt = 0
    do while (time < end_time)

        do kz = 1, nx(3)
            do jy = 1, nx(2)
                do ix = 1, nx(1)
                    y = tracers(ix, jy, :, kz)
                    p = args(ix, jy, :, kz)
                    call rkc_integrate(time_derivative, time, time + dt_save, y, p)
                    tracers(ix, jy, :, kz) = y
                end do
            end do
        end do

        time = time + dt_save

        print *, 'saving output', nt
        call save_tracers(time_in_days=.true.)
        nt = nt + 1

    end do

    ! Finalization -------------------------------------------------------------
    close (save_unit)
    deallocate (tracers, args, y_0, y, p_0, p)

contains ! ---------------------------------------------------------------------

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

        write (save_unit, fmt) io_time, tracers(1, 1, :, 1)
    end subroutine save_tracers

end program chem_ode_miniapp
