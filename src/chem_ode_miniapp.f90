program chem_ode_miniapp
    !! ADD PROGRAM DOCSTRING(S)
    !!

    use chemistry, only: carbonate_chem_type
    use rkc_integrator, only: rkc_type

    implicit none ! ------------------------------------------------------------

    character(len=*), parameter :: input_file = "user_inputs.nml"
    character(len=128) :: save_name
        !! output filename
    real :: dt_save = 1e99
        !! data output intervals, [s]
    real :: start_time = 0.0, end_time = 1e-5, time = 0.0
        !! time integration variables, [s]
    real :: temperature = 25.0, salinity = 35.0
        !! temperature [deg C], and salinity [units]
    integer :: nx(3)
        !! 3D size of domain
    integer :: nflat, npts

    integer :: nt, save_unit, nml_unit
    integer :: ix, jy, kz

    type(carbonate_chem_type), target :: chem
    type(rkc_type) :: solver

    real, allocatable, target :: y_3d(:, :, :, :), y_0(:)
        !! 3D reacting scalars state vector and 0D initial condition
    real, allocatable, target :: p_3d(:, :, :, :), p_0(:)
        !! 3D non-reacting scalars vector (e.g., temperature, salinity, etc.)
        !! and it's 0D initial condition
    real, pointer, contiguous :: y_1d(:, :), p_1d(:, :)

    namelist /params/ start_time, end_time, save_name, dt_save, nx, nflat
    namelist /carbonate_ic/ temperature, salinity, y_0

    ! Configuration and Setup --------------------------------------------------
    ! Read namelists from input file
    print *, 'Reading input file'
    open (newunit=nml_unit, file=input_file, status="old")
    read (nml_unit, nml=params)
    rewind (nml_unit)

    allocate (y_3d(nx(1), nx(2), nx(3), chem%nscl), y_0(chem%nscl))
    allocate (p_3d(nx(1), nx(2), nx(3), chem%narg), p_0(chem%narg))

    read (nml_unit, nml=carbonate_ic)
    close (nml_unit)
    p_0(1) = temperature
    p_0(2) = salinity

    print *, 'Filling 3D initial conditions'
    !TODO: Add perturbations to the ICs, like sinusoids or random noise, so that
    !      each spatial point solves a slightly different trajectory in state space
    do kz = 1, nx(3)
        do jy = 1, nx(2)
            do ix = 1, nx(1)
                y_3d(ix, jy, kz, :) = y_0
                p_3d(ix, jy, kz, :) = p_0
            end do
        end do
    end do

    print *, 'Initializing carbonate chemistry object'
    select case(nflat)
    case(1)
        npts = nx(1)
        ! y_1d and p_1d must be contiguous, so in this case, allocate their own
        ! memory instead of using them as pointers
        allocate(y_1d(npts, chem%nscl), p_1d(npts, chem%narg))

    case(2)
        npts = nx(1) * nx(2)
        ! y_1d and p_1d must be contiguous, so in this case, allocate their own
        ! memory instead of using them as pointers
        allocate(y_1d(npts, chem%nscl), p_1d(npts, chem%narg))

    case(3)
        npts = product(nx)
        y_1d(1:npts, 1:chem%nscl) => y_3d ! pointing to whole array is contiguous
        p_1d(1:npts, 1:chem%narg) => p_3d
    end select
    call chem%initialize(npts, y_1d, p_1d)

    print *, 'Initializing RKC integrator object'
    call solver%initialize(chem)

    ! Open file for saving tracer history
    open (newunit=save_unit, file=trim(adjustl(save_name)), action="write", status="replace")

    ! Save initial conditions
    print *, 'Saving initial conditions'
    call save_tracers(time_in_days=.true.)

    ! Time integration loop ----------------------------------------------------
    nt = 1
    do while (time < end_time)

        select case(nflat)
        case(1)
            do kz = 1, nx(3)
                do jy = 1, nx(2)
                    ! solve 3d field one x-vector at a time
                    y_1d(:, :) = y_3d(:, jy, kz, :)
                    p_1d(:, :) = p_3d(:, jy, kz, :)
                    call solver%integrate(time, time + dt_save)
                    y_3d(:, jy, kz, :) = y_1d
                end do
            end do
        case(2)
            do kz = 1, nx(3)
                ! solve 3d field one xy-plane at a time
                ! copy operation requires call to reshape
                y_1d(:, :) = reshape(y_3d(:, :, kz, :), [npts, chem%nscl])
                p_1d(:, :) = reshape(p_3d(:, :, kz, :), [npts, chem%narg])
                call solver%integrate(time, time + dt_save)
                y_3d(:, :, kz, :) = reshape(y_1d, [nx(1), nx(2), chem%nscl])
            end do
        case(3)
            ! directly solve whole 3d field with no external copying!
            call solver%integrate(time, time + dt_save)
        end select

        time = time + dt_save

        print *, 'Saving output', nt
        call save_tracers(time_in_days=.true.)
        nt = nt + 1

    end do

    ! Finalization -------------------------------------------------------------
    print *, 'Finalizing the program.'
    close (save_unit)
    call solver%destroy()
    call chem%destroy()
    deallocate (y_3d, p_3d, y_0, p_0)

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

        write (str_nscl, '(I0)') chem%nscl + 1 ! +1 for time
        fmt = '('//trim(str_nscl)//'ES15.5)'

        write (save_unit, fmt) io_time, y_3d(1, 1, 1, :)
    end subroutine save_tracers

end program chem_ode_miniapp
