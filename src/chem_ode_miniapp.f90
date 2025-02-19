program chem_ode_miniapp
    !! ADD PROGRAM DOCSTRING(S)
    !!
    use mpi
    use chemistry, only: initialize_chemistry, compute_chemistry
    use integrators, only: initialize_integrator, finalize_integrator, &
                           solve_interval

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
    integer :: nx(3), nx_loc(3)
        !! 3D size of domain
    integer :: nscl, nargs
    integer :: nt, save_unit, nml_unit
    integer :: ix, jy, kz
       !! MPI variables
    integer:: rank, nprocs, ierr, px, py, px_rank, py_rank, comm2d
    

    real, allocatable :: tracers(:, :, :, :), tracers_loc(:, :, :, :), y_0(:), y(:)
        !! 3D reacting scalars state vector and 0D initial condition
    real, allocatable:: args(:, :, :, :), args_loc(:, :, :, :), p_0(:), p(:)
        !! 3D non-reacting scalars vector (e.g., temperature, salinity, etc.)
        !! and it's 0D initial condition

    namelist /params/ integrator, start_time, end_time, save_name, dt_save, &
                      nx, model, temperature, salinity
    namelist /carbonate_ic/ y_0
    namelist /npzd_ic/ y_0
    
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

    ! breaking down the matrices 
    px = floor(sqrt(real(nprocs)))
    do while (mod(nprocs, px) /= 0)
        px = px - 1
    end do
    py = nprocs / px

        ! Create a spatial dimensions
    !CALL MPI_Cart_create(MPI_COMM_WORLD, 2, [px, py], [.TRUE., .TRUE.], .TRUE., comm2d)
    !CALL MPI_Cart_coords(comm2d, rank, 2, [px_rank, py_rank], ierr)

    nx_loc(1) = nx(1) / px
    nx_loc(2) = nx(2) / py
    nx_loc(3) = nx(3)

    ! Configuration and Setup --------------------------------------------------
    ! Read namelists from input file
    open (newunit=nml_unit, file=input_file, status="old")
    read (nml_unit, nml=params)
    rewind(nml_unit)

    ! Initialize chemistry, which associates the `compute_chemistry` pointer
    print *, 'chem model = ', model
    call initialize_chemistry(trim(model), nscl, nargs)

    ! NOTE: allocation of nscl/nargs as intermediate dimension before
    ! z-direction is how NCAR-LES does it currently. This is sure
    ! to be inefficient and should be changed as part of testing.
    ! DON'T FORGET TO CHANGE SAVE_TRACERS AS WELL!
    if (rank == 0) then
        allocate(tracers(nx(1), nx(2), nscl, nx(3)))
        allocate(args(nx(1), nx(2), nargs, nx(3)))
    end if
    
    allocate(tracers_loc(nx_loc(1), nx_loc(2), nscl, nx_loc(3)), y_0(nscl), y(nscl))
    allocate(args_loc(nx_loc(1), nx_loc(2), nargs, nx_loc(3)), p_0(nargs), p(nargs))

    ! Read in the chemical initial condition from the input file
    if (model == 'carbonate') then
        read (nml_unit, nml=carbonate_ic)
        p_0(1) = temperature
        p_0(2) = salinity
    else if (model == 'npzd') then
        read (nml_unit, nml=npzd_ic)
        p_0(1) = temperature
    end if

    ! CHANGE FOR LOOP FOR MPI
    do kz = 1, nx_loc(3)
        do jy = 1, nx_loc(2)
            do ix = 1, nx_loc(1)
                tracers_loc(ix, jy, :, kz) = y_0
                args_loc(ix, jy, :, kz) = p_0
            end do
        end do
    end do

    close (nml_unit)

    ! Initialize the ODE solver, which associates the `solve_interval` pointer
    print *, 'integrator = ', integrator
    call initialize_integrator(integrator, rhs_wrapped, y_0, 1e-8, 1e-6, 1e-10)

    ! Open file for saving tracer history
    open (newunit=save_unit, file=trim(adjustl(save_name)), action="write", status="replace")

    ! Save initial conditions
    ! MAKE SURE THIS TAKES THE WHOLE MATRIX AND AVERAGES IT EVERY TIME STEP IT SAVES
    call save_tracers(time_in_days=.true.)

    ! Time integration loop ----------------------------------------------------
    nt = 0
    do while (time < end_time)
        ! CHANGE FOR LOOP FOR MPI
        do kz = 1, nx_loc(3)
            do jy = 1, nx_loc(2)
                do ix = 1, nx_loc(1)
                    y = tracers_loc(ix, jy, :, kz)
                    p = args_loc(ix, jy, :, kz)
                    call solve_interval(time, time + dt_save, y, p)
                    tracers_loc(ix, jy, :, kz) = y
                end do
            end do
        end do

        call MPI_Gather(tracers_loc, nx_loc(1)*nx_loc(2)*nx_loc(3)*nscl, MPI_REAL, tracers, nx_loc(1)*nx_loc(2)*nx_loc(3)*nscl, MPI_REAL, 0, MPI_COMM_WORLD, ierr)

        if (rank == 0) then
            print *, 'saving output', nt
            call save_tracers(time_in_days=.true.)
        end if

        time = time + dt_save

        nt = nt + 1

    end do

    ! Finalization -------------------------------------------------------------
    close (save_unit)
    call finalize_integrator()
    deallocate (tracers, args, y_0, y, p_0, p, tracers_loc, args_loc)

    call MPI_FINALIZE(ierr)

contains ! ---------------------------------------------------------------------

    subroutine rhs_wrapped(t, y, ydot, p)
        real, intent(in) :: t, y(:)
        real, intent(inout) :: ydot(:)
        real, intent(in), optional :: p(:)
        associate( t => t ); end associate ! suppress unused dummy argument warning
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
        real :: io_tracer(nscl) ! the tracer value at output
        character(len=2) :: str_nscl ! string representation of the integer `nscl`, up to 99
        character(len=:), allocatable :: fmt ! text I/O format string
        integer :: i ! loop index

        io_time = time ! io_time currently in seconds
        if (present(time_in_days) .and. time_in_days) io_time = io_time / SEC_PER_DAY

        do i = 1, nscl
            io_tracer(i) = sum(tracers(1:nx(1), 1:nx(2), i, 1:nx(3))) / (nx(1) * nx(2) * nx(3))
        end do

        write (str_nscl, '(I0)') nscl + 1 ! +1 for time
        fmt = '('//trim(str_nscl)//'ES15.5)'

        write (save_unit, fmt) io_time, io_tracer
    end subroutine save_tracers

end program chem_ode_miniapp
