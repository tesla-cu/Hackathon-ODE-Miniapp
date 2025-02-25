program chem_ode_miniapp
    !! ADD PROGRAM DOCSTRING(S)
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
    integer :: ixl, ixg, jyl, kzl, ip
       !! MPI variables
    integer :: rank, nprocs, ierr, px, py, px_rank, py_rank, comm2d
    integer :: dims(2), coords(2)
    logical :: periods(2)
    real :: linear_ramp
    real, allocatable :: tracers(:, :, :, :), y_0(:), y(:)
        !! 3D reacting scalars state vector and 0D initial condition
    real, allocatable:: args(:, :, :, :), p_0(:), p(:)
        !! 3D non-reacting scalars vector (e.g., temperature, salinity, etc.)
        !! and it's 0D initial condition

    namelist /params/ integrator, start_time, end_time, save_name, dt_save, &
                      nx, model, temperature, salinity
    namelist /carbonate_ic/ y_0
    namelist /npzd_ic/ y_0

    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)

    dims = [0, 0]
    periods = [.true., .true.]

    call MPI_Dims_create(nprocs, 2, dims, ierr)
    call MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, .true., comm2d, ierr)
    call MPI_Cart_coords(comm2d, rank, 2, coords, ierr)

    px = dims(1)
    py = dims(2)
    px_rank = coords(1)
    py_rank = coords(2)

    ! Configuration and Setup --------------------------------------------------
    ! Read namelists from input file
    open (newunit=nml_unit, file=input_file, status="old")
    read (nml_unit, nml=params)
    rewind(nml_unit)

    ! Initialize chemistry, which associates the `compute_chemistry` pointer
    !print *, 'chem model = ', model
    call initialize_chemistry(trim(model), nscl, nargs)
    
    nx_loc(1) = nx(1) / px
    nx_loc(2) = nx(2) / py
    nx_loc(3) = nx(3)

    allocate(tracers(nx_loc(1), nx_loc(2), nscl, nx_loc(3)), y_0(nscl), y(nscl))
    allocate(args(nx_loc(1), nx_loc(2), nargs, nx_loc(3)), p_0(nargs), p(nargs))

    ! Read in the chemical initial condition from the input file
    if (model == 'carbonate') then
        read (nml_unit, nml=carbonate_ic)
        p_0(1) = temperature
        p_0(2) = salinity
    else if (model == 'npzd') then
        read (nml_unit, nml=npzd_ic)
        p_0(1) = temperature
    end if

    close (nml_unit)
    
    ! Add some pt-to-pt variations, added 'l' to end of indices to be extra clear
    do kzl = 1, nx_loc(3)
        do jyl = 1, nx_loc(2)
            do ixl = 1, nx_loc(1)
                ! ramp all initial conditions from 80% to 120% of nominal value
                ! across the entirety of the x-dimension
                !ixg = nx_loc(1)*px_rank + ixl ! px_rank goes from 0 to px-1
                !linear_ramp = 0.8 + 0.4 * real(ixg-1)/real(nx(1)-1) ! ixg/nx(1) goes from 0.0 to 1.0
                tracers(ixl, jyl, :, kzl) = y_0! * linear_ramp
                args(ixl, jyl, :, kzl) = p_0! * linear_ramp
            end do
        end do
    end do
    ! Initialize the ODE solver, which associates the `solve_interval` pointer
    call initialize_integrator(integrator, rhs_wrapped, y_0, 1e-8, 1e-6, 1e-10)

    ! Open file on ROOT for saving tracer history
    if (rank == 0) open (newunit=save_unit, file=trim(adjustl(save_name)), action="write", status="replace")

    ! Compute the averages and save
    if (rank == 0) print *, 'saving initial conditions'
    call save_tracers(time_in_days=.true., verbose=.true.)

    ! Time integration loop ----------------------------------------------------
    nt = 0
    do while (time < end_time)
        print *, 'in while loop'
        ! CHANGE FOR LOOP FOR MPI
        do kzl = 1, nx_loc(3)
            do jyl = 1, nx_loc(2)
                do ixl = 1, nx_loc(1)
                    y = tracers(ixl, jyl, :, kzl)
                    p = args(ixl, jyl, :, kzl)
                    call solve_interval(time, time + dt_save, y, p)
                    print *, 'finished numerical solver for timestep ', nt
                    tracers(ixl, jyl, :, kzl) = y
                end do
            end do
        end do
        nt = nt + 1
        time =  time + dt_save
        if (rank == 0) print *, 'saving output', nt
        call save_tracers(time_in_days=.true.)

    end do

    ! Finalization -------------------------------------------------------------
    close (save_unit)
    call finalize_integrator()
    deallocate (tracers, args, y_0, y, p_0, p)

    call MPI_FINALIZE(ierr)

contains ! ---------------------------------------------------------------------

    subroutine rhs_wrapped(t, y, ydot, p)
        real, intent(in) :: t, y(:)
        real, intent(inout) :: ydot(:)
        real, intent(in), optional :: p(:)
        associate( t => t ); end associate ! suppress unused dummy argument warning
        call compute_chemistry(y, ydot, p)
    end subroutine rhs_wrapped

    subroutine save_tracers(time_in_days, verbose)
        !! DOCSTRING
        implicit none
        logical, intent(in), optional :: time_in_days, verbose
            !! Optional T/F whether to output time in days instead of seconds

        real, parameter :: SEC_PER_DAY = 86400.0

        ! Local variables
        real :: io_time ! the time at output, optionally converted to days
        real :: io_tracer1(nscl) ! the tracer value at output
        real :: io_tracer2(nscl) ! the tracer value at output
        character(len=2) :: str_nscl ! string representation of the integer `nscl`, up to 99
        character(len=:), allocatable :: fmt ! text I/O format string
        integer :: i ! loop index

        io_time = time ! io_time currently in seconds
        if (present(time_in_days) .and. time_in_days) io_time = io_time / SEC_PER_DAY

        do i = 1, nscl
            io_tracer1(i) = sum(tracers(:, :, i, :))
        end do
        !> This will only change `io_tracer` on rank 0, but printing out to prove it...
        call MPI_Reduce(io_tracer1, io_tracer2, nscl, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        if (present(verbose) .and. verbose) print *, 'Rank ', rank, 'has io_tracer(1) = ', io_tracer2(1)

        write (str_nscl, '(I0)') nscl + 1 ! +1 for time
        fmt = '(A10,'//trim(str_nscl)//'ES15.5)'

        if (rank == 0) then
            write (save_unit, fmt) 'averages: ', io_time, io_tracer2 / product(nx) ! convert sum to average
            write (save_unit, fmt) 'first pt: ', io_time, tracers(1, 1, :, 1)
        end if
    end subroutine save_tracers

end program chem_ode_miniapp
