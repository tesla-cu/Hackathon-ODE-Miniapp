program chem_ode_miniapp
    !! ADD PROGRAM DOCSTRING(S)
    use mpi
    use chemistry, only: time_derivative, nscl, nargs
    use miniapp_rkc, only: initialize_rkc, rkc_integrate

    implicit none ! ------------------------------------------------------------
    ! subroutines for gpu
    !$acc routine (save_tracers) 
    !$acc routine (time_derivative) 
    !$acc routine (initialize_rkc)  
    !$acc routine (rkc_integrate) 

    character(len=*), parameter :: input_file = "user_inputs.nml"

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
    integer :: nt, save_unit, nml_unit
    integer :: ixl, ixg, jyl, kzl, ip, jyg, kzg
    real :: linear_x, linear_y, linear_z, exp_z
       !! MPI variables
    integer :: rank, nprocs, ierr, px, py, px_rank, py_rank, comm2d
    integer :: dims(2), coords(2)
    logical :: periods(2)
    real, allocatable :: tracers(:, :, :, :)
        !! 3D reacting scalars state vector and 0D initial condition
    real, allocatable :: args(:, :, :, :)
        !! 3D non-reacting scalars vector (e.g., temperature, salinity, etc.)
        !! and it's 0D initial condition
    real :: y_0(nscl), p_0(nargs), y(nscl), p(nargs)

    ! variables for gpu
    !$acc declare create(time, dt_save, start_time, end_time, nx, nx_loc, tracers, args, y_0, p_0, y, p)

    namelist /params/ start_time, end_time, save_name, dt_save, &
        nx, temperature, salinity, y_0
    
    call MPI_INIT(ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
   
    ! Assign each MPI process a GPU 
    !$acc set device_num(rank) 

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
    rewind (nml_unit)
    
    nx_loc(1) = nx(1) / px
    nx_loc(2) = nx(2) / py
    nx_loc(3) = nx(3)
    allocate(tracers(nx_loc(1), nx_loc(2), nscl, nx_loc(3)))
    allocate(args(nx_loc(1), nx_loc(2), nargs, nx_loc(3)))

    ! Read in the chemical initial condition from the input file
    read (nml_unit, nml=params)
    p_0(1) = temperature
    p_0(2) = salinity
    close (nml_unit)
    
    ! Add some pt-to-pt variations, added 'l' to end of indices to be extra clear
    !$acc parallel loop collapse(3) copyin(y_0, p_0) copyout(tracers, args)
    do kzl = 1, nx_loc(3)
        !kzg = nx_loc(3) + kzl
        !linear_z = 0.8 + 0.4 * real(kzg-1)/real(nx(3)-1)
        !exp_z = exp(-real(kzg-1)/real(nx(3)-1)) ! -z decay from 0m to -100m 
        do jyl = 1, nx_loc(2)
            !jyg = nx_loc(2)*py_rank + jyl ! px_rank goes from 0 to px-1
            !linear_y = 0.8 + 0.4 * real(jyg-1)/real(nx(2)-1) ! ixg/nx(1) goes from 0.0 to 1.0
            do ixl = 1, nx_loc(1)
                ! ramp all initial conditions from 80% to 120% of nominal value
                ! across the entirety of the x-dimension
                !ixg = nx_loc(1)*px_rank + ixl ! px_rank goes from 0 to px-1
                !linear_x = 0.8 + 0.4 * real(ixg-1)/real(nx(1)-1) ! ixg/nx(1) goes from 0.0 to 1.0
                tracers(ixl, jyl, :, kzl) = y_0(:)
                args(ixl, jyl, :, kzl) = p_0(:)
            end do
        end do
    end do
    !$acc end parallel loop
    !$acc exit data delete(y_0, p_0)
    ! Initialize the ODE solver, which associates the `solve_interval` pointer
    call initialize_rkc(1e-6, 1e-10)

    ! Open file on ROOT for saving tracer history
    if (rank == 0) open (newunit=save_unit, file=trim(adjustl(save_name)), action="write", status="replace")

    ! Compute the averages and save
    call save_tracers(time_in_days=.true.)

    ! Time integration loop ----------------------------------------------------
    nt = 0
    do while (time < end_time)
        ! CHANGE FOR LOOP FOR MPI
        !$acc parallel loop collapse(3) copyin(tracers, args, p, y)
        do kzl = 1, nx_loc(3)
            do jyl = 1, nx_loc(2)
                do ixl = 1, nx_loc(1)
                    p = args(ixl, jyl, :, kzl)
                    y = tracers(ixl, jyl, :, kzl)
                    call rkc_integrate(time_derivative, time, time + dt_save, y, p)
                    tracers(ixl, jyl, :, kzl) = y
                end do
            end do
        end do
        !$acc end parallel loop copyout(tracers)
        nt = nt + 1
        time =  time + dt_save
        !$acc kernels
        if (rank == 0) print *, 'saving output', nt
        !$acc end kernels
        call save_tracers(time_in_days=.true.)

    end do

    ! Finalization -------------------------------------------------------------
    close (save_unit)
    deallocate (tracers, args)

    call MPI_FINALIZE(ierr)
contains ! ---------------------------------------------------------------------
    subroutine save_tracers(time_in_days)
        !$acc routine seq 
        !! DOCSTRING
        implicit none
        logical, intent(in), optional :: time_in_days
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
        
        write (str_nscl, '(I0)') nscl + 1 ! +1 for time
        fmt = '(A10,'//trim(str_nscl)//'ES15.5)'

        if (rank == 0) then
            write (save_unit, fmt) 'averages: ', io_time, io_tracer2 / product(nx) ! convert sum to average
            write (save_unit, fmt) 'first pt: ', io_time, tracers(1, 1, :, 1)
        end if
    end subroutine save_tracers

end program chem_ode_miniapp
