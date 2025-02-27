program chem_ode_miniapp
    !! ADD PROGRAM DOCSTRING(S)
    !!
    use mpi
    use iso_fortran_env, only: DP => real64, LI => int64
    use chemistry, only: time_derivative, nscl, nargs
    use miniapp_rkc, only: initialize_rkc, rkc_integrate

    implicit none ! ------------------------------------------------------------
    ! subroutines for gpu

    character(len=*), parameter :: input_file = "./test/user_inputs.nml"
    character(len=128) :: save_name
        !! output filenames
    real :: dt_save = 1e99
        !! data output intervals, [s]
    real :: start_time = 0.0, end_time = 1e-5, time_track = 0.0
        !! time integration variables, [s]
    real :: temperature = 25.0, salinity = 35.0
        !! temperature [deg C], and salinity [units]
    integer :: nx(3), nx_loc(3)
        !! 3D size of domain
    integer :: nt, save_unit, nml_unit, nflat, npts
    integer :: ix, jy, kz, k
    real :: linear_x, linear_y, linear_z, exp_z
       !! MPI variables
    integer :: rank, nprocs, ierr, px, py, px_rank, py_rank, comm2d !, MPI_COMM_WORLD
    integer :: dims(2), coords(2)
    logical :: periods(2)
    real, allocatable :: tracers(:, :, :, :)
        !! 3D reacting scalars state vector and 0D initial condition
    real, allocatable :: args(:, :, :, :)
        !! 3D non-reacting scalars vector (e.g., temperature, salinity, etc.)
        !! and it's 0D initial condition
    real :: y_0(nscl), p_0(nargs), y(nscl), p(nargs)

    integer(LI) :: c0, c1, cr
    real(DP)    :: rate
    character(len=10) :: clock_time

    ! variables for gpu
    !$acc declare create(time_track, dt_save, start_time, end_time, nx, nx_loc, tracers, args, y_0, p_0, y, p)

    namelist /params/ start_time, end_time, save_name, dt_save, &
        nx, nflat, temperature, salinity, y_0
    
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

    nx_loc(1) = nx(1) / px
    nx_loc(2) = nx(2) / py
    nx_loc(3) = nx(3)

    call date_and_time(time=clock_time)
    call system_clock(count_rate=cr)
    rate = real(cr, DP)
    call system_clock(c0)

    write(*, '(A)') '------------------------------------------------'//   &
                    '------------------------------------------------'
    write(*, '(A)') 'MINIAPP started at '//             &
                    clock_time(1:2)//':'//clock_time(3:4)//':'//           &
                    clock_time(5:10)//new_line('a')

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
    allocate (tracers(nscl, nx_loc(1), nx_loc(2), nx_loc(3)))
    allocate (args(nargs, nx_loc(1), nx_loc(2), nx_loc(3)))

    !TODO: Add perturbations to the ICs, like sinusoids or random noise, so that
    !      each spatial point solves a slightly different trajectory in state space
    do kz = 1, nx_loc(3)
        do jy = 1, nx_loc(2)
            do ix = 1, nx_loc(1)
                tracers(:, ix, jy, kz) = y_0
                args(:, ix, jy, kz) = p_0
            end do
        end do
    end do

    close (nml_unit)

    call initialize_rkc(1e-6, 1e-10)

    ! Open file on ROOT for saving tracer history
    if (rank == 0) open (newunit=save_unit, file=trim(adjustl(save_name)), action="write", status="replace")

    ! Compute the averages and save
    call save_tracers(time_in_days=.true.)

!$acc enter data copyin(tracers,args)
    ! Time integration loop ----------------------------------------------------
    nt = 0
    do while (time_track < end_time)

        select case(nflat)
        case(0)
            write(*,*) "Inside Case 0"
!$acc parallel
!$acc loop gang vector collapse(3) private(y,p)
            do kz = 1, nx_loc(3)
                do jy = 1, nx_loc(2)
                    do ix = 1, nx_loc(1)
                      do k=1,nscl
                        y(k) = tracers(k, ix, jy, kz) ! these are not contiguous arrays, must be copied!
                      enddo
                      do k=1,nargs
                        p(k) = args(k, ix, jy, kz) ! these are not contiguous arrays, must be copied!
                      enddo
                        call rkc_integrate(time_track, time_track + dt_save, y, p, npts, nscl, nargs)
                      do k=1,nscl
                        tracers(k, ix, jy, kz) = y(k)
                      enddo
                    end do
                end do
            end do
!$acc end parallel

        case(1)
            write(*,*) "Inside Case 1"
            do kz = 1, nx_loc(3)
                do jy = 1, nx_loc(2)
                    y(:) = reshape(tracers(:, :, jy, kz), [npts*nscl])
                    p(:) = reshape(args(:, :, jy, kz), [npts*nargs])
                    call rkc_integrate(time_track, time_track + dt_save, y, p, npts, nscl, nargs)
                    tracers(:, :, jy, kz) = reshape(y, [nscl, nx_loc(1)])
                end do
            end do

        case(2)
            write(*,*) "Inside Case 2"
            do kz = 1, nx_loc(3)
                y(:) = reshape(tracers(:, :, :, kz), [npts*nscl])
                p(:) = reshape(args(:, :, :, kz), [npts*nargs])
                call rkc_integrate(time_track, time_track + dt_save, y, p, npts, nscl, nargs)
                tracers(:, :, :, kz) = reshape(y, [nscl, nx_loc(1), nx_loc(2)])
            end do

        case(3)
            write(*,*) "Inside Case 3"
            y(:) = reshape(tracers, [npts*nscl])
            p(:) = reshape(args, [npts*nargs])
            call rkc_integrate(time_track, time_track + dt_save, y, p, npts, nscl, nargs)
            tracers(:, :, :, :) = reshape(y, [nscl, nx_loc(1), nx_loc(2), nx_loc(3)])

        end select

        time_track = time_track + dt_save

!$acc update host(tracers)
        print *, 'saving output', nt
        call save_tracers(time_in_days=.true.)

    end do

    ! Finalization -------------------------------------------------------------
    close (save_unit)
!$acc exit data delete(tracers,args)
    deallocate (tracers, args)

    call system_clock(c1)
    call date_and_time(time=clock_time)

    write(*, '(A)') '------------------------------------------------'//   &
                    '------------------------------------------------'
    write(*, '(A)') 'MINIAPP finished at '//             &
                    clock_time(1:2)//':'//clock_time(3:4)//':'//           &
                    clock_time(5:10)//new_line('a')

    WRITE(*,*) "system_clock: ", (c1 - c0) / rate

    call MPI_FINALIZE(ierr)
contains ! ---------------------------------------------------------------------
    subroutine save_tracers(time_in_days)
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

        io_time = time_track ! io_time currently in seconds
        if (present(time_in_days) .and. time_in_days) io_time = io_time / SEC_PER_DAY

        do i = 1, nscl
            io_tracer1(i) = sum(tracers(i, :, :, :))
        end do
        !> This will only change `io_tracer` on rank 0, but printing out to prove it...
        call MPI_Reduce(io_tracer1, io_tracer2, nscl, MPI_REAL8, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
        
        write (str_nscl, '(I0)') nscl + 1 ! +1 for time
        fmt = '(A10,'//trim(str_nscl)//'ES15.5)'

        if (rank == 0) then
            write (save_unit, fmt) 'averages: ', io_time, io_tracer2 / product(nx_loc) ! convert sum to average
            write (save_unit, fmt) 'first pt: ', io_time, tracers(:, 1, 1, 1)
        end if
    end subroutine save_tracers

end program chem_ode_miniapp
