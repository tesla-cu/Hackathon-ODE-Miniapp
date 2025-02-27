program chem_ode_miniapp
    !! ADD PROGRAM DOCSTRING(S)
    !!
    use iso_fortran_env, only: DP => real64, LI => int64
    use chemistry, only: nscl, nargs
    use miniapp_rkc, only: initialize_rkc, rkc_integrate

    implicit none ! ------------------------------------------------------------

    character(len=*), parameter :: input_file = "./test/user_inputs.nml"

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
    integer :: nt, save_unit, nml_unit, nflat, npts
    integer :: ix, jy, kz, k

    real, allocatable :: tracers(:, :, :, :)
        !! 3D reacting scalars state vector
    real, allocatable :: args(:, :, :, :)
        !! 3D non-reacting scalars vector (e.g., temperature, salinity, etc.)
    real :: y_0(nscl), y(nscl)
    real :: p_0(nargs), p(nargs)

    integer(LI) :: c0, c1, cr
    real(DP)    :: rate
    character(len=10) :: clock_time

    namelist /params/ start_time, end_time, save_name, dt_save, &
        nx, nflat, temperature, salinity, y_0

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
    close (nml_unit)

    p_0(1) = temperature
    p_0(2) = salinity

    ! NOTE: allocation of nscl/nargs as intermediate dimension before
    ! z-direction is how NCAR-LES does it currently. This is sure
    ! to be inefficient and should be changed as part of testing.
    ! DON'T FORGET TO CHANGE SAVE_TRACERS AS WELL!
    allocate (tracers(nscl, nx(1), nx(2), nx(3)))
    allocate (args(nargs, nx(1), nx(2), nx(3)))

    !TODO: Add perturbations to the ICs, like sinusoids or random noise, so that
    !      each spatial point solves a slightly different trajectory in state space
    do kz = 1, nx(3)
        do jy = 1, nx(2)
            do ix = 1, nx(1)
                tracers(:, ix, jy, kz) = y_0
                args(:, ix, jy, kz) = p_0
            end do
        end do
    end do

    select case(nflat)
    case(0)
        npts = 1

    case(1)
        npts = nx(1)

    case(2)
        npts = nx(1) * nx(2)

    case(3)
        npts = product(nx)
    end select
!!!    allocate (y(npts*nscl), p(npts*nargs))

    call initialize_rkc(1e-6, 1e-10)

    ! Save initial conditions
    open (newunit=save_unit, file=trim(adjustl(save_name)), action="write", status="replace")
    call save_tracers(time_in_days=.true.)

!$acc enter data copyin(tracers,args)
    ! Time integration loop ----------------------------------------------------
    nt = 0
    do while (time < end_time)

        select case(nflat)
        case(0)
            write(*,*) "Inside Case 0"
!$acc parallel
!$acc loop gang vector collapse(3) private(y,p)
            do kz = 1, nx(3)
                do jy = 1, nx(2)
                    do ix = 1, nx(1)
                      do k=1,nscl
                        y(k) = tracers(k, ix, jy, kz) ! these are not contiguous arrays, must be copied!
                      enddo
                      do k=1,nargs
                        p(k) = args(k, ix, jy, kz) ! these are not contiguous arrays, must be copied!
                      enddo
                        call rkc_integrate(time, time + dt_save, y, p, npts, nscl, nargs)
                      do k=1,nscl
                        tracers(k, ix, jy, kz) = y(k)
                      enddo
                    end do
                end do
            end do
!$acc end parallel

        case(1)
            write(*,*) "Inside Case 1"
            do kz = 1, nx(3)
                do jy = 1, nx(2)
                    y(:) = reshape(tracers(:, :, jy, kz), [npts*nscl])
                    p(:) = reshape(args(:, :, jy, kz), [npts*nargs])
                    call rkc_integrate(time, time + dt_save, y, p, npts, nscl, nargs)
                    tracers(:, :, jy, kz) = reshape(y, [nscl, nx(1)])
                end do
            end do

        case(2)
            write(*,*) "Inside Case 2"
            do kz = 1, nx(3)
                y(:) = reshape(tracers(:, :, :, kz), [npts*nscl])
                p(:) = reshape(args(:, :, :, kz), [npts*nargs])
                call rkc_integrate(time, time + dt_save, y, p, npts, nscl, nargs)
                tracers(:, :, :, kz) = reshape(y, [nscl, nx(1), nx(2)])
            end do

        case(3)
            write(*,*) "Inside Case 3"
            y(:) = reshape(tracers, [npts*nscl])
            p(:) = reshape(args, [npts*nargs])
            call rkc_integrate(time, time + dt_save, y, p, npts, nscl, nargs)
            tracers(:, :, :, :) = reshape(y, [nscl, nx(1), nx(2), nx(3)])

        end select

        time = time + dt_save

!$acc update host(tracers)
        print *, 'saving output', nt
        call save_tracers(time_in_days=.true.)
        nt = nt + 1

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

        write (save_unit, fmt) io_time, tracers(:, 1, 1, 1)
    end subroutine save_tracers

end program chem_ode_miniapp
