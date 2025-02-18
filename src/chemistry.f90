module chemistry
    implicit none
    private
    public :: initialize_chemistry, compute_chemistry

    real, parameter :: SEC_PER_DAY = 86400.0

    abstract interface
        subroutine time_derivative(tracers, dcdt, args)
            real, intent(in) :: tracers(:), args(:)
            real, intent(inout) :: dcdt(:)
        end subroutine time_derivative
    end interface
    procedure(time_derivative), pointer, protected :: compute_chemistry => null()

contains

    subroutine initialize_chemistry(model, nscl, nargs)
        implicit none

        character(len=*), intent(in) :: model
        integer, intent(out) :: nscl, nargs

        if (model == 'carbonate') then
            nscl = 6
            nargs = 2

            compute_chemistry => dcdt_carbonate

        else if (model == 'NPZD') then
            nscl = 4
            nargs = 1

            compute_chemistry => dcdt_npzd

        else
            stop 'ERROR: chemistry model not recognized'
        end if

    end subroutine initialize_chemistry

    subroutine dcdt_carbonate(c, dcdt, args)
        real, intent(in) :: c(:), args(:)
        real, intent(inout) :: dcdt(:)

        real :: K1s, K2s, Kw, Kb, Rgas, salt, temp, H_qss
        real :: a1, a2, a3, a4, a5, a6, a7
        real :: b1, b2, b3, b4, b5, b6, b7

        temp = args(1) + 273.15
        salt = args(2)

        K1s = exp(-2307.1266 / temp + 2.83655 - 1.5529413 * log(temp) + &
                  (-4.0484 / temp - 0.20760841) * (salt**0.5) + 0.08468345 * salt - &
                  0.00654208 * (salt**1.5) + log(1.0 - 0.001005 * salt)) * (1.0e6)
        K2s = exp(-3351.6106 / temp - 9.226508 - 0.2005743 * log(temp) + &
                  (-23.9722 / temp - 0.106901773) * (salt**0.5) + 0.1130822 * salt - &
                  0.00846934 * (salt**1.5) + log(1.0 - 0.001005 * salt)) * (1.0e6)
        Kw = exp(148.96502 - 13847.26 / temp - 23.65218 * log(temp) + &
                 (118.67 / temp - 5.977 + 1.0495 * log(temp)) * (salt**0.5) - &
                 0.01615 * salt) * (1.0e6) !(DoE, 1994)
        Kb = exp((-8966.9 - 2890.53 * (salt**0.5) - 77.942 * salt + &
                  1.728 * (salt**1.5) - 0.0996 * (salt**2)) / temp &
                 + 148.0248 + 137.1942 * (salt**0.5) + 1.62142 * salt - &
                 (24.4344 + 25.085 * (salt**0.5) + 0.2474 * salt) * log(temp) + &
                 0.053105 * (salt**0.5) * temp) * (1.0e6) !(Dickson, 1990)
        Rgas = 0.0083143

        a1 = exp(1246.98 - 6.19 * (10.0**4) / temp - 183.0 * log(temp))
        a2 = (4.7e7) * exp(-23.3 / (Rgas * temp)) / (1.0e6)
        a3 = (5.0e10) / (1.0e6)
        a4 = (6.0e9) / (1.0e6)
        a5 = (1.4e-3) * (1.0e6)
        a6 = (4.58e10) * exp(-(20.8 / (Rgas * temp))) / (1.0e6)
        a7 = (3.05e10) * exp(-(20.8 / (Rgas * temp))) / (1.0e6)
        b1 = a1 / K1s
        b2 = (Kw * a2 / K1s) * (1.0e6)
        b3 = a3 * K2s
        b4 = (a4 * Kw / K2s) * (1.0e6)
        b5 = (a5 / Kw) / (1.0e6)
        b6 = (a6 * Kw / Kb) * (1.0e6)
        b7 = a7 * K2s / Kb

        H_qss = (a1 * c(1) + b3 * c(2) + a5) / (b1 * c(2) + a3 * c(3) + b5 * c(6))

        dcdt(1) = b1 * c(2) * H_qss + b2 * c(2) - a1 * c(1) - a2 * c(1) * c(6)

        dcdt(2) = a1 * c(1) + a2 * c(1) * c(6) - b1 * c(2) * H_qss - b2 * c(2) &
                  + a3 * c(3) * H_qss - b3 * c(2) - a4 * c(2) * c(6) + b4 * c(3) &
                  + a7 * c(3) * c(4) - b7 * c(5) * c(2)

        dcdt(3) = -a3 * c(3) * H_qss + b3 * c(2) + a4 * c(2) * c(6) - b4 * c(3) &
                  - a7 * c(3) * c(4) + b7 * c(5) * c(2)

        dcdt(4) = -a6 * c(4) * c(6) + b6 * c(5) - a7 * c(3) * c(4) + b7 * c(5) * c(2)

        dcdt(5) = a6 * c(4) * c(6) - b6 * c(5) + a7 * c(3) * c(4) - b7 * c(5) * c(2)

        dcdt(6) = b2 * c(2) - a2 * c(1) * c(6) - a4 * c(2) * c(6) + b4 * c(3) + a5 &
                  - b5 * H_qss * c(6) - a6 * c(4) * c(6) + b6 * c(5)

    end subroutine dcdt_carbonate

    subroutine dcdt_npzd(tracers, dcdt, args)
        ! NPZ from P Franks 1986 recommended by Nikki Lovenduski
        ! Parameters
        real, intent(in) :: tracers(:), args(:)
        real, intent(inout) :: dcdt(:)
        real :: vp, intensity, temp, light_intensity
        real :: kn = 1.0        ! umolN/l
        real :: rm = 1.0        ! 1/d
        real :: death_rate_zoo = 0.2    ! 1/d
        real :: lambda = 0.2    ! umolN/l
        real :: death_rate_phyto = 0.1  ! 1/d
        real :: alpha = 0.3
        real :: beta = 0.6
        real :: phi = 0.4 ! 1/d
        real :: r_npzd = 0.15  ! 1/d
        real :: a_npz = 0.6 ! 1/d
        real :: b_npz = 1.066
        real :: c_npz = 1.0
        real :: P, Z, N, D

        P = tracers(1)
        Z = tracers(2)
        N = tracers(3)
        D = tracers(4)
        temp = args(1)

        light_intensity = 1.0
        !intensity = rm * P * lambda !Mayzaud-Poulet (1/d)
        intensity = rm !Ivlev (1/d)
        vp = (a_npz * b_npz**(c_npz * temp)) !from Eppley 1972 (1/d)

        dcdt(1) = vp * (N / (kn + N)) * light_intensity * P &
                  - intensity * (1.0 - exp(-lambda * P)) * Z &
                  - death_rate_phyto * P - r_npzd * P

        dcdt(2) = beta * intensity * (1.0 - exp(-lambda * P)) * Z - death_rate_zoo * Z

        dcdt(3) = -vp * (N / (kn + N)) * light_intensity * P &
                  + alpha * intensity * (1.0 - exp(-lambda * P)) * Z + death_rate_phyto * P &
                  + death_rate_zoo * Z + phi * D

        dcdt(4) = r_npzd * P + (1 - alpha - beta) * intensity * (1.0 - exp(-lambda * P)) * Z &
                  - phi * D

        dcdt = dcdt / SEC_PER_DAY

    end subroutine dcdt_npzd

end module chemistry
