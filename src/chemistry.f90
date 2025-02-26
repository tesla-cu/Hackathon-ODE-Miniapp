module chemistry
    implicit none
    public

    integer, parameter :: nscl = 6, nargs = 2
    real, parameter :: SEC_PER_DAY = 86400.0

contains
    subroutine time_derivative(t, y, ydot, p)
        !$acc rountine seq
        real, intent(in) :: t
        real, intent(in) :: y(nscl), p(nargs)
        real, intent(inout) :: ydot(nscl)

        real :: K1s, K2s, Kw, Kb, Rgas, salt, temp, H_qss
        real :: a1, a2, a3, a4, a5, a6, a7
        real :: b1, b2, b3, b4, b5, b6, b7

        associate (t => t); end associate ! suppress unused dummy argument warning

        temp = p(1) + 273.15
        salt = p(2)

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

        H_qss = (a1 * y(1) + b3 * y(2) + a5) / (b1 * y(2) + a3 * y(3) + b5 * y(6))

        ydot(1) = b1 * y(2) * H_qss + b2 * y(2) - a1 * y(1) - a2 * y(1) * y(6)

        ydot(2) = a1 * y(1) + a2 * y(1) * y(6) - b1 * y(2) * H_qss - b2 * y(2) &
                  + a3 * y(3) * H_qss - b3 * y(2) - a4 * y(2) * y(6) + b4 * y(3) &
                  + a7 * y(3) * y(4) - b7 * y(5) * y(2)

        ydot(3) = -a3 * y(3) * H_qss + b3 * y(2) + a4 * y(2) * y(6) - b4 * y(3) &
                  - a7 * y(3) * y(4) + b7 * y(5) * y(2)

        ydot(4) = -a6 * y(4) * y(6) + b6 * y(5) - a7 * y(3) * y(4) + b7 * y(5) * y(2)

        ydot(5) = a6 * y(4) * y(6) - b6 * y(5) + a7 * y(3) * y(4) - b7 * y(5) * y(2)

        ydot(6) = b2 * y(2) - a2 * y(1) * y(6) - a4 * y(2) * y(6) + b4 * y(3) + a5 &
                  - b5 * H_qss * y(6) - a6 * y(4) * y(6) + b6 * y(5)

    end subroutine time_derivative

end module chemistry
