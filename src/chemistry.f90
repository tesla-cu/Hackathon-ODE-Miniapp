module chemistry
    implicit none
    public

    integer, parameter :: nscl = 6, nargs = 2
    real, parameter :: SEC_PER_DAY = 86400.0

contains

    pure subroutine time_derivative(t, y, ydot, p)
        real, intent(in) :: t
        real, intent(in) :: y(0:) ! explicitly setting lbound inside this routine
        real, intent(out) :: ydot(0:) ! same
        real, intent(in), optional :: p(0:) ! same

        real :: c(nscl), dcdt(nscl)
        integer :: npts, ipt, ic, iarg
        real :: K1s, K2s, Kw, Kb, Rgas, salt, temp, H_qss
        real :: a1, a2, a3, a4, a5, a6, a7
        real :: b1, b2, b3, b4, b5, b6, b7

        associate (t => t); end associate
        npts = size(y) / nscl

        !$acc parallel loop private(dcdt, c)
        do ipt = 0, npts-1
            ic = ipt*nscl
            iarg = ipt*nargs
            c(:) = y(ic:ic+nscl-1)
            dcdt(:) = ydot(ic:ic+nscl-1)
            temp = p(iarg) + 273.15
            salt = p(iarg+1)

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

            ydot(ic:ic+nscl-1) = dcdt
        end do
    end subroutine time_derivative

end module chemistry
