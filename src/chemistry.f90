module chemistry
    implicit none

contains

    pure subroutine compute_chemistry(c_3d, dcdt_3d, p_3d)
        real, contiguous, intent(in) :: c_3d(:, :, :, :), p_3d(:, :, :, :)
        real, contiguous, intent(inout) :: dcdt_3d(:, :, :, :)

        real :: K1s, K2s, Kw, Kb, Rgas, S, T, H_qss !, logT, invT, sqrtS
        real :: a1, a2, a3, a4, a5, a6, a7
        real :: b1, b2, b3, b4, b5, b6, b7
        integer :: ix, jy, kz, nx, ny, nz, nc

        nx = size(c_3d, dim=1)
        ny = size(c_3d, dim=2)
        nz = size(c_3d, dim=3)
        nc = size(c_3d, dim=4)

        do kz = 1, nz
            do jy = 1, ny
                do ix = 1, nx
                    associate(c => c_3d(ix, jy, kz, :), &
                              args => p_3d(ix, jy, kz, :), &
                              dcdt => dcdt_3d(ix, jy, kz, :) &
                              )

                        T = args(1) + 273.15
                        S = args(2)

                        K1s = exp( &
                                (-2307.1266 / T + 2.83655) &
                                - 1.5529413 * log(T) &
                                + (-4.0484 / T - 0.20760841) * (S**0.5) &
                                + 0.08468345 * S &
                                - 0.00654208 * (S**1.5) &
                                + log(1.0 - 0.001005 * S) &
                            ) * (1.0e6)
                        K2s = exp( &
                                (-3351.6106 / T - 9.226508) &
                                - 0.2005743 * log(T) &
                                + (-23.9722 / T - 0.106901773) * (S**0.5) &
                                + 0.1130822 * S &
                                - 0.00846934 * (S**1.5) &
                                + log(1.0 - 0.001005 * S) &
                            ) * (1.0e6)
                        Kw = exp( &
                                (-13847.26 / T + 148.96502) &
                                - 23.65218 * log(T) &
                                + (118.67 / T - 5.977 + 1.0495 * log(T)) * (S**0.5) &
                                - 0.01615 * S &
                            ) * (1.0e6) !(DoE, 1994)
                        Kb = exp( &
                                (-8966.9 - 2890.53 * (S**0.5) &
                                - 77.942 * S &
                                + 1.728 * (S**1.5) &
                                - 0.0996 * (S**2) &
                                ) / T &
                                + (148.0248 + 137.1942 * (S**0.5)) &
                                + 1.62142 * S &
                                - (24.4344 + 25.085 * (S**0.5) + 0.2474 * S) * log(T) &
                                + (0.053105 * (S**0.5) * T) &
                            ) * (1.0e6) !(Dickson, 1990)
                        Rgas = 0.0083143

                        a1 = exp(1246.98 - 6.19 * (10.0**4) / T - 183.0 * log(T))
                        a2 = (4.7e7) * exp(-23.3 / (Rgas * T)) / (1.0e6)
                        a3 = (5.0e10) / (1.0e6)
                        a4 = (6.0e9) / (1.0e6)
                        a5 = (1.4e-3) * (1.0e6)
                        a6 = (4.58e10) * exp(-(20.8 / (Rgas * T))) / (1.0e6)
                        a7 = (3.05e10) * exp(-(20.8 / (Rgas * T))) / (1.0e6)
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

                        dcdt(3) = - a3 * c(3) * H_qss + b3 * c(2) + a4 * c(2) * c(6) - b4 * c(3) &
                                  - a7 * c(3) * c(4) + b7 * c(5) * c(2)

                        dcdt(4) = - a6 * c(4) * c(6) + b6 * c(5) - a7 * c(3) * c(4) + b7 * c(5) * c(2)

                        dcdt(5) = a6 * c(4) * c(6) - b6 * c(5) + a7 * c(3) * c(4) - b7 * c(5) * c(2)

                        dcdt(6) = b2 * c(2) - a2 * c(1) * c(6) - a4 * c(2) * c(6) + b4 * c(3) + a5 &
                                  - b5 * H_qss * c(6) - a6 * c(4) * c(6) + b6 * c(5)

                    end associate
                end do
            end do
        end do

    end subroutine compute_chemistry

end module chemistry
