module chemistry
    use integrand, only: integrand_type, IK, RK
    implicit none
    private
    public :: carbonate_chem_type

    type, extends(integrand_type) :: carbonate_chem_type
        private
        integer(IK), public :: nscl = 6, narg = 2
        integer(IK), public :: npts
        real(RK), contiguous, pointer, public :: species(:, :) => null()
            !! objects may or may not "own" this data, but symmetric assignment forces a data copy
        real(RK), contiguous, pointer, public :: args(:, :) => null()
            !! objects may or may not "own" this data, and symmetric assignment associates pointers

    contains
        procedure, pass(self), public :: initialize
        procedure, pass(self), public :: destroy
        procedure, pass(self), public :: size
        procedure, pass(self), public :: state
        procedure, pass(self), public :: d_dt => compute_chemistry !< Time derivative
        ! procedure, pass(lhs) :: local_error
        ! +
        procedure, pass(lhs) :: integrand_add_integrand !< `+` operator.
        procedure, pass(lhs) :: integrand_add_real      !< `+ real` operator.
        procedure, pass(rhs) :: real_add_integrand      !< `real +` operator.
        ! *
        procedure, pass(lhs) :: integrand_multiply_integrand   !< `*` operator.
        procedure, pass(lhs) :: integrand_multiply_real        !< `* real` operator.
        procedure, pass(rhs) :: real_multiply_integrand        !< `real *` operator.
        procedure, pass(lhs) :: integrand_multiply_real_scalar !< `* real_scalar` operator.
        procedure, pass(rhs) :: real_scalar_multiply_integrand !< `real_scalar *` operator.
        ! -
        procedure, pass(lhs) :: integrand_sub_integrand !< `-` operator.
        procedure, pass(lhs) :: integrand_sub_real      !< `- real` operator.
        procedure, pass(rhs) :: real_sub_integrand      !< `real -` operator.
        ! =
        procedure, pass(lhs) :: integrand_eq_integrand !< `=` operator.
        procedure, pass(lhs) :: integrand_eq_real      !< `= real` operator.
    end type carbonate_chem_type

contains

    subroutine initialize(self, npts, species, args)
        !< Initialize integrand.
        class(carbonate_chem_type), intent(inout) :: self
        integer(IK), intent(in) :: npts
        real(RK), contiguous, target, intent(in), optional :: species(:, :)
        real(RK), contiguous, target, intent(in), optional :: args(:, :)

        call self%destroy()
        self%npts = npts

        ! first associate pointers with optional inputs
        if (present(species)) self%species => species
        if (present(args)) self%args => args

        ! then, if not associated, allocate memory for the pointers
        if (.not. associated(self%species)) allocate(self%species(npts, self%nscl))
        if (.not. associated(self%args)) allocate(self%args(npts, self%narg))
    end subroutine initialize

    subroutine destroy(self)
        !< Destroy integrand.
        !< WARNING: DOUBLE CHECK THIS FOR MEMORY LEAKING, CONFUSED ABOUT STANDARD!
        class(carbonate_chem_type), intent(inout) :: self

        self%npts = -1
        if (associated(self%species)) nullify(self%species)
        if (associated(self%args)) nullify(self%args)
    end subroutine destroy

    pure function size(self)
        class(carbonate_chem_type), intent(in) :: self
        integer(IK) :: size
        size = self%npts * self%nscl
    end function

    function state(self)
        !! Reference a 1D pointer to the 2D species array
        class(carbonate_chem_type), intent(in) :: self
        real(RK), contiguous, pointer :: state(:)
        state(1:self%size()) => self%species
    end function

    pure function compute_chemistry(self, t) result(dState_dt)
        class(carbonate_chem_type), intent(in) :: self !< Integrand object.
        real(RK), intent(in), optional :: t    !< time of integration (unused)
        real(RK), allocatable, target :: dState_dt(:)   !< Result as 1D array

        real(RK), contiguous, pointer :: dSpecies_dt(:, :) !< Result reshaped as 2D array
        real :: K1s, K2s, Kw, Kb, Rgas, Salinity, Temperature, H_qss
        real :: a1, a2, a3, a4, a5, a6, a7
        real :: b1, b2, b3, b4, b5, b6, b7
        integer :: ipt
        logical :: time_check

        time_check = present(t) ! suppresses unused argument warning
        allocate(dState_dt(self%size()))
        dSpecies_dt(1:self%npts, 1:self%nscl) => dState_dt ! reshape dState_dt via pointer reference

        do ipt = 1, self%npts
            associate(c => self%species(ipt, :), &
                      args => self%args(ipt, :), &
                      dcdt => dSpecies_dt(ipt, :) & ! does the correct indexing into opr for us
                    )

                Temperature = args(1) + 273.15
                Salinity = args(2)

                K1s = exp( &
                        (-2307.1266 / Temperature + 2.83655) &
                        - 1.5529413 * log(Temperature) &
                        + (-4.0484 / Temperature - 0.20760841) * (Salinity**0.5) &
                        + 0.08468345 * Salinity &
                        - 0.00654208 * (Salinity**1.5) &
                        + log(1.0 - 0.001005 * Salinity) &
                        ) * (1.0e6)
                K2s = exp( &
                        (-3351.6106 / Temperature - 9.226508) &
                        - 0.2005743 * log(Temperature) &
                        + (-23.9722 / Temperature - 0.106901773) * (Salinity**0.5) &
                        + 0.1130822 * Salinity &
                        - 0.00846934 * (Salinity**1.5) &
                        + log(1.0 - 0.001005 * Salinity) &
                        ) * (1.0e6)
                Kw = exp( &
                        (-13847.26 / Temperature + 148.96502) &
                        - 23.65218 * log(Temperature) &
                        + (118.67 / Temperature - 5.977 + 1.0495 * log(Temperature)) * (Salinity**0.5) &
                        - 0.01615 * Salinity &
                        ) * (1.0e6) !(DoE, 1994)
                Kb = exp( &
                        (-8966.9 - 2890.53 * (Salinity**0.5) &
                        - 77.942 * Salinity &
                        + 1.728 * (Salinity**1.5) &
                        - 0.0996 * (Salinity**2) &
                        ) / Temperature &
                        + (148.0248 + 137.1942 * (Salinity**0.5)) &
                        + 1.62142 * Salinity &
                        - (24.4344 + 25.085 * (Salinity**0.5) + 0.2474 * Salinity) * log(Temperature) &
                        + (0.053105 * (Salinity**0.5) * Temperature) &
                        ) * (1.0e6) !(Dickson, 1990)
                Rgas = 0.0083143

                a1 = exp(1246.98 - 6.19 * (10.0**4) / Temperature - 183.0 * log(Temperature))
                a2 = (4.7e7) * exp(-23.3 / (Rgas * Temperature)) / (1.0e6)
                a3 = (5.0e10) / (1.0e6)
                a4 = (6.0e9) / (1.0e6)
                a5 = (1.4e-3) * (1.0e6)
                a6 = (4.58e10) * exp(-(20.8 / (Rgas * Temperature))) / (1.0e6)
                a7 = (3.05e10) * exp(-(20.8 / (Rgas * Temperature))) / (1.0e6)
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

            end associate
        end do

    end function compute_chemistry

    ! operator(+) implementations ----------------------------------------------
    pure function integrand_add_integrand(lhs, rhs) result(opr)
        !! `+` operator.
        class(carbonate_chem_type), intent(in) :: lhs !< Left hand side.
        class(integrand_type), intent(in) :: rhs !< Right hand side.
        real(RK), allocatable :: opr(:) !< Operator result.
        integer :: iv, jp, ks ! indices into vector, points, species
        allocate(opr(lhs%size()))

        select type (rhs)
        class is (carbonate_chem_type)
            iv = 0
            do ks = 1, lhs%nscl
                do jp = 1, lhs%npts
                    iv = iv + 1
                    opr(iv) = lhs%species(jp, ks) + rhs%species(jp, ks)
                end do
            end do
        end select
    end function integrand_add_integrand

    pure function integrand_add_real(lhs, rhs) result(opr)
        !! `+ real_array` operator.
        class(carbonate_chem_type), intent(in) :: lhs !< Left hand side.
        real(RK), intent(in) :: rhs(1:) !< Right hand side.
        real(RK), allocatable :: opr(:) !< Operator result.
        integer :: iv, jp, ks ! indices into vector, points, species
        allocate(opr(lhs%size()))

        iv = 0
        do ks = 1, lhs%nscl
            do jp = 1, lhs%npts
                iv = iv + 1
                opr(iv) = lhs%species(jp, ks) + rhs(iv)
            end do
        end do
    end function integrand_add_real

    pure function real_add_integrand(lhs, rhs) result(opr)
        !! `real_array +` operator.
        real(RK), intent(in) :: lhs(1:) !< Left hand side.
        class(carbonate_chem_type), intent(in) :: rhs     !< Left hand side.
        real(RK), allocatable :: opr(:) !< Operator result.
        integer :: iv, jp, ks ! indices into vector, points, species
        allocate(opr(rhs%size()))

        iv = 0
        do ks = 1, rhs%nscl
            do jp = 1, rhs%npts
                iv = iv + 1
                opr(iv) = lhs(iv) + rhs%species(jp, ks)
            end do
        end do
    end function real_add_integrand

    ! *
    pure function integrand_multiply_integrand(lhs, rhs) result(opr)
        !! `*` operator.
        class(carbonate_chem_type), intent(in) :: lhs     !< Left hand side.
        class(integrand_type), intent(in) :: rhs     !< Right hand side.
        real(RK), allocatable :: opr(:) !< Operator result.
        integer :: iv, jp, ks ! indices into vector, points, species
        allocate(opr(lhs%size()))

        select type (rhs)
        class is (carbonate_chem_type)
            iv = 0
            do ks = 1, lhs%nscl
                do jp = 1, lhs%npts
                    iv = iv + 1
                    opr(iv) = lhs%species(jp, ks) * rhs%species(jp, ks)
                end do
            end do
        end select
    end function integrand_multiply_integrand

    pure function integrand_multiply_real(lhs, rhs) result(opr)
        !! `* real_array` operator.
        class(carbonate_chem_type), intent(in) :: lhs !< Left hand side.
        real(RK), intent(in) :: rhs(1:) !< Right hand side.
        real(RK), allocatable :: opr(:) !< Operator result.
        integer :: iv, jp, ks ! indices into vector, points, species
        allocate(opr(lhs%size()))

        iv = 0
        do ks = 1, lhs%nscl
            do jp = 1, lhs%npts
                iv = iv + 1
                opr(iv) = lhs%species(jp, ks) * rhs(iv)
            end do
        end do
    end function integrand_multiply_real

    pure function real_multiply_integrand(lhs, rhs) result(opr)
        !! `real_array *` operator.
        real(RK), intent(in) :: lhs(1:) !< Left hand side.
        class(carbonate_chem_type), intent(in) :: rhs     !< Left hand side.
        real(RK), allocatable :: opr(:) !< Operator result.
        integer :: iv, jp, ks ! indices into vector, points, species
        allocate(opr(rhs%size()))

        iv = 0
        do ks = 1, rhs%nscl
            do jp = 1, rhs%npts
                iv = iv + 1
                opr(iv) = lhs(iv) * rhs%species(jp, ks)
            end do
        end do
    end function real_multiply_integrand

    pure function integrand_multiply_real_scalar(lhs, rhs) result(opr)
        !! `* real_scalar` operator.
        class(carbonate_chem_type), intent(in) :: lhs !< Left hand side.
        real(RK), intent(in) :: rhs !< Right hand side.
        real(RK), allocatable :: opr(:) !< Operator result.
        integer :: iv, jp, ks ! indices into vector, points, species
        allocate(opr(lhs%size()))

        iv = 0
        do ks = 1, lhs%nscl
            do jp = 1, lhs%npts
                iv = iv + 1
                opr(iv) = lhs%species(jp, ks) * rhs
            end do
        end do
    end function integrand_multiply_real_scalar

    pure function real_scalar_multiply_integrand(lhs, rhs) result(opr)
        !! `real_scalar *` operator.
        real(RK), intent(in) :: lhs !< Left hand side.
        class(carbonate_chem_type), intent(in) :: rhs     !< Left hand side.
        real(RK), allocatable :: opr(:)  !< Operator result.
        integer :: iv, jp, ks ! indices into vector, points, species
        allocate(opr(rhs%size()))

        iv = 0
        do ks = 1, rhs%nscl
            do jp = 1, rhs%npts
                iv = iv + 1
                opr(iv) = lhs * rhs%species(jp, ks)
            end do
        end do
    end function real_scalar_multiply_integrand

    ! -
    pure function integrand_sub_integrand(lhs, rhs) result(opr)
        !! `-` operator.
        class(carbonate_chem_type), intent(in) :: lhs !< Left hand side.
        class(integrand_type), intent(in) :: rhs !< Right hand side.
        real(RK), allocatable :: opr(:)  !< Operator result.
        integer :: iv, jp, ks ! indices into vector, points, species
        allocate(opr(lhs%size()))

        select type (rhs)
        class is (carbonate_chem_type)
            iv = 0
            do ks = 1, lhs%nscl
                do jp = 1, lhs%npts
                    iv = iv + 1
                    opr(iv) = lhs%species(jp, ks) - rhs%species(jp, ks)
                end do
            end do
        end select
    end function integrand_sub_integrand

    pure function integrand_sub_real(lhs, rhs) result(opr)
        !! `- real_array` operator.
        class(carbonate_chem_type), intent(in) :: lhs !< Left hand side.
        real(RK), intent(in) :: rhs(1:) !< Right hand side.
        real(RK), allocatable :: opr(:) !< Operator result.
        integer :: iv, jp, ks ! indices into vector, points, species
        allocate(opr(lhs%size()))

        iv = 0
        do ks = 1, lhs%nscl
            do jp = 1, lhs%npts
                iv = iv + 1
                opr(iv) = lhs%species(jp, ks) - rhs(iv)
            end do
        end do
    end function integrand_sub_real

    pure function real_sub_integrand(lhs, rhs) result(opr)
        !! `real_array -` operator.
        real(RK), intent(in) :: lhs(1:) !< Left hand side.
        class(carbonate_chem_type), intent(in) :: rhs !< Left hand side.
        real(RK), allocatable :: opr(:)  !< Operator result.
        integer :: iv, jp, ks ! indices into vector, points, species
        allocate(opr(rhs%size()))

        iv = 0
        do ks = 1, rhs%nscl
            do jp = 1, rhs%npts
                iv = iv + 1
                opr(iv) = lhs(iv) - rhs%species(jp, ks)
            end do
        end do
    end function real_sub_integrand

    ! =
    subroutine integrand_eq_integrand(lhs, rhs)
        !< `=` operator.
        ! Intrinsic assignment would associate the species pointers, not the
        ! underlying data, but to use two objects inside an integrator with an
        ! update line like `y = y_new`, then we must force a copy of the species
        ! data from the RHS to the LHS of the assignment operator.
        class(carbonate_chem_type), intent(inout) :: lhs !< Left hand side.
        class(integrand_type), intent(in) :: rhs !< Right hand side.

        select type (rhs)
        class is (carbonate_chem_type)
            lhs%npts = rhs%npts
            lhs%args => rhs%args
            ! force the copy of species data, or nullify
            if (associated(rhs%species)) then
                if (.not. associated(lhs%species)) allocate(lhs%species(lhs%npts, lhs%nscl))
                lhs%species(:, :) = rhs%species
            else
                nullify (lhs%species)
            end if
        end select
    end subroutine integrand_eq_integrand

    pure subroutine integrand_eq_real(lhs, rhs)
        !< Assign a real to an integrand field.
        class(carbonate_chem_type), intent(inout) :: lhs !< Left hand side.
        real(RK), intent(in) :: rhs(1:) !< Right hand side.
        integer :: iv, jp, ks ! indices into vector, points, species

        iv = 0
        do ks = 1, lhs%nscl
            do jp = 1, lhs%npts
                iv = iv + 1
                lhs%species(jp, ks) = rhs(iv)
            end do
        end do
    end subroutine integrand_eq_real

end module chemistry
