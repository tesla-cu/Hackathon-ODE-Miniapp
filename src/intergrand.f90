module integrand
!< Define the abstract type `integrand` for solving ordinary differential equations of variable size and shape

!< Copyright (C) 2022, Stefano Zhagi, et al.
!< <https://github.com/Fortran-FOSS-Programmers/FOODIE>

!< This program is free software: you can redistribute it and/or modify
!< it under the terms of the GNU General Public License as published by
!< the Free Software Foundation, either version 3 of the License, or
!< (at your option) any later version.

!< This program is distributed in the hope that it will be useful,
!< but WITHOUT ANY WARRANTY; without even the implied warranty of
!< MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!< GNU General Public License for more details.

!< You should have received a copy of the GNU General Public License
!< along with this program.  If not, see <http://www.gnu.org/licenses/>.

    implicit none
    private
    public :: integrand_type, IK, RK

    integer, parameter :: IK = kind(1), RK = kind(1.0) ! default integer and real for now

    type, abstract :: integrand_type
    private
    !< Abstract type for a system of ODEs of variable size and shape
    contains
        ! public deferred procedures that concrete integrand-field must implement
        procedure(integrand_dimension), pass(self), deferred, public :: size !< Return integrand size.
        procedure(integrand_state), pass(self), deferred, public :: state !< Return integrand state array.
        procedure(time_derivative), pass(self), deferred, public :: d_dt !< Time derivative, residuals.

        ! operators
        ! procedure(local_error_operator), pass(lhs), deferred :: local_error !< `||integrand - integrand||` operator.
        ! generic, public :: operator(.lterror.) => local_error !< Estimate local truncation error.
        ! +
        procedure(symmetric_operator), pass(lhs), deferred :: integrand_add_integrand !< `+` operator.
        procedure(integrand_op_real), pass(lhs), deferred :: integrand_add_real      !< `+ real` operator.
        procedure(real_op_integrand), pass(rhs), deferred :: real_add_integrand      !< `real +` operator.
        generic, public :: operator(+) => &
            integrand_add_integrand, integrand_add_real, real_add_integrand !< Overloading `+` operator.
        ! *
        procedure(symmetric_operator), pass(lhs), deferred :: integrand_multiply_integrand   !< `*` operator.
        procedure(integrand_op_real), pass(lhs), deferred :: integrand_multiply_real        !< `* real` operator.
        procedure(real_op_integrand), pass(rhs), deferred :: real_multiply_integrand        !< `real *` operator.
        procedure(integrand_op_real_scalar), pass(lhs), deferred :: integrand_multiply_real_scalar !< `* real_scalar` operator.
        procedure(real_scalar_op_integrand), pass(rhs), deferred :: real_scalar_multiply_integrand !< `real_scalar *` operator.
        generic, public :: operator(*) => &
            integrand_multiply_integrand, &
            integrand_multiply_real, &
            real_multiply_integrand, &
            integrand_multiply_real_scalar, &
            real_scalar_multiply_integrand !< Overloading `*` operator.
        ! -
        procedure(symmetric_operator), pass(lhs), deferred :: integrand_sub_integrand !< `-` operator.
        procedure(integrand_op_real), pass(lhs), deferred :: integrand_sub_real      !< `- real` operator.
        procedure(real_op_integrand), pass(rhs), deferred :: real_sub_integrand      !< `real -` operator.
        generic, public :: operator(-) => &
            integrand_sub_integrand, integrand_sub_real, real_sub_integrand !< Overloading `-` operator.
        ! =
        procedure(symmetric_assignment), pass(lhs), deferred :: integrand_eq_integrand !< `=` operator.
        procedure(assignment_from_real), pass(lhs), deferred :: integrand_eq_real      !< `= real` operator.
        generic, public :: assignment(=) => &
            integrand_eq_integrand, integrand_eq_real !< Overloading `=` assignament.
    end type integrand_type

    abstract interface
        !< Abstract type bound procedures necessary for implementing a concrete extension of [[integrand_type]].

        pure function integrand_dimension(self)
            !< Return integrand dimension.
            import :: integrand_type, IK
            class(integrand_type), intent(in) :: self !< Integrand.
            integer(IK) :: integrand_dimension !< Integrand dimension.
        end function integrand_dimension

        function integrand_state(self)
            !< Return integrand dimension.
            import :: integrand_type, RK
            class(integrand_type), intent(in) :: self !< Integrand.
            real(RK), contiguous, pointer :: integrand_state(:) !< Integrand state.
        end function integrand_state

        function time_derivative(self, t) result(dState_dt)
            !< Time derivative function of integrand class, i.e. the residuals function.
            import :: integrand_type, RK
            class(integrand_type), intent(in) :: self !< Integrand field.
            real(RK), intent(in), optional :: t !< Time.
            real(RK), allocatable :: dState_dt(:) !< Result of the time derivative function of integrand field.
        end function time_derivative

        ! operators
        function local_error_operator(lhs, rhs) result(error)
            !< Estimate local truncation error between 2 solution approximations.
            import :: integrand_type, RK
            class(integrand_type), intent(in) :: lhs   !< Left hand side.
            class(integrand_type), intent(in) :: rhs   !< Right hand side.
            real(RK) :: error !< Error estimation.
        end function local_error_operator

        pure function integrand_op_real(lhs, rhs) result(operator_result)
            !< Asymmetric type operator `integrand.op.real`.
            import :: integrand_type, RK
            class(integrand_type), intent(in) :: lhs                !< Left hand side.
            real(RK), intent(in) :: rhs(1:)            !< Right hand side.
            real(RK), allocatable :: operator_result(:) !< Operator result.
        end function integrand_op_real

        pure function real_op_integrand(lhs, rhs) result(operator_result)
            !< Asymmetric type operator `real.op.integrand`.
            import :: integrand_type, RK
            class(integrand_type), intent(in) :: rhs                !< Right hand side.
            real(RK), intent(in) :: lhs(1:)            !< Left hand side.
            real(RK), allocatable :: operator_result(:) !< Operator result.
        end function real_op_integrand

        pure function integrand_op_real_scalar(lhs, rhs) result(operator_result)
            !< Asymmetric type operator `integrand.op.real`.
            import :: integrand_type, RK
            class(integrand_type), intent(in) :: lhs    !< Left hand side.
            real(RK), intent(in) :: rhs                 !< Right hand side.
            real(RK), allocatable :: operator_result(:) !< Operator result.
        end function integrand_op_real_scalar

        pure function real_scalar_op_integrand(lhs, rhs) result(operator_result)
            !< Asymmetric type operator `real.op.integrand`.
            import :: integrand_type, RK
            real(RK), intent(in) :: lhs                !< Left hand side.
            class(integrand_type), intent(in) :: rhs                !< Right hand side.
            real(RK), allocatable :: operator_result(:) !< Operator result.
        end function real_scalar_op_integrand

        pure function symmetric_operator(lhs, rhs) result(operator_result)
            !< Symmetric type operator integrand.op.integrand.
            import :: integrand_type, RK
            class(integrand_type), intent(in) :: lhs                !< Left hand side.
            class(integrand_type), intent(in) :: rhs                !< Right hand side.
            real(RK), allocatable :: operator_result(:) !< Operator result.
        end function symmetric_operator

        subroutine symmetric_assignment(lhs, rhs)
            !< Symmetric assignment integrand = integrand.
            import :: integrand_type
            class(integrand_type), intent(inout) :: lhs !< Left hand side.
            class(integrand_type), intent(in) :: rhs !< Right hand side.
        end subroutine symmetric_assignment

        pure subroutine assignment_from_real(lhs, rhs)
            !< Asymmetric assignment integrand = real.
            import :: integrand_type, RK
            class(integrand_type), intent(inout) :: lhs     !< Left hand side.
            real(RK), intent(in) :: rhs(1:) !< Right hand side.
        end subroutine assignment_from_real
    end interface

end module integrand
