module integrator
!< Define the abstract type [[integrator_type]] of FOODIE ODE integrators.

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
    public :: integrator_type

    integer, parameter :: IK = kind(1) ! default integer

    type, abstract :: integrator_type
    !< Abstract type for the numerical integration of a system of ODEs of variable size and shape
        integer(IK) :: error = 0 !< Error status code.
        character(len=:), allocatable :: error_message !< Error message, hopefully meaningful.
    contains
        ! public methods
        procedure, pass(self) :: check_error      !< Check for error occurrencies.
        procedure, pass(self) :: description      !< Return informative integrator description.
        procedure, pass(self) :: destroy_abstract !< Destroy only members of abstract [[integrator_type]] type.
        procedure, pass(self) :: trigger_error    !< Trigger an error.
        ! deferred methods
        procedure(assignment_interface), pass(lhs), deferred :: integr_assign_integr !< Operator `=`.
        ! operators
        generic :: assignment(=) => integr_assign_integr !< Overload `=`.
    end type integrator_type

    abstract interface
        subroutine assignment_interface(lhs, rhs)
            !< Operator `=`.
            import :: integrator_type
            class(integrator_type), intent(inout) :: lhs !< Left hand side.
            class(integrator_type), intent(in) :: rhs !< Right hand side.
        end subroutine assignment_interface
    end interface

contains

    ! public methods
    subroutine check_error(self, is_severe)
        !< Check for error occurencies.
        !<
        !< If `is_severe=.true.` an stop is called.
        class(integrator_type), intent(in) :: self      !< Integrator.
        logical, intent(in), optional :: is_severe !< Flag to activate severe failure, namely errors trigger a stop.

        if (self%error /= 0) then
            if (allocated(self%error_message)) then
                write (stderr, '(A)') 'error: '//self%error_message
            else
                write (stderr, '(A)') 'error: an obscure error occurred!'
            end if
            write (stderr, '(A,I4)') 'error code: ', self%error
            if (present(is_severe)) then
                if (is_severe) stop
            end if
        end if
    end subroutine check_error

    elemental subroutine destroy_abstract(self)
        !< Destroy only members of abstract [[integrator_type]] type.
        class(integrator_type), intent(inout) :: self !< Integrator.

        if (allocated(self%description_)) deallocate (self%description_)
        self%error = 0
        if (allocated(self%error_message)) deallocate (self%error_message)
    end subroutine destroy_abstract

    subroutine trigger_error(self, error, error_message, is_severe)
        !< Check for error occurencies.
        !<
        !< If `is_severe=.true.` an stop is called.
        class(integrator_type), intent(inout) :: self          !< Integrator.
        integer(IK), intent(in) :: error         !< Error status code.
        character(len=*), intent(in), optional :: error_message !< Error message, hopefully meaningful.
        logical, intent(in), optional :: is_severe     !< Flag to activate severe faliure, namely errors trigger a stop.

        self%error = error
        if (present(error_message)) self%error_message = error_message
        call self%check_error(is_severe=is_severe)
    end subroutine trigger_error

end module integrator
