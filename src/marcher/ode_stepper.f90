!>
!! @file   ode_stepper.f90
!! @author Maciej Maliborski <maciej.maliborski@gmail.com>
!! @date   Fri Mar 23 19:25:11 2012
!!
!! @brief  ODE stepper class.
!!
module class_ode_stepper

   use constants_module
   use class_platonic
   use class_ode_system

   private

   !> Generic stepper object
   !!
   type, public, abstract, extends(platonic) :: ode_stepper
      !> Dimension of ODE system to solve (implies the amount of workspace)
      integer :: dim = 0
   contains
      !> Applies the stepping function to the system of equations, defined
      !! by sys, using the step-size h and advances the system from time t
      !! and state y(t) to time t+h and state y(t+h).
      procedure(apply), deferred :: apply
      !> Resets the stepper object (cleans workspace but not frees).
      !! It should be used whenever the next stepper use will not be
      !! a continuation of a previous step.
      procedure(reset), deferred :: reset
   end type ode_stepper

   interface

      subroutine apply( this, sys, y, t, h, yerr, dydt_in, dydt_out, error )
         import :: ode_stepper, ode_system
         class(ode_stepper), intent(inout) :: this
         class(ode_system), intent(inout) :: sys
         real, pointer, contiguous, intent(inout) :: y(:)
         real, intent(in) :: t
         real, intent(inout) :: h
         !> Optional arguments
         real, optional, pointer, contiguous, intent(inout) :: yerr(:)
         real, optional, pointer, contiguous, intent(in)  :: dydt_in(:)
         real, optional, pointer, contiguous, intent(inout) :: dydt_out(:)
         integer, optional, intent(out) :: error
      end subroutine apply

      subroutine reset( this, error )
         import :: ode_stepper
         class(ode_stepper), intent(inout) :: this
         integer, optional, intent(out) :: error
      end subroutine reset

   end interface

end module class_ode_stepper
