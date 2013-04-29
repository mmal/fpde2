!>
!! @file   ode_system.f90
!! @author Maciej Maliborski <maciej.maliborski@gmail.com>
!! @date   Fri Mar 23 19:36:27 2012
!!
!! @brief  ODE system class.
!!
module class_ode_system

   use logger_module

   private

   public :: jac_interface, msm_interface

   !> Generic ode system object
   !!
   type, public, abstract, extends(named) :: ode_system
      !> ODE differential equation jacobian function pointer
      procedure(jac_interface), pointer, pass :: jac => null()
      !> ODE differential equation mass matrix function pointer
      procedure(msm_interface), pointer, pass :: msm => null()
      !> dimension of the system (number of equations)
      integer :: dim = 0
   contains
      !> ODE differential equation right hand side function pointer
      procedure(fun_interface), deferred :: fun
   end type ode_system

   abstract interface

      subroutine fun_interface( self, t, y, dydt, error )
         import :: ode_system
         class(ode_system), intent(inout) :: self
         real, intent(in) :: t
         real, intent(in) :: y(:)
         real, intent(out) :: dydt(:)
         integer, optional :: error
      end subroutine fun_interface

      !> Pointer methods interfaces
      !!
      subroutine jac_interface( self, t, y, dfdy, dfdt, error )
         import :: ode_system
         class(ode_system) :: self
         real, intent(in) :: t
         real, intent(in) :: y(:)
         real, intent(out) :: dfdy(:,:)
         real, intent(out) :: dfdt(:)
         integer, optional :: error
      end subroutine jac_interface

      subroutine msm_interface( self, t, y, m, error )
         import :: ode_system
         class(ode_system) :: self
         real, intent(in) :: t
         real, intent(in) :: y(:)
         real, intent(out) :: m(:,:)
         integer, optional :: error
      end subroutine msm_interface

   end interface

end module class_ode_system
