!>
!! @file   ode_system.f90
!! @author Maciej Maliborski <maciej.maliborski@gmail.com>
!! @date   Fri Mar 23 19:36:27 2012
!!
!! @brief  ODE system class.
!!
module class_ode_system

   use constants_module
   use logger_module

   private

   public :: fun_interface, jac_interface

   abstract interface
      subroutine fun_interface( t, y, dydt, params, status )
         real, intent(in) :: t
         real, pointer, contiguous, intent(in) :: y(:)
         real, pointer, contiguous, intent(out) :: dydt(:)
         class(*) :: params
         integer, optional :: status
      end subroutine fun_interface

      subroutine jac_interface( t, y, dfdy, dfdt, params, status )
         real, intent(in) :: t
         real, pointer, contiguous, intent(in) :: y(:)
         real, pointer, contiguous, intent(out) :: dfdy(:,:)
         real, pointer, contiguous, intent(out) :: dfdt(:)
         class(*) :: params
         integer, optional :: status
      end subroutine jac_interface

      subroutine msm_interface( t, y, m, params, status )
         real, intent(in) :: t
         real, pointer, contiguous, intent(in) :: y(:)
         real, pointer, contiguous, intent(out) :: m(:)
         class(*) :: params
         integer, optional :: status
      end subroutine msm_interface

   end interface

   type, public, extends(named) :: ode_system

      !> ODE differential equation right hand side function pointer
      procedure(fun_interface), pointer, nopass :: fun => null()
      !> ODE differential equation jacobian function pointer
      procedure(jac_interface), pointer, nopass :: jac => null()
      !> ODE differential equation mass matrix function pointer
      procedure(msm_interface), pointer, nopass :: msm => null()
      !> dimension of the system (number of equations)
      integer :: dim = 0
      !> ODE parameters required to calculate the right side and jacobian
      !! of the given ODE
      class(*), pointer :: params
   contains
      !> Class ode_system initialization method
      procedure :: init
   end type ode_system

contains

   !> Class ode_system initialization method.
   !!
   !! @param sys ode_system object
   !! @param fun[in] right hand side of differential equation
   !! @param jac[in] right hand side of jacobian function
   !! @param dim[in] dimension of a system (number of equations)
   !! @param params[in] parameters required to calculate the
   !! right side and jacobian of the ODE
   !!
   subroutine init( sys, fun, jac, msm, dim, params )
      class(ode_system) :: sys
      procedure(fun_interface), pointer, intent(in) :: fun
      procedure(jac_interface), pointer, intent(in), optional :: jac
      procedure(msm_interface), pointer, intent(in), optional :: msm
      integer, intent(in) :: dim
      class(*), target, intent(in), optional :: params

      sys % dim = dim

      if ( sys % dim .le. 0 ) then
         call sys%log(FPDE_LOG_ERROR, &
              "ODE system dimension passed in sys%init cannot be <= 0")
      end if

      if( associated( fun ) ) then
         sys % fun => fun
      else
         call sys%log(FPDE_LOG_ERROR, &
              "ODE system no function passed in sys%init")
      end if

      if ( present( jac ) .and. associated( jac ) ) then
         sys % jac => jac
      end if

      if ( present( msm ) .and. associated( msm ) ) then
         sys % msm => msm
      end if

      if ( present( params ) ) then
         sys % params => params
      end if

      sys % status = FPDE_STATUS_OK

   end subroutine init

end module class_ode_system
