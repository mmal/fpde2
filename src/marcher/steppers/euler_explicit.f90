module class_ode_stepper_euler_explicit

   use constants_module
   use class_ode_system
   use class_ode_stepper

   private

   type, public, extends(ode_stepper) :: ode_stepper_euler_explicit
      real, pointer, contiguous :: k(:)
   contains
      procedure :: init
      procedure :: apply
      procedure :: reset
      procedure :: free
   end type ode_stepper_euler_explicit

contains

   subroutine init( self, error )
      class(ode_stepper_euler_explicit), intent(inout) :: self
      integer, optional, intent(out) :: error
      !> Local variables
      integer :: err
      if ( self%dim <= 0 ) then
         err = FPDE_STATUS_ERROR
         call self%loge("incorrect stepper dimension")
         return
      end if
      allocate(self%k(self%dim), stat=err)
      if ( err == 0 ) then
         err = FPDE_STATUS_OK
      else
         err = FPDE_STATUS_ERROR
         call self%loge("memory allocation failed")
      end if
      if ( present(error) ) error = err
   end subroutine init


   subroutine apply( self, sys, y, t, h, yerr, dydt_in, dydt_out, error )
      class(ode_stepper_euler_explicit), intent(inout) :: self
      class(ode_system), intent(inout) :: sys
      real, intent(inout) :: y(:)
      real, intent(in) :: t
      real, intent(inout) :: h
      !> Optional arguments
      real, optional, intent(inout) :: yerr(:)
      real, optional, intent(in)  :: dydt_in(:)
      real, optional, intent(inout) :: dydt_out(:)
      integer, optional, intent(out) :: error
      !> Local variables
      integer :: err
      err = FPDE_STATUS_OK
      !> compute derivatives
      call sys%fun(t, y, self % k, sys % status)
      if ( sys % status /= FPDE_STATUS_OK ) then
         ! get the error status
         err = sys % status
         ! call the logger
         call self%loge("ode%fun failed")
         ! since input error is untouched we do not need
         ! to restore it
      else
         y = y + h * self % k !> apply Euler formula
      end if
      if ( present(error) ) error = err
   end subroutine apply


   subroutine reset( self, error )
      class(ode_stepper_euler_explicit), intent(inout) :: self
      integer, optional, intent(out) :: error
      !> Local variables
      integer :: err
      err = FPDE_STATUS_OK
      if ( associated(self % k) ) then
         self % k = 0.0
      else
         call self%loge("cannot reset when not initialized")
         err = FPDE_STATUS_ERROR
      end if
      if ( present(error) ) error = err
   end subroutine reset


   subroutine free( self, error )
      class(ode_stepper_euler_explicit) :: self
      integer, optional, intent(out) :: error
      !> Local variables
      integer :: err
      err = FPDE_STATUS_OK
      deallocate(self%k, stat=err)
      if (err /= 0 ) then
         err = FPDE_STATUS_ERROR
         call self%loge("memory deallocation failed")
      end if
      if(present(error)) error = err
   end subroutine free

end module class_ode_stepper_euler_explicit
