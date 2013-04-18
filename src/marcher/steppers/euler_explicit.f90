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

   subroutine init( p, error )
      class(ode_stepper_euler_explicit), intent(inout) :: p
      integer, optional, intent(out) :: error
      !> Local variables
      integer :: err
      associate(this => p)
        allocate(this%k(this%dim), stat=err)
        if (err == 0 ) then
           err = FPDE_STATUS_OK
        else
           err = FPDE_STATUS_ERROR
           call this%loge("memory allocation failed")
        end if
      end associate
      if ( present(error) ) error = err
   end subroutine init


   subroutine apply( this, sys, y, t, h, yerr, dydt_in, dydt_out, error )
      class(ode_stepper_euler_explicit), intent(inout) :: this
      class(ode_system), intent(inout) :: sys
      real, pointer, contiguous, intent(inout) :: y(:)
      real, intent(in) :: t
      real, intent(inout) :: h
      !> Optional arguments
      real, optional, pointer, contiguous, intent(inout) :: yerr(:)
      real, optional, pointer, contiguous, intent(in)  :: dydt_in(:)
      real, optional, pointer, contiguous, intent(inout) :: dydt_out(:)
      integer, optional, intent(out) :: error
      !> Local variables
      integer :: err
      err = FPDE_STATUS_OK
      !> compute derivatives
      call sys%fun(t, y, this % k, sys % params, sys % status)
      if ( sys % status /= FPDE_STATUS_OK ) then
         ! get the error status
         err = sys % status
         ! call the logger
         call this%loge("ode%fun failed")
         ! since input error is untouched we do not need
         ! to restore it
      else
         y = y + h * this % k !> apply Euler formula
      end if
      if ( present(error) ) error = err
   end subroutine apply


   subroutine reset( this, error )
      class(ode_stepper_euler_explicit), intent(inout) :: this
      integer, optional, intent(out) :: error
      !> Local variables
      integer :: err
      err = FPDE_STATUS_OK
      if ( associated(this % k) ) then
         this % k = 0.0
      else
         call this%loge("cannot reset when not initialized")
         err = FPDE_STATUS_ERROR
      end if
      if ( present(error) ) error = err
   end subroutine reset


   subroutine free( p, error )
      class(ode_stepper_euler_explicit) :: p
      integer, optional, intent(out) :: error
      !> Local variables
      integer :: err
      associate( this => p )
        deallocate(this%k, stat=err)
        if (err == 0 ) then
           err = FPDE_STATUS_OK
        else
           err = FPDE_STATUS_ERROR
           call this%loge("memory deallocation failed")
        end if
      end associate
      if(present(error)) error = FPDE_STATUS_OK
   end subroutine free

end module class_ode_stepper_euler_explicit
