module class_ode_marcher_fixed_step

   use constants_module
   use class_ode_system
   use class_ode_stepper
   use class_ode_marcher

   private

   type, public, extends(ode_marcher) :: ode_marcher_fixed_step
      real, pointer, contiguous :: y0(:)
      integer :: count = 0
   contains
      procedure :: init
      procedure :: apply
      procedure :: reset
      procedure :: free
   end type ode_marcher_fixed_step

contains

   subroutine init( self, error )
      class(ode_marcher_fixed_step), intent(inout) :: self
      integer, optional, intent(out) :: error
      !> local variables
      integer :: err, n
      err = FPDE_STATUS_OK
      n = self % dim

      if ( n .le. 0 ) then
         call self%loge("dimension less then zero")
         err = FPDE_STATUS_ERROR
         goto 01
      end if

      if ( .not. associated( self % s) ) then
         call self%loge("stepper not set")
         err = FPDE_STATUS_ERROR
         goto 01
      else
         call self%logd("stepper '"//trim(self%s%name)//"' is set")
      end if

      if ( err == FPDE_STATUS_OK ) then
         allocate( self%y0(n), stat=err )
         if( err == 0 ) then
            err = FPDE_STATUS_OK
            call self%logd("memory allocated")
         else
            err = FPDE_STATUS_ERROR
            call self%loge("memory allocation failed")
         end if
      end if
01    if ( present(error) ) error = err
   end subroutine init


   subroutine apply( self, sys, y, t, t1, h, error )
      class(ode_marcher_fixed_step), intent(inout) :: self
      class(ode_system) :: sys
      real, intent(inout) :: t
      real, intent(in) :: t1
      real, optional, intent(inout) :: h
      real, pointer, intent(inout) :: y(:)
      integer, optional, intent(out) :: error
      !> local variables
      integer :: err
      real :: t0, dt
      err = FPDE_STATUS_OK
      ! necessary copies
      t0=t
      dt=t1-t0
      ! check compatibility of dimensions
      if ( self % dim /= self % s % dim ) then
         err = FPDE_STATUS_ERROR
         call self%loge("different dimensions of stepper and marcher")
         goto 02
      end if
      ! check integration direction
      if ( (dt<0.0 .and. h>0.0) .or. (dt>0.0 .and. h<0.0) ) then
         call self%loge("marching direction does not match integration interval")
         err = FPDE_STATUS_ERROR
         goto 02
      end if
      ! backup input vector
      self % y0 = y
      ! perform one integration step
      call self%s%apply( sys=sys, y=y, t=t0, h=h, error=err )
      ! check stepper exit status
      if ( err /= FPDE_STATUS_OK ) then
         ! restore input vector
         y = self % y0
         call self%loge("stepper failed")
      else
         self % count = self % count + 1
      end if
      ! set error flag and return
02    if ( present(error) ) error = err
      return
   end subroutine apply


   subroutine reset( self, error )
      class(ode_marcher_fixed_step), intent(inout) :: self
      integer, optional, intent(out) :: error
      !> local variables
      integer :: err
      if ( present(error) ) error = err
   end subroutine reset


   subroutine free( self, error )
      class(ode_marcher_fixed_step), intent(inout) :: self
      integer, optional, intent(out) :: error
      !> local variables
      integer :: err
      deallocate(self%y0, stat=err)
      if (err == 0 ) then
         err = FPDE_STATUS_OK
      else
         err = FPDE_STATUS_ERROR
         call self%loge("memory deallocation failed")
      end if
      if(present(error)) error = err
   end subroutine free

end module class_ode_marcher_fixed_step
