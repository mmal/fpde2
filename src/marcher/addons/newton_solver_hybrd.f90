module class_newton_solver_hybrd

   use constants_module
   use class_platonic
   use class_newton_solver
   use fenl_hybrd_interface

   private

   type, public, extends(newton_solver) :: newton_solver_hybrd
      real, allocatable, dimension(:) :: wa1, wa2, wa3, wa4, r, &
           & qtf, diag, fjac
   contains
      procedure :: init
      procedure :: solve
      procedure :: free
   end type newton_solver_hybrd

contains

   subroutine init( self, error )
      class(newton_solver_hybrd), intent(inout) :: self
      integer, optional, intent(out) :: error
      !> local variables
      integer :: n, err
      err = FPDE_STATUS_OK
      self%name = "hybrd newton solver"
      if ( self%dim <= 0 ) then
         err = FPDE_STATUS_ERROR
         call self%loge("incorrect dimension")
         return
      end if
      n = self%dim
      allocate( self%wa1(n), self%wa2(n), self%wa3(n), self%wa4(n), &
           & self%r(n*(n+1)/2), self%qtf(n), self%diag(n), self%fjac(n*n), stat=err )
      if ( err /= 0 ) then
         err = FPDE_STATUS_ERROR
         call self%loge("memory allocation failed")
      end if
      if( present(error) ) error = err
   end subroutine init


   subroutine solve( self, f, xvec, fvec, error )
      use fenl_hybrd_interface
      class(newton_solver_hybrd), intent(inout) :: self
      real, intent(inout) :: xvec(:)
      real, intent(out) :: fvec(:)
      integer, optional, intent(out) :: error
      procedure(f_) :: f
      !> local variables
      integer :: n, maxfev, ml, mu, mode, nprint, lr, index, info, nfev
      real :: xtol, epsfcn, factor
      !> user specified parameters
      n = self%dim
      maxfev = self%maxfev
      xtol = self%tol
      !>      ml(mu) is a nonnegative integer input variable which specifies
      !!        the number of subdiagonals (superdiagonals) within the band of the
      !!        jacobian matrix. if the jacobian is not banded, set
      !!        ml(mu) to at least n - 1.
      ml = n-1
      mu = n-1
      !>      epsfcn is an input variable used in determining a suitable
      !!        step length for the forward-difference approximation. this
      !!        approximation assumes that the relative errors in the
      !!        functions are of the order of epsfcn. if epsfcn is less
      !!        than the machine precision, it is assumed that the relative
      !!        errors in the functions are of the order of the machine
      !!        precision.
      epsfcn = 0.0
      !>      mode is an integer input variable. if mode = 1, the
      !!        variables will be scaled internally. if mode = 2,
      !!        the scaling is specified by the input diag. other
      !!        values of mode are equivalent to mode = 1.
      mode = 2
      self%diag = 1.0
      !> less important arguments of hybrd
      info = 0
      nprint = 0
      lr = (n*(n + 1))/2
      index = 6*n + lr
      factor = 100.0
      !> call the hybrd solver
      call hybrd(fcn,n,xvec,fvec,xtol,maxfev,ml,mu,epsfcn,self%diag,mode, &
           &   factor,nprint,info,nfev,self%fjac,n,self%r,lr, &
           &   self%qtf,self%wa1,self%wa2,self%wa3,self%wa4)
      if (info == 1 ) then
         call self%logd("algorithm estimates that the relative error is at most [tol]")
         info = FPDE_STATUS_OK
      else
         select case (info)
         case(0)
            call self%loge("improper input parameters")
         case(2)
            call self%loge("number of calls to [fcn] has reached or exceeded [maxfev]")
         case(3)
            call self%loge("[tol] is too small, no further improvement in the approximate solution [x] is possible")
         case(4:5)
            call self%loge("iteration is not making good progress")
         end select
         info = FPDE_STATUS_ERROR
      end if
      if( present(error) ) error = info
   contains
      subroutine fcn(n, x, fvec, iflag)
         integer, intent(in) :: n
         real, intent(in) :: x(n)
         real, intent(out) :: fvec(n)
         integer, intent(inout) :: iflag
         call f(x=x, y=fvec, error=iflag)
      end subroutine fcn
   end subroutine solve


   subroutine free( self, error )
      class(newton_solver_hybrd), intent(inout) :: self
      integer, optional, intent(out) :: error
      !> local variables
      integer :: err
      err = FPDE_STATUS_OK
      deallocate( self%wa1, self%wa2, self%wa3, self%wa4, &
           & self%r, self%qtf, self%diag, self%fjac, stat=err )
      if ( err /= 0 ) then
         err = FPDE_STATUS_ERROR
         call self%loge("memory allocation failed")
      end if
      if( present(error) ) error = err
   end subroutine free

end module class_newton_solver_hybrd
