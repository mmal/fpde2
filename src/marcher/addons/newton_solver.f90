module class_newton_solver

   use class_platonic

   private

   public :: f_

   type, public, abstract, extends(platonic) :: newton_solver
      integer :: dim, maxfev
      real :: tol
   contains
      procedure(solve), deferred :: solve
   end type newton_solver

   interface
      subroutine f_(x, y, error)
         real, intent(in) :: x(:)
         real, intent(out) :: y(:)
         integer, optional, intent(out) :: error
      end subroutine f_

      subroutine solve( self, f, xvec, fvec, error )
         import :: newton_solver, f_
         class(newton_solver), intent(inout) :: self
         real, intent(inout) :: xvec(:)
         real, intent(out) :: fvec(:)
         integer, optional, intent(out) :: error
         procedure(f_) :: f
      end subroutine solve
   end interface

end module class_newton_solver
