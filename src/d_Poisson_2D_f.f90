subroutine Poisson_2D_double_precision(Lx,Ly,Nnx,Nny,lambda,f,chii)

! Include modules defined by mkl_poisson.f90 and mkl_dfti.f90 header files
use mkl_poisson

implicit none
integer, intent(in) :: Nnx,Nny
real(8), intent(in) :: Lx,Ly,lambda
real(8), dimension(Nnx,Nny), intent(in)  :: f
real(8), dimension(Nnx,Nny), intent(out) :: chii

double precision pi
parameter(pi=3.14159265358979324D0)

integer nx, ny, stat
integer ipar(128)
double precision ax, bx, ay, by, q
double precision bd_ax(Nny), bd_bx(Nny), bd_ay(Nnx), bd_by(Nnx)
double precision dpar(13*(Nnx-1)/2+7)
type(DFTI_DESCRIPTOR), pointer :: xhandle
character(4) BCtype

nx=Nnx-1
ny=Nny-1

! Defining the rectangular domain 0<x<Lx, 0<y<Ly for 2D Poisson Solver
ax=0.0D0
bx=Lx
ay=0.0D0
by=Ly

!*******************************************************************************
! Setting the coefficient q to 0.
! Note that this is the way to use Helmholtz Solver to solve Poisson problem!
!*******************************************************************************
q=(lambda*lambda)

!******************************************************************************
! NOTE THE SIGN CONVENTION USED FOR THE POISSON EQUATION!!!!!!!!!!!!!!!!!!!!!!!
!******************************************************************************

!******************************************************************************
! BOUNDARY CONDITIONS
!******************************************************************************
BCtype = 'DDDD'
!bd_ax(:) = 0.0D0
!bd_bx(:) = 0.0D0
!bd_ay(:) = 0.0D0
!bd_by(:) = 0.0D0

bd_ax(:) = chii(1,:)
bd_bx(:) = chii(Nnx,:)
bd_ay(:) = chii(:,1)
bd_by(:) = chii(:,Nny)
! NOTE:
!print *, ' The SOLUTION depends NOT on boundary values of the rhs for Dirichlet !!!!!!!!!!!!!!'
!f(1,:) = 1.d10
!write(*,*) f(1,:)
!print *, ' The SOLUTION depends on boundary values of the rhs for Neumann !!!!!!!!!!!!!!!!!!!!'
!f(:,1) = 1.d10
!write(*,*) f(:,1)


! Initializing ipar array to make it free from garbage
ipar(:)=0

! Initializing simple data structures of Poisson Library for 2D Poisson Solver
call d_init_Helmholtz_2D(ax, bx, ay, by, nx, ny, BCtype, q, ipar, dpar, stat)
if (stat.ne.0) stop 111

! Initializing complex data structures of Poisson Library for 2D Poisson Solver
! NOTE: Right-hand side f may be altered after the Commit step. If you want
! to keep it, you should save it in another memory location!
call d_commit_Helmholtz_2D(f, bd_ax, bd_bx, bd_ay, bd_by, xhandle, ipar, dpar, stat)
if (stat.ne.0) stop 222

! Computing the approximate solution of 2D Poisson problem
! NOTE: Boundary data stored in the arrays bd_ax, bd_bx, bd_ay, bd_by
! should not be changed between the Commit step and the subsequent call to
! the Solver routine! Otherwise the results may be wrong.
call d_Helmholtz_2D(f, bd_ax, bd_bx, bd_ay, bd_by, xhandle, ipar, dpar, stat)
if (stat.ne.0) stop 333

! Cleaning the memory used by xhandle
call free_Helmholtz_2D(xhandle, ipar, stat)
if (stat.ne.0) stop 444
! Now we can use xhandle to solve another 2D Poisson problem

chii = f

! Free MKL memory if any was allocated
call mkl_free_buffers

end 

