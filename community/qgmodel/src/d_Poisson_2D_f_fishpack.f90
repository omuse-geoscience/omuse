subroutine Poisson_2D_double_precision(Lx,Ly,Nnx,Nny,lambda,f,chii)

! Include modules defined by mkl_poisson.f90 and mkl_dfti.f90 header files
! use mkl_poisson

implicit none
integer, intent(in) :: Nnx,Nny
real(8), intent(in) :: Lx,Ly,lambda
real(8), dimension(Nnx,Nny), intent(inout)  :: f
real(8), dimension(Nnx,Nny), intent(out) :: chii

double precision pi
parameter(pi=3.14159265358979324D0)

integer nx, ny, stat
integer AA
DOUBLE PRECISION ptrb
double precision ax, bx, ay, by, q
double precision bd_ax(Nny), bd_bx(Nny), bd_ay(Nnx), bd_by(Nnx)
!double precision dpar(13*(Nnx-1)/2+7)
DOUBLE PRECISION, allocatable,dimension (:) :: workarea
!type(DFTI_DESCRIPTOR), pointer :: xhandle
character(4) BCtype

nx=Nnx-1
ny=Nny-1
AA = (1./log(2.))*log(1.0 * Nnx)
ALLOCATE(workarea(4*(Nnx) + (13 + nx)*(Nny)))
! Defining the rectangular domain 0<x<Lx, 0<y<Ly for 2D Poisson Solver
ax=0.0D0
bx=Lx
ay=0.0D0
by=Ly

!*******************************************************************************
! Setting the coefficient q to 0.
! Note that this is the way to use Helmholtz Solver to solve Poisson problem!
!*******************************************************************************
q=-(lambda*lambda)

!******************************************************************************
! NOTE THE SIGN CONVENTION USED FOR THE POISSON EQUATION!!!!!!!!!!!!!!!!!!!!!!!
!******************************************************************************

!******************************************************************************
! BOUNDARY CONDITIONS
!******************************************************************************
BCtype = 'DDDD'
bd_ax(:) = 0.0D0
bd_bx(:) = 0.0D0
bd_ay(:) = 0.0D0
bd_by(:) = 0.0D0

f(1,:) = -chii(1,:)
f(Nnx,:) = -chii(Nnx,:)
f(:,1) = -chii(:,1)
f(:,Nny) = -chii(:,Nny)
! NOTE:
!print *, ' The SOLUTION depends NOT on boundary values of the rhs for Dirichlet !!!!!!!!!!!!!!'
!f(1,:) = 1.d10
!write(*,*) f(1,:)
!print *, ' The SOLUTION depends on boundary values of the rhs for Neumann !!!!!!!!!!!!!!!!!!!!'
!f(:,1) = 1.d10
!write(*,*) f(:,1)
f = -1.0d0 * f

ptrb = 0.0d0
CALL HWSCRT(&
&    ax, bx, nx, 1, bd_ax, bd_bx, &
&    ay, by, ny, 1, bd_ay, bd_by, &
&    q, &
&    f, nx+1,&
&    ptrb, stat, workarea &
&)


!print *, f
DEALLOCATE(workarea)
! Now we can use xhandle to solve another 2D Poisson problem
chii = f

! Free MKL memory if any was allocated
!call mkl_free_buffers

end 


subroutine mkl_set_dynamic(input)
implicit none
integer, intent(in) :: input

end subroutine

subroutine mkl_set_num_threads(input)
implicit none
integer, intent(in) :: input

end subroutine

