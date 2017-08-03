subroutine chi(tau,A_H,R_H,Lx,Ly,lambda0,lambda1,e111,phi1z0, &
             & Nm,Nx,Ny,dx,H,rho,beta0,err_tol,max_it,relax_coef, &
             & psi,boundary,vis_bot_curr,vis_bot_prev,vis_lat_prev, &
             & chi_prev,chii, boundaries)

implicit none

type boundary_grid
  real(8), allocatable,dimension (:,:,:) :: psi
  real(8), allocatable,dimension (:,:,:) :: chi
end type

integer, intent(in) :: boundary(4)
integer, intent(in) :: Nx,Ny,Nm,max_it
real(8), intent(in) :: tau,A_H,R_H,H,rho,beta0,Lx,Ly,lambda0,lambda1,e111,phi1z0, &
                     & dx,err_tol,relax_coef
real(8), dimension(Nm,Nx,Ny), intent(in)  :: psi,chi_prev,vis_bot_curr,vis_bot_prev,vis_lat_prev
real(8), dimension(Nm,Nx,Ny), intent(out) :: chii
type(boundary_grid) :: boundaries(4)

integer :: m,i,j
real(8) :: d
real(8), dimension(Nx,Ny)    :: rhs,jac_term01,jac_term0LP0,jac_term1LP1,jac_term0LP1,jac_term1LP0
real(8), dimension(Nm,Nx,Ny) :: wind_term,beta_term

call wind(Nm,Nx,Ny,wind_term)
call beta(Nm,Nx,Ny,psi,beta_term)

d = 1./dx

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if( Nm==1) then
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
call jacobian(Nx,Ny,psi(1,:,:),vis_bot_curr(1,:,:),jac_term0LP0)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!      calculate right hand side
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rhs(:,:) = 0.d0
rhs = -beta0*0.5*d*beta_term(1,:,:) &
    & +(tau/(rho*H))*0.5*d*wind_term(1,:,:) &
    & -0.25*d*d*d*d*jac_term0LP0 &                                            ! one extra d*d due to vis in jacobian.f 
    & -R_H*d*d*vis_bot_prev(1,:,:) &
    & +A_H*d*d*d*d*vis_lat_prev(1,:,:)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! BOUNDARY VALUES -- necessary for the POISSON SOLVER !!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chii(1,1,:)=0
chii(1,Nx,:)=0
chii(1,:,1)=0
chii(1,:,Ny)=0

if( boundary(1) == 1) then
  rhs(1,:)   = 0.d0
else if( boundary(1) == 0) then
  rhs(1,1)   = 0.d0
  rhs(1,Ny)   = 0.d0
  rhs(1,2:Ny-1)  = 2.d0*chi_prev(1,2,2:Ny-1)
else if(boundary(1) == 2) then
!tbd: fix beta_term and wind_term and jac_term and possibly vis terms (atm not necessary)
  chii(1,1,1:Ny)=boundaries(1)%chi(1,1,1:Ny)
endif

if( boundary(2) == 1) then
  rhs(Nx,:)  = 0.d0
else if( boundary(2) == 0) then
  rhs(Nx,1)  = 0.d0
  rhs(Nx,Ny) = 0.d0
  rhs(Nx,2:Ny-1) = 2.d0*chi_prev(1,Nx-1,2:Ny-1)
else if(boundary(2) == 2) then
! tbd
  chii(1,Nx,1:Ny)=boundaries(2)%chi(1,Nx,1:Ny)
endif

if( boundary(3) == 1) then
  rhs(:,1)   = 0.d0
else if( boundary(3) == 0) then
  rhs(1,1)   = 0.d0
  rhs(Nx,1)  = 0.d0
  rhs(2:Nx-1,1)  = 2.d0*chi_prev(1,2:Nx-1,2)
else if(boundary(3) == 2) then
! tbd
  chii(1,1:Nx,1)=boundaries(3)%chi(1,1:Nx,1)
endif

if( boundary(4) == 1) then
  rhs(:,Ny)  = 0.d0
else if( boundary(4) == 0) then
  rhs(1,Ny)  = 0.d0
  rhs(Nx,Ny) = 0.d0
  rhs(2:Nx-1,Ny) = 2.d0*chi_prev(1,2:Nx-1,Ny-1)
else if(boundary(4) == 2) then
! tbd
  chii(1,1:Nx,Ny)=boundaries(4)%chi(1,1:Nx,Ny)
endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!      solve the poisson equation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rhs = -rhs
call Poisson_2D_double_precision(Lx,Ly,Nx,Ny,lambda0,rhs,chii(1,:,:))
!******************************************************************************
! NOTE THE SIGN CONVENTION USED FOR THE POISSON EQUATION!!!!!!!!!!!!!!!!!!!!!!!
!******************************************************************************

else if( Nm==2) then
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
call jacobian(Nx,Ny,psi(1,:,:),vis_bot_curr(1,:,:),jac_term0LP0)
call jacobian(Nx,Ny,psi(2,:,:),vis_bot_curr(2,:,:),jac_term1LP1)
call jacobian(Nx,Ny,psi(2,:,:),vis_bot_curr(1,:,:),jac_term1LP0)
call jacobian(Nx,Ny,psi(1,:,:),vis_bot_curr(2,:,:),jac_term0LP1)
call jacobian(Nx,Ny,psi(1,:,:),psi(2,:,:),jac_term01)

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!      calculate right hand side of barotropic mode
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rhs(:,:) = 0.d0
rhs = -beta0*0.5*d*beta_term(1,:,:) &
    & +(tau/(rho*H))*0.5*d*wind_term(1,:,:) &
    & -0.25*d*d*d*d*jac_term0LP0 &                                               !one extra d*d due to vis in jacobian.f 
    & -0.25*d*d*d*d*jac_term1LP1 &                                               !one extra d*d due to vis in jacobian.f 
    & -R_H*d*d*vis_bot_prev(1,:,:) &
    & +A_H*d*d*d*d*vis_lat_prev(1,:,:)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! BOUNDARY VALUES -- necessary for the POISSON SOLVER !!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chii(1,1,:)=0
chii(1,Nx,:)=0
chii(1,:,1)=0
chii(1,:,Ny)=0

if( boundary(1) == 1) then
  rhs(1,:)   = 0.d0
else if( boundary(1) == 0) then
  rhs(1,1)   = 0.d0
  rhs(1,Ny)   = 0.d0
  rhs(1,2:Ny-1)  = 2.d0*chi_prev(1,2,2:Ny-1)
else if(boundary(1) == 2) then
!tbd
  chii(1,1,1:Nx)=boundaries(1)%chi(1,1,1:Nx)
endif

if( boundary(2) == 1) then
  rhs(Nx,:)  = 0.d0
else if( boundary(2) == 0) then
  rhs(Nx,1)  = 0.d0
  rhs(Nx,Ny) = 0.d0
  rhs(Nx,2:Ny-1) = 2.d0*chi_prev(1,Nx-1,2:Ny-1)
else if(boundary(2) == 2) then
! tbd
  chii(1,Nx,1:Ny)=boundaries(2)%chi(1,Nx,1:Ny)
endif

if( boundary(3) == 1) then
  rhs(:,1)   = 0.d0
else if( boundary(3) == 0) then
  rhs(1,1)   = 0.d0
  rhs(Nx,1)  = 0.d0
  rhs(2:Nx-1,1)  = 2.d0*chi_prev(1,2:Nx-1,2)
else if(boundary(3) == 2) then
! tbd
  chii(1,1:Nx,1)=boundaries(3)%chi(1,1:Nx,1)
endif

if( boundary(4) == 1) then
  rhs(:,Ny)  = 0.d0
else if( boundary(4) == 0) then
  rhs(1,Ny)  = 0.d0
  rhs(Nx,Ny) = 0.d0
  rhs(2:Nx-1,Ny) = 2.d0*chi_prev(1,2:Nx-1,Ny-1)
else if(boundary(4) == 2) then
! tbd
  chii(1,1:Nx,Ny)=boundaries(4)%chi(1,1:Nx,Ny)
endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!      solve the helmholtz equation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rhs = -rhs
call Poisson_2D_double_precision(Lx,Ly,Nx,Ny,lambda0,rhs,chii(1,:,:))
!******************************************************************************
! NOTE THE SIGN CONVENTION USED FOR THE POISSON EQUATION!!!!!!!!!!!!!!!!!!!!!!!
!******************************************************************************

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!      calculate right hand side of baroclinic mode
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rhs(:,:) = 0.d0
rhs = -beta0*0.5*d*beta_term(2,:,:) &
    & +(tau/(rho*H))*0.5*d*wind_term(2,:,:)*phi1z0 &
    & -0.25*d*d*d*d*jac_term1LP0 &                                               !one extra d*d due to vis in jacobian.f 
    & -0.25*d*d*d*d*jac_term0LP1 &                                               !one extra d*d due to vis in jacobian.f 
    & +0.25*d*d*jac_term01*(lambda1**2) &                                    !one extra d*d due to vis in jacobian.f 
    & -0.25*d*d*d*d*jac_term1LP1*e111 &                                          !one extra d*d due to vis in jacobian.f 
    & -R_H*d*d*vis_bot_prev(2,:,:) &
    & +A_H*d*d*d*d*vis_lat_prev(2,:,:)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! BOUNDARY VALUES -- necessary for the POISSON SOLVER !!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
chii(2,1,:)=0
chii(2,Nx,:)=0
chii(2,:,1)=0
chii(2,:,Ny)=0

if( boundary(1) == 1) then
  rhs(1,:)   = 0.d0
else if( boundary(1) == 0) then
  rhs(1,1)   = 0.d0
  rhs(1,Ny)   = 0.d0
  rhs(1,2:Ny-1)  = 2.d0*chi_prev(2,2,2:Ny-1)
else if(boundary(1) == 2) then
!tbd
  chii(2,1,1:Nx)=boundaries(1)%chi(2,1,1:Nx)
endif

if( boundary(2) == 1) then
  rhs(Nx,:)  = 0.d0
else if( boundary(2) == 0) then
  rhs(Nx,1)  = 0.d0
  rhs(Nx,Ny) = 0.d0
  rhs(Nx,2:Ny-1) = 2.d0*chi_prev(2,Nx-1,2:Ny-1)
else if(boundary(2) == 2) then
! tbd
  chii(2,Nx,1:Ny)=boundaries(2)%chi(2,Nx,1:Ny)
endif

if( boundary(3) == 1) then
  rhs(:,1)   = 0.d0
else if( boundary(3) == 0) then
  rhs(1,1)   = 0.d0
  rhs(Nx,1)  = 0.d0
  rhs(2:Nx-1,1)  = 2.d0*chi_prev(2,2:Nx-1,2)
else if(boundary(3) == 2) then
! tbd
  chii(2,1:Nx,1)=boundaries(3)%chi(2,1:Nx,1)
endif

if( boundary(4) == 1) then
  rhs(:,Ny)  = 0.d0
else if( boundary(4) == 0) then
  rhs(1,Ny)  = 0.d0
  rhs(Nx,Ny) = 0.d0
  rhs(2:Nx-1,Ny) = 2.d0*chi_prev(2,2:Nx-1,Ny-1)
else if(boundary(4) == 2) then
! tbd
  chii(2,1:Nx,Ny)=boundaries(4)%chi(2,1:Nx,Ny)
endif

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!      solve the helmholtz equation
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rhs = -rhs
call Poisson_2D_double_precision(Lx,Ly,Nx,Ny,lambda1,rhs,chii(2,:,:))
!******************************************************************************
! NOTE THE SIGN CONVENTION USED FOR THE POISSON EQUATION!!!!!!!!!!!!!!!!!!!!!!!
!******************************************************************************
end if

end

