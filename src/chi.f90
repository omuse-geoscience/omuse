subroutine chi(tau,A_H,R_H,Lx,Ly,lambda0,lambda1,e111,phi1z0, &
             & Nm,Nx,Ny,dx,H,rho,beta0,err_tol,max_it,relax_coef, &
             & psi,free_slip,vis_bot_curr,vis_bot_prev,vis_lat_prev, &
             & chi_prev,chii)

implicit none
integer, intent(in) :: Nx,Ny,Nm,max_it,free_slip
real(8), intent(in) :: tau,A_H,R_H,H,rho,beta0,Lx,Ly,lambda0,lambda1,e111,phi1z0, &
                     & dx,err_tol,relax_coef
real(8), dimension(Nm,Nx,Ny), intent(in)  :: psi,chi_prev,vis_bot_curr,vis_bot_prev,vis_lat_prev
real(8), dimension(Nm,Nx,Ny), intent(out) :: chii

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
if( free_slip == 1) then
 rhs(1,:)   = 0.d0
 rhs(Nx,:)  = 0.d0
 rhs(:,1)   = 0.d0
 rhs(:,Ny)  = 0.d0
else if( free_slip == 0) then
 rhs(1,1)   = 0.d0
 rhs(Nx,1)  = 0.d0
 rhs(1,Ny)  = 0.d0
 rhs(Nx,Ny) = 0.d0
 do j=2,Ny-1
  rhs(1,j)  = 2.d0*chi_prev(1,1+1,j)
 end do
 do j=2,Ny-1
  rhs(Nx,j) = 2.d0*chi_prev(1,Nx-1,j)
 end do
 do i=2,Nx-1
  rhs(i,1)  = 2.d0*chi_prev(1,i,1+1)
 end do
 do i=2,Nx-1
  rhs(i,Ny) = 2.d0*chi_prev(1,i,Ny-1)
 end do
end if

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
if( free_slip == 1) then
 rhs(1,:)   = 0.d0
 rhs(Nx,:)  = 0.d0
 rhs(:,1)   = 0.d0
 rhs(:,Ny)  = 0.d0
else if( free_slip == 0) then
 rhs(1,1)   = 0.d0
 rhs(Nx,1)  = 0.d0
 rhs(1,Ny)  = 0.d0
 rhs(Nx,Ny) = 0.d0
 do j=2,Ny-1
  rhs(1,j)  = 2.d0*chi_prev(1,1+1,j)
 end do
 do j=2,Ny-1
  rhs(Nx,j) = 2.d0*chi_prev(1,Nx-1,j)
 end do
 do i=2,Nx-1
  rhs(i,1)  = 2.d0*chi_prev(1,i,1+1)
 end do
 do i=2,Nx-1
  rhs(i,Ny) = 2.d0*chi_prev(1,i,Ny-1)
 end do
end if

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
    & +0.25*d*d*d*d*jac_term01*(lambda1**2) &                                    !one extra d*d due to vis in jacobian.f 
    & -0.25*d*d*d*d*jac_term1LP1*e111 &                                          !one extra d*d due to vis in jacobian.f 
    & -R_H*d*d*vis_bot_prev(2,:,:) &
    & +A_H*d*d*d*d*vis_lat_prev(2,:,:)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! BOUNDARY VALUES -- necessary for the POISSON SOLVER !!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if( free_slip == 1) then
 rhs(1,:)   = 0.d0
 rhs(Nx,:)  = 0.d0
 rhs(:,1)   = 0.d0
 rhs(:,Ny)  = 0.d0
else if( free_slip == 0) then
 rhs(1,1)   = 0.d0
 rhs(Nx,1)  = 0.d0
 rhs(1,Ny)  = 0.d0
 rhs(Nx,Ny) = 0.d0
 do j=2,Ny-1
  rhs(1,j)  = 2.d0*chi_prev(2,1+1,j)
 end do
 do j=2,Ny-1
  rhs(Nx,j) = 2.d0*chi_prev(2,Nx-1,j)
 end do
 do i=2,Nx-1
  rhs(i,1)  = 2.d0*chi_prev(2,i,1+1)
 end do
 do i=2,Nx-1
  rhs(i,Ny) = 2.d0*chi_prev(2,i,Ny-1)
 end do
end if

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

