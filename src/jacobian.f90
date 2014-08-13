subroutine jacobian(Nx,Ny,psi,vis,jac_term)

implicit none
integer, intent(in) :: Nx,Ny
real(8), dimension(Nx,Ny), intent(in)  :: psi,vis
real(8), dimension(Nx,Ny), intent(out) :: jac_term

integer :: i,j

jac_term(:,:) = 0.d0

do i=2,Nx-1
 do j=2,Ny-1

!     Arakawa Jacobian
  jac_term(i,j) = ((psi(i+1,j)-psi(i-1,j))*(vis(i,j+1)-vis(i,j-1)) &
                & -(psi(i,j+1)-psi(i,j-1))*(vis(i+1,j)-vis(i-1,j)) &
                & +psi(i+1,j)*(vis(i+1,j+1)-vis(i+1,j-1)) &
                & -psi(i-1,j)*(vis(i-1,j+1)-vis(i-1,j-1)) &
                & -psi(i,j+1)*(vis(i+1,j+1)-vis(i-1,j+1)) &
                & +psi(i,j-1)*(vis(i+1,j-1)-vis(i-1,j-1)) &
                & +vis(i,j+1)*(psi(i+1,j+1)-psi(i-1,j+1)) &
                & -vis(i,j-1)*(psi(i+1,j-1)-psi(i-1,j-1)) &
                & -vis(i+1,j)*(psi(i+1,j+1)-psi(i+1,j-1)) &
                & +vis(i-1,j)*(psi(i-1,j+1)-psi(i-1,j-1))) &
                & *0.33333333

 end do
end do

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! NOTE: For THIS Jacobian the BOUNDARY VALUES vanish everywhere,
!
! for BOTH free-slip and no-slip boundary conditions!!!!!!!!!!!!!!!!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! we ignore the interface boundaries here for the moment
! (vis and psi know about the boundaries ofcourse)
! since the boundary values of the jacobian are of no consequence atm
! otherwise we would need to add vis on the boundary as boundary condition!!

end


