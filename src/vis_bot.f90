subroutine vis_bot(Nm,Nx,Ny,boundary,psi,vis_bott)

implicit none
integer, intent(in) :: boundary(4)
integer, intent(in) :: Nx,Ny,Nm
real(8), dimension(Nm,Nx,Ny), intent(in)  :: psi
real(8), dimension(Nm,Nx,Ny), intent(out) :: vis_bott

integer :: i,j,m

vis_bott(:,:,:) = 0.d0

do m=1,Nm
 do i=2,Nx-1
  do j=2,Ny-1

   vis_bott(m,i,j) = psi(m,i+1,j)+psi(m,i-1,j) &
                   &+psi(m,i,j+1)+psi(m,i,j-1)-4.d0*psi(m,i,j)

  end do
 end do
end do

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! BOUNDARY VALUES -- necessary for the JACOBIAN !!!
! but note: if you only have first or second order derivatives (e.g., the Laplacian only),
! then free slip is your default boundary condition!!!!!!!!!!!!!!!!!!!!
! (Because you cannot prescribe the second velocity component.)
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if( boundary(1) == 1) then
 ! nothing to do  
else if( boundary(1) == 0) then
  do m=1,Nm
   vis_bott(m,1,2:Ny-1) = 2.d0*psi(m,2,2:Ny-1)
  end do
else if(boundary(1) == 2) then
!tbd
endif

if( boundary(2) == 1) then
 ! nothing to do
else if( boundary(2) == 0) then
  do m=1,Nm
   vis_bott(m,Nx,2:Ny-1) = 2.d0*psi(m,Nx-1,2:Ny-1)
  end do
else if(boundary(2) == 2) then
! tbd
endif

if( boundary(3) == 1) then
 ! nothing to do
else if( boundary(3) == 0) then
  do m=1,Nm
   vis_bott(m,2:Nx-1,1) = 2.d0*psi(m,2:Nx-1,2)
  end do
else if(boundary(3) == 2) then
! tbd
endif

if( boundary(4) == 1) then
 ! nothing to do
else if( boundary(4) == 0) then
  do m=1,Nm
   vis_bott(m,2:Nx-1,Ny) = 2.d0*psi(m,2:Nx-1,Ny-1)
  end do
else if(boundary(4) == 2) then
! tbd
endif

end

