subroutine vis_bot(Nm,Nx,Ny,free_slip,psi,vis_bott)

implicit none
integer, intent(in) :: Nx,Ny,Nm,free_slip
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
if( free_slip == 1 ) then

 ! nothing to do

else if( free_slip == 0 ) then

 do m=1,Nm
  i=1
  do j=2,Ny-1
   vis_bott(m,i,j) = 2.d0*psi(m,i+1,j)
  end do

  i=Nx
  do j=2,Ny-1
   vis_bott(m,i,j) = 2.d0*psi(m,i-1,j)
  end do

  j=1
  do i=2,Nx-1
   vis_bott(m,i,j) = 2.d0*psi(m,i,j+1)
  end do

  j=Ny
  do i=2,Nx-1
   vis_bott(m,i,j) = 2.d0*psi(m,i,j-1)
  end do

 end do

end if

end

