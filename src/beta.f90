subroutine beta(Nm,Nx,Ny,psi,beta_term)

implicit none
integer, intent(in) :: Nx,Ny,Nm
real(8), dimension(Nm,Nx,Ny), intent(in)  :: psi
real(8), dimension(Nm,Nx,Ny), intent(out) :: beta_term

integer :: i,j,m

beta_term(:,:,:) = 0.d0

do m=1,Nm
 do i=2,Nx-1
  do j=2,Ny-1

   beta_term(m,i,j) = psi(m,i+1,j)-psi(m,i-1,j)

  end do
 end do
end do
!write(*,*) beta_term

end

