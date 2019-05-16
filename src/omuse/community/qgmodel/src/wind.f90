subroutine wind(Nm,Nx,Ny,windy)

implicit none
integer, intent(in) :: Nx,Ny,Nm
real(8), dimension(Nm,Nx,Ny), intent(out) :: windy

integer :: i,j,m
real(8) :: pi = 3.14159265358979d0
real(8), dimension(Nx,Ny) :: tau

tau(:,:)     = 0.d0
windy(:,:,:) = 0.d0

do i=1,Nx
 do j=1,Ny

  tau(i,j) = cos(2.*pi*((j-1.)/(Ny-1.)-0.5))+2.*sin(pi*((j-1.)/(Ny-1.)-0.5))
!  tau(i,j) = -1./(2.*pi)*cos(2.*pi*(j-1.)/(Ny-1.))
!  tau(i,j) = -cos(pi*(j-1.)/(Ny-1.))
!  tau(i,j) = -cos(pi*(j-1.5)/(ny-2.))                                             ! what they had in the code
!  tau(i,j) = -1./pi*cos(pi*(j-1.)/(Ny-1.))                                         ! what I used
!  tau(i,j) = -1./pi*cos(pi*(j-1.)/(ny-1.))*sin(pi*(i-1.)/(ny-1.))                 ! Veronis (1966) ????

 end do
end do

do m=1,Nm
 do i=2,Nx-1
  do j=2,Ny-1

  windy(m,i,j) = -tau(i,j+1)+tau(i,j-1)                                            ! include the whole curl, i.e., also x-derivative for generality ?????

  end do
 end do
end do

!      write(*,*) tau
!      write(*,*) windy

end


!      real(8), allocatable, dimension (:,:) :: wind_term
!      allocate (wind_term(Nx,Ny))
!      call wind(Nx,Ny,wind_term)


