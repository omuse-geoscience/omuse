subroutine vis_lat(Nm,Nx,Ny,boundary,psi,vis_latt,boundaries)

implicit none

type boundary_grid
  real(8), allocatable,dimension (:,:,:) :: psi
  real(8), allocatable,dimension (:,:,:) :: chi
end type

integer, intent(in) :: boundary(4)
integer, intent(in) :: Nx,Ny,Nm
real(8), dimension(Nm,Nx,Ny), intent(in)  :: psi
real(8), dimension(Nm,Nx,Ny), intent(out) :: vis_latt
type(boundary_grid) :: boundaries(4)
integer :: fac(4)

integer :: i,j,m

vis_latt(:,:,:) = 0.d0

do m=1,Nm
 do i=3,Nx-2
  do j=3,Ny-2

   vis_latt(m,i,j) = psi(m,i+2,j)+psi(m,i-2,j)+psi(m,i,j+2)+psi(m,i,j-2)+ &
                   &2*psi(m,i+1,j+1)+2*psi(m,i+1,j-1)+2*psi(m,i-1,j+1)+2*psi(m,i-1,j-1)- &
                   &8*psi(m,i+1,j)-8*psi(m,i-1,j)-8*psi(m,i,j+1)-8*psi(m,i,j-1)+20*psi(m,i,j)

  end do
 end do
end do

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! IMPLEMENT BOUNDARY CONDITIONS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fac=0
do i=1,4
  if(boundary(i)==1) fac(i)=-1
  if(boundary(i)==0) fac(i)=1
enddo

do m=1,Nm

  i=2
  do j=3,Ny-2
  
   vis_latt(m,i,j) = psi(m,i+2,j)+fac(1)*psi(m,2,j)+psi(m,i,j+2)+psi(m,i,j-2)+ &
                   & 2*psi(m,i+1,j+1)+2*psi(m,i+1,j-1)+2*psi(m,i-1,j+1)+2*psi(m,i-1,j-1)- &
                   & 8*psi(m,i+1,j)-8*psi(m,i-1,j)-8*psi(m,i,j+1)-8*psi(m,i,j-1)+20*psi(m,i,j)
  
  end do
  
  i=Nx-1
  do j=3,Ny-2
  
   vis_latt(m,i,j) =fac(2)*psi(m,Nx-1,j)+psi(m,i-2,j)+psi(m,i,j+2)+psi(m,i,j-2)+ &
                   & 2*psi(m,i+1,j+1)+2*psi(m,i+1,j-1)+2*psi(m,i-1,j+1)+2*psi(m,i-1,j-1)- &
                   & 8*psi(m,i+1,j)-8*psi(m,i-1,j)-8*psi(m,i,j+1)-8*psi(m,i,j-1)+20*psi(m,i,j)
  
  end do
  
  j=2
  do i=3,Nx-2
  
   vis_latt(m,i,j) = psi(m,i+2,j)+psi(m,i-2,j)+psi(m,i,j+2)+fac(3)*psi(m,i,2)+ &
                   & 2*psi(m,i+1,j+1)+2*psi(m,i+1,j-1)+2*psi(m,i-1,j+1)+2*psi(m,i-1,j-1)- &
                   & 8*psi(m,i+1,j)-8*psi(m,i-1,j)-8*psi(m,i,j+1)-8*psi(m,i,j-1)+20*psi(m,i,j)
  
  end do
  
  j=Ny-1
  do i=3,Nx-2
  
   vis_latt(m,i,j) = psi(m,i+2,j)+psi(m,i-2,j)+fac(4)*psi(m,i,Ny-1)+psi(m,i,j-2)+ &
                   & 2*psi(m,i+1,j+1)+2*psi(m,i+1,j-1)+2*psi(m,i-1,j+1)+2*psi(m,i-1,j-1)- &
                   & 8*psi(m,i+1,j)-8*psi(m,i-1,j)-8*psi(m,i,j+1)-8*psi(m,i,j-1)+20*psi(m,i,j)
  
  end do
  
  i=2
  j=2
  
   vis_latt(m,i,j) = psi(m,i+2,j)+fac(1)*psi(m,2,j)+psi(m,i,j+2)+fac(3)*psi(m,i,2)+ &
                   & 2*psi(m,i+1,j+1)+2*psi(m,i+1,j-1)+2*psi(m,i-1,j+1)+2*psi(m,i-1,j-1)- &
                   & 8*psi(m,i+1,j)-8*psi(m,i-1,j)-8*psi(m,i,j+1)-8*psi(m,i,j-1)+20*psi(m,i,j)
  
  i=Nx-1
  j=2
  
   vis_latt(m,i,j) =fac(2)*psi(m,Nx-1,j)+psi(m,i-2,j)+psi(m,i,j+2)+fac(3)*psi(m,i,2)+ &
                   & 2*psi(m,i+1,j+1)+2*psi(m,i+1,j-1)+2*psi(m,i-1,j+1)+2*psi(m,i-1,j-1)- &
                   & 8*psi(m,i+1,j)-8*psi(m,i-1,j)-8*psi(m,i,j+1)-8*psi(m,i,j-1)+20*psi(m,i,j)
  
  i=2
  j=Ny-1
  
   vis_latt(m,i,j) = psi(m,i+2,j)+fac(1)*psi(m,2,j)+fac(4)*psi(m,i,Ny-1)+psi(m,i,j-2)+ &
                   & 2*psi(m,i+1,j+1)+2*psi(m,i+1,j-1)+2*psi(m,i-1,j+1)+2*psi(m,i-1,j-1)- &
                   & 8*psi(m,i+1,j)-8*psi(m,i-1,j)-8*psi(m,i,j+1)-8*psi(m,i,j-1)+20*psi(m,i,j)
  
  i=Nx-1
  j=Ny-1
  
   vis_latt(m,i,j) =fac(2)*psi(m,Nx-1,j)+psi(m,i-2,j)+fac(4)*psi(m,i,Ny-1)+psi(m,i,j-2)+ &
                   & 2*psi(m,i+1,j+1)+2*psi(m,i+1,j-1)+2*psi(m,i-1,j+1)+2*psi(m,i-1,j-1)- &
                   & 8*psi(m,i+1,j)-8*psi(m,i-1,j)-8*psi(m,i,j+1)-8*psi(m,i,j-1)+20*psi(m,i,j)

end do

! fix interface boundaries
! add missing boundary terms

if(boundary(1)==2) then
  do m=1,Nm
    i=2
    do j=2,Ny-1    
      vis_latt(m,i,j) = vis_latt(m,i,j)+boundaries(1)%psi(m,i-2,j)
    enddo
  enddo
endif

if(boundary(2)==2) then
  do m=1,Nm
    i=Nx-1
    do j=2,Ny-1    
      vis_latt(m,i,j) = vis_latt(m,i,j)+boundaries(2)%psi(m,i+2,j)
    enddo
  enddo
endif

if(boundary(3)==2) then
  do m=1,Nm
    j=2
    do i=2,Nx-1    
      vis_latt(m,i,j) = vis_latt(m,i,j)+boundaries(3)%psi(m,i,j-2)
    enddo
  enddo
endif

if(boundary(4)==2) then
  do m=1,Nm
    j=Ny-1
    do i=2,Nx-1    
      vis_latt(m,i,j) = vis_latt(m,i,j)+boundaries(4)%psi(m,i,j+2)
    enddo
  enddo
endif

end


