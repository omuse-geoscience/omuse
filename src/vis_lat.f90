subroutine vis_lat(Nm,Nx,Ny,free_slip,psi,vis_latt)

implicit none
integer, intent(in) :: Nx,Ny,Nm,free_slip
real(8), dimension(Nm,Nx,Ny), intent(in)  :: psi
real(8), dimension(Nm,Nx,Ny), intent(out) :: vis_latt

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
if( free_slip == 1) then

 do m=1,Nm
 i=2
 do j=3,Ny-2

   vis_latt(m,i,j) = psi(m,i+2,j)-psi(m,2,j)+psi(m,i,j+2)+psi(m,i,j-2)+ &
                   & 2*psi(m,i+1,j+1)+2*psi(m,i+1,j-1)+2*psi(m,i-1,j+1)+2*psi(m,i-1,j-1)- &
                   & 8*psi(m,i+1,j)-8*psi(m,i-1,j)-8*psi(m,i,j+1)-8*psi(m,i,j-1)+20*psi(m,i,j)

 end do

 i=Nx-1
 do j=3,Ny-2

   vis_latt(m,i,j) =-psi(m,Nx-1,j)+psi(m,i-2,j)+psi(m,i,j+2)+psi(m,i,j-2)+ &
                   & 2*psi(m,i+1,j+1)+2*psi(m,i+1,j-1)+2*psi(m,i-1,j+1)+2*psi(m,i-1,j-1)- &
                   & 8*psi(m,i+1,j)-8*psi(m,i-1,j)-8*psi(m,i,j+1)-8*psi(m,i,j-1)+20*psi(m,i,j)

 end do

 j=2
 do i=3,Nx-2

   vis_latt(m,i,j) = psi(m,i+2,j)+psi(m,i-2,j)+psi(m,i,j+2)-psi(m,i,2)+ &
                   & 2*psi(m,i+1,j+1)+2*psi(m,i+1,j-1)+2*psi(m,i-1,j+1)+2*psi(m,i-1,j-1)- &
                   & 8*psi(m,i+1,j)-8*psi(m,i-1,j)-8*psi(m,i,j+1)-8*psi(m,i,j-1)+20*psi(m,i,j)

 end do

 j=Ny-1
 do i=3,Nx-2

   vis_latt(m,i,j) = psi(m,i+2,j)+psi(m,i-2,j)-psi(m,i,Ny-1)+psi(m,i,j-2)+ &
                   & 2*psi(m,i+1,j+1)+2*psi(m,i+1,j-1)+2*psi(m,i-1,j+1)+2*psi(m,i-1,j-1)- &
                   & 8*psi(m,i+1,j)-8*psi(m,i-1,j)-8*psi(m,i,j+1)-8*psi(m,i,j-1)+20*psi(m,i,j)

 end do

 i=2
 j=2

   vis_latt(m,i,j) = psi(m,i+2,j)-psi(m,2,j)+psi(m,i,j+2)-psi(m,i,2)+ &
                   & 2*psi(m,i+1,j+1)+2*psi(m,i+1,j-1)+2*psi(m,i-1,j+1)+2*psi(m,i-1,j-1)- &
                   & 8*psi(m,i+1,j)-8*psi(m,i-1,j)-8*psi(m,i,j+1)-8*psi(m,i,j-1)+20*psi(m,i,j)

 i=Nx-1
 j=2

   vis_latt(m,i,j) =-psi(m,Nx-1,j)+psi(m,i-2,j)+psi(m,i,j+2)-psi(m,i,2)+ &
                   & 2*psi(m,i+1,j+1)+2*psi(m,i+1,j-1)+2*psi(m,i-1,j+1)+2*psi(m,i-1,j-1)- &
                   & 8*psi(m,i+1,j)-8*psi(m,i-1,j)-8*psi(m,i,j+1)-8*psi(m,i,j-1)+20*psi(m,i,j)

 i=2
 j=Ny-1

   vis_latt(m,i,j) = psi(m,i+2,j)-psi(m,2,j)-psi(m,i,Ny-1)+psi(m,i,j-2)+ &
                   & 2*psi(m,i+1,j+1)+2*psi(m,i+1,j-1)+2*psi(m,i-1,j+1)+2*psi(m,i-1,j-1)- &
                   & 8*psi(m,i+1,j)-8*psi(m,i-1,j)-8*psi(m,i,j+1)-8*psi(m,i,j-1)+20*psi(m,i,j)

 i=Nx-1
 j=Ny-1

   vis_latt(m,i,j) =-psi(m,Nx-1,j)+psi(m,i-2,j)-psi(m,i,Ny-1)+psi(m,i,j-2)+ &
                   & 2*psi(m,i+1,j+1)+2*psi(m,i+1,j-1)+2*psi(m,i-1,j+1)+2*psi(m,i-1,j-1)- &
                   & 8*psi(m,i+1,j)-8*psi(m,i-1,j)-8*psi(m,i,j+1)-8*psi(m,i,j-1)+20*psi(m,i,j)

 end do

else if( free_slip ==0) then

 do m=1,Nm
 i=2
 do j=3,Ny-2

   vis_latt(m,i,j) = psi(m,i+2,j)+psi(m,2,j)+psi(m,i,j+2)+psi(m,i,j-2)+ &
                   & 2*psi(m,i+1,j+1)+2*psi(m,i+1,j-1)+2*psi(m,i-1,j+1)+2*psi(m,i-1,j-1)- &
                   & 8*psi(m,i+1,j)-8*psi(m,i-1,j)-8*psi(m,i,j+1)-8*psi(m,i,j-1)+20*psi(m,i,j)

 end do

 i=Nx-1
 do j=3,Ny-2

   vis_latt(m,i,j) = psi(m,Nx-1,j)+psi(m,i-2,j)+psi(m,i,j+2)+psi(m,i,j-2)+ &
                   & 2*psi(m,i+1,j+1)+2*psi(m,i+1,j-1)+2*psi(m,i-1,j+1)+2*psi(m,i-1,j-1)- &
                   & 8*psi(m,i+1,j)-8*psi(m,i-1,j)-8*psi(m,i,j+1)-8*psi(m,i,j-1)+20*psi(m,i,j)
       
 end do

 j=2
 do i=3,Nx-2

   vis_latt(m,i,j) = psi(m,i+2,j)+psi(m,i-2,j)+psi(m,i,j+2)+psi(m,i,2)+ &
                   & 2*psi(m,i+1,j+1)+2*psi(m,i+1,j-1)+2*psi(m,i-1,j+1)+2*psi(m,i-1,j-1)- &
                   & 8*psi(m,i+1,j)-8*psi(m,i-1,j)-8*psi(m,i,j+1)-8*psi(m,i,j-1)+20*psi(m,i,j)

 end do

 j=Ny-1
 do i=3,Nx-2

   vis_latt(m,i,j) = psi(m,i+2,j)+psi(m,i-2,j)+psi(m,i,Ny-1)+psi(m,i,j-2)+ &
                   & 2*psi(m,i+1,j+1)+2*psi(m,i+1,j-1)+2*psi(m,i-1,j+1)+2*psi(m,i-1,j-1)- &
                   & 8*psi(m,i+1,j)-8*psi(m,i-1,j)-8*psi(m,i,j+1)-8*psi(m,i,j-1)+20*psi(m,i,j)

 end do

 i=2
 j=2

   vis_latt(m,i,j) = psi(m,i+2,j)+psi(m,2,j)+psi(m,i,j+2)+psi(m,i,2)+ &
                   & 2*psi(m,i+1,j+1)+2*psi(m,i+1,j-1)+2*psi(m,i-1,j+1)+2*psi(m,i-1,j-1)- &
                   & 8*psi(m,i+1,j)-8*psi(m,i-1,j)-8*psi(m,i,j+1)-8*psi(m,i,j-1)+20*psi(m,i,j)

 i=Nx-1
 j=2

   vis_latt(m,i,j) = psi(m,Nx-1,j)+psi(m,i-2,j)+psi(m,i,j+2)+psi(m,i,2)+ &
                   & 2*psi(m,i+1,j+1)+2*psi(m,i+1,j-1)+2*psi(m,i-1,j+1)+2*psi(m,i-1,j-1)- &
                   & 8*psi(m,i+1,j)-8*psi(m,i-1,j)-8*psi(m,i,j+1)-8*psi(m,i,j-1)+20*psi(m,i,j)

 i=2
 j=Ny-1

   vis_latt(m,i,j) = psi(m,i+2,j)+psi(m,2,j)+psi(m,i,Ny-1)+psi(m,i,j-2)+ &
                   & 2*psi(m,i+1,j+1)+2*psi(m,i+1,j-1)+2*psi(m,i-1,j+1)+2*psi(m,i-1,j-1)- &
                   & 8*psi(m,i+1,j)-8*psi(m,i-1,j)-8*psi(m,i,j+1)-8*psi(m,i,j-1)+20*psi(m,i,j)

 i=Nx-1
 j=Ny-1

   vis_latt(m,i,j) = psi(m,Nx-1,j)+psi(m,i-2,j)+psi(m,i,Ny-1)+psi(m,i,j-2)+ &
                   & 2*psi(m,i+1,j+1)+2*psi(m,i+1,j-1)+2*psi(m,i-1,j+1)+2*psi(m,i-1,j-1)- &
                   & 8*psi(m,i+1,j)-8*psi(m,i-1,j)-8*psi(m,i,j+1)-8*psi(m,i,j-1)+20*psi(m,i,j)

 end do

else

  write(*,*) 'Problem with boundary conditions'
  stop 999

end if

!write(*,*) vis_lat

end


