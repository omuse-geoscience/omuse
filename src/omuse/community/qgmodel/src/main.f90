program main

implicit none
integer    :: i,counter,savecounter,save_num,max_it,free_slip, &
            & Nx,Ny,Nt,Nm,restart,restart_num
real(8)    :: tau,A_H,R_H,H,rho,beta0,Lx,Ly,lambda0,lambda1,e111,phi1z0, &
            & dx,dy,dt,T,t_curr,err_tol,relax_coef
real(8), allocatable,dimension (:,:)   :: max_psi
real(8), allocatable,dimension (:,:,:) :: psi_1,psi_2,inter,chii,chi_prev, &
                                        & vis_bot_prev,vis_bot_curr, &
                                        & vis_lat_prev,vis_lat_curr
character(len=25) :: filename

integer :: boundary(4)

call mkl_set_dynamic(0)
call mkl_set_num_threads( 16 )
 
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     PARAMETER SPECIFICATION
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!Read the physical parameters of the system:
open(25,file='parameter_physical.dat',status='old',action='read')
 read(25,*) tau
 read(25,*) A_H
 read(25,*) R_H
 read(25,*) lambda0
 read(25,*) lambda1
 read(25,*) e111
 read(25,*) phi1z0
 read(25,*) H
 read(25,*) rho
 read(25,*) beta0
 read(25,*) Nm
close(25)
!call parameter_physical(Re,delta_M,delta_S)
open(29,file='data/parameter_physical.dat',status='replace',action='write')
 write(29,*) 'tau     =' ,tau
 write(29,*) 'A_H     =' ,A_H
 write(29,*) 'R_H     =' ,R_H
 write(29,*) 'lambda0 =' ,lambda0
 write(29,*) 'lambda1 =' ,lambda1
 write(29,*) 'e111    =' ,e111
 write(29,*) 'phi1z0  =' ,phi1z0
 write(29,*) 'H       =' ,H
 write(29,*) 'rho     =' ,rho
 write(29,*) 'beta0   =' ,beta0
 write(29,*) 'Nm      =' ,Nm
close(29)

!Read the numerical parameters of the system:
open(25,file='parameter_numerical.dat',status='old',action='read')
 read(25,*) Lx
 read(25,*) Ly
 read(25,*) dx
 read(25,*) dy
 read(25,*) dt
 read(25,*) T
 read(25,*) savecounter
 read(25,*) err_tol
 read(25,*) max_it
 read(25,*) relax_coef
 read(25,*) free_slip
 boundary(1:4)=free_slip
 read(25,*) restart
 read(25,*) restart_num
close(25)
open(29,file='data/parameter_numerical.dat',status='replace',action='write')
 write(29,*) 'Lx = ' ,Lx
 write(29,*) 'Ly = ' ,Ly
 write(29,*) 'dx = ' ,dx
 write(29,*) 'dy = ' ,dy
 write(29,*) 'dt = ' ,dt
 write(29,*) 'T  = ' ,T
 write(29,*) 'savecounter  = ' ,savecounter
 write(29,*) 'err_tol      = ' ,err_tol
 write(29,*) 'max_it       = ' ,max_it
 write(29,*) 'relax_coef   = ' ,relax_coef
 write(29,*) 'free_slip    = ' ,free_slip                         ! 1 : free-slip, 0 : no-slip
 write(29,*) 'restart      = ' ,restart                           ! 1 : yes,       0 : no
 write(29,*) 'restart_num  = ' ,restart_num
close(29)

Nx = Lx/dx + 1
Ny = Ly/dy + 1
!Nt = T/dt + 1
Nt = INT(T/dt)/savecounter + 1

!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     INITIALIZE ARRAYS
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allocate (psi_1(Nm,Nx,Ny))
psi_1(:,:,:)        = 0.d0
allocate (psi_2(Nm,Nx,Ny))
psi_2(:,:,:)        = 0.d0
allocate (inter(Nm,Nx,Ny))
inter(:,:,:)        = 0.d0
allocate (chii(Nm,Nx,Ny))
chii(:,:,:)         = 0.d0
allocate (chi_prev(Nm,Nx,Ny))
chi_prev(:,:,:)     = 0.d0
allocate (vis_bot_prev(Nm,Nx,Ny))
vis_bot_prev(:,:,:) = 0.d0
allocate (vis_bot_curr(Nm,Nx,Ny))
vis_bot_curr(:,:,:) = 0.d0
allocate (vis_lat_prev(Nm,Nx,Ny))
vis_lat_prev(:,:,:) = 0.d0
allocate (vis_lat_curr(Nm,Nx,Ny))
vis_lat_curr(:,:,:) = 0.d0
allocate (max_psi(Nm,Nt))
max_psi(:,:) = 0.d0

if( restart==1 ) then
  write(filename,'(a,i5.5,a)')'psi1_1_',restart_num,'.dat'
  open(25,file='data/'//filename,status='old',action='read')
   do i=1,Nx
    read(25,'(1000000e25.16)',advance='yes') psi_1(1,i,:)
   end do
  close(25)
  write(filename,'(a,i5.5,a)')'psi2_1_',restart_num,'.dat'
  open(25,file='data/'//filename,status='old',action='read')
   do i=1,Nx
    read(25,'(1000000e25.16)',advance='yes') psi_2(1,i,:)
   end do
  close(25)
  if( Nm==2 ) then
   write(filename,'(a,i5.5,a)')'psi1_2_',restart_num,'.dat'
   open(25,file='data/'//filename,status='old',action='read')
    do i=1,Nx
     read(25,'(1000000e25.16)',advance='yes') psi_1(2,i,:)
    end do
   close(25)
   write(filename,'(a,i5.5,a)')'psi2_2_',restart_num,'.dat'
   open(25,file='data/'//filename,status='old',action='read')
    do i=1,Nx
     read(25,'(1000000e25.16)',advance='yes') psi_2(2,i,:)
    end do
   close(25)
  end if
 call vis_bot(Nm,Nx,Ny,boundary,psi_2,vis_bot_curr)
 call vis_lat(Nm,Nx,Ny,boundary,psi_2,vis_lat_curr)
end if
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
!     START TIME LOOP
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t_curr   = 0.d0
counter  = 0
save_num = 0

do while ( t_curr < T )

 t_curr = t_curr + dt

! update viscosity
 vis_bot_prev  = vis_bot_curr
 call vis_bot(Nm,Nx,Ny,boundary,psi_1,vis_bot_curr)
 vis_lat_prev  = vis_lat_curr
 call vis_lat(Nm,Nx,Ny,boundary,psi_1,vis_lat_curr)

!find chi, and take a step
 call chi(tau,A_H,R_H,Lx,Ly,lambda0,lambda1,e111,phi1z0, &
     &    Nm,Nx,Ny,dx,H,rho,beta0,err_tol,max_it,relax_coef, &
     &    psi_1,boundary,vis_bot_curr,vis_bot_prev,vis_lat_prev,chi_prev, &
     &    chii)
!Robert-Asselin filter
 inter = psi_2 + 2.d0*dt*chii
 psi_1 = psi_1 + 0.1d0*(inter-2.d0*psi_1+psi_2)
 psi_2 = inter
!!!
 chi_prev = chii
 counter = counter + 1
! max_psi(counter+save_num*savecounter+1) = maxval(psi_2)

! do exactly the same thing again with opposite psi arrays
! (specials: leap-frog time stepping, and using previous viscosity)
 t_curr = t_curr + dt

 vis_bot_prev  = vis_bot_curr
 call vis_bot(Nm,Nx,Ny,boundary,psi_2,vis_bot_curr)
 vis_lat_prev  = vis_lat_curr
 call vis_lat(Nm,Nx,Ny,boundary,psi_2,vis_lat_curr)
 call chi(tau,A_H,R_H,Lx,Ly,lambda0,lambda1,e111,phi1z0, &
     &    Nm,Nx,Ny,dx,H,rho,beta0,err_tol,max_it,relax_coef, &
     &    psi_2,boundary,vis_bot_curr,vis_bot_prev,vis_lat_prev,chi_prev, &
     &    chii)
!Robert-Asselin filter
 inter = psi_1 + 2.d0*dt*chii
 psi_2 = psi_2 + 0.1d0*(inter-2.d0*psi_2+psi_1)
 psi_1 = inter
!!!
 chi_prev = chii
 counter = counter + 1
! max_psi(counter+save_num*savecounter+1) = maxval(psi_1)

! write out psi to a file so it can later be plotted
 if(counter==savecounter) then
  save_num = save_num + 1

  write(filename,'(a,i5.5,a)')'psi1_1_',save_num,'.dat'
  open(29,file='data/'//filename,status='replace',action='write')
   do i=1,Nx
    write(29,'(1000000e25.16)',advance='yes') psi_1(1,i,:)
   end do
  close(29)
  write(filename,'(a,i5.5,a)')'psi2_1_',save_num,'.dat'
  open(29,file='data/'//filename,status='replace',action='write')
   do i=1,Nx
    write(29,'(1000000e25.16)',advance='yes') psi_2(1,i,:)
   end do
  close(29)
  max_psi(1,save_num+1) = maxval(psi_1(1,:,:))
  open(29,file='data/max_psi1.dat',status='replace',action='write')
    write(29,'(1000000000e25.16)',advance='yes') max_psi(1,:)
  close(29)

  if( Nm==2 ) then
   write(filename,'(a,i5.5,a)')'psi1_2_',save_num,'.dat'
   open(29,file='data/'//filename,status='replace',action='write')
    do i=1,Nx
     write(29,'(1000000e25.16)',advance='yes') psi_1(2,i,:)
    end do
   close(29)
   write(filename,'(a,i5.5,a)')'psi2_2_',save_num,'.dat'
   open(29,file='data/'//filename,status='replace',action='write')
    do i=1,Nx
     write(29,'(1000000e25.16)',advance='yes') psi_2(2,i,:)
    end do
   close(29)
   max_psi(2,save_num+1) = maxval(psi_1(2,:,:))
   open(29,file='data/max_psi2.dat',status='replace',action='write')
     write(29,'(1000000000e25.16)',advance='yes') max_psi(2,:)
   close(29)
  end if

  if( maxval(psi_1(1,:,:))==0.0 .AND. maxval(psi_1(2,:,:))==0.0 ) then
   open(29,file='error_message.dat',status='replace',action='write')
    write(29,*) 'Run crashed!'
   close(29)
   STOP
  end if

  counter = 0
 end if

end do

end program main

