module qgmodel

  implicit none
  integer    :: i,counter,savecounter,save_num,max_it, &
              & Nx,Ny,Nt,Nm !,restart,restart_num
  real(8)    :: tau,A_H,R_H,H,rho,beta0,Lx,Ly,lambda0,lambda1,e111,phi1z0, &
              & dx,dy,dt,T,t_curr,err_tol,relax_coef
  real(8)    :: begin_time, wind_sigma, ra_alpha, raw_alpha
  real(8), allocatable,dimension (:,:)   :: max_psi
  real(8), allocatable,dimension (:,:,:) :: psi_1,psi_2,inter,chii,chi_prev, &
                                          & vis_bot_prev,vis_bot_curr, &
                                          & vis_lat_prev,vis_lat_curr
  real(8), allocatable,dimension(:,:) :: tau_x
  integer :: interface_wind
  type boundary_grid
    real(8), allocatable,dimension (:,:,:) :: psi
    real(8), allocatable,dimension (:,:,:) :: chi
  end type
  
  integer,parameter :: nbc=2
  type(boundary_grid) :: boundaries(4)
  
  character(len=25) :: filename
  integer :: boundary(4)

  character(len=15) :: timestep_method="leapfrog"

contains

function initialize_code() result(ret)
  integer :: ret
  ret=0
  
  call mkl_set_dynamic(0)
  call mkl_set_num_threads( 4 )

  call default_numerical_parameters()
  call default_physical_parameters()
    
end function

subroutine default_numerical_parameters()
  Lx             = 4000000.d0
  Ly             = 4000000.d0
  dx             = 10000.d0
  dy             = 10000.d0
  dt             = 1800.d0
  T              = 86400.d0
  savecounter    = 48
  err_tol        = 1.0d-6
  max_it         = 20000
  relax_coef     = 1.7d0
  boundary(1:4)  = 1   ! 1=free_slip, 0=no_slip, 2=interface
  
  ra_alpha       = 0.1d0
  raw_alpha      = 1.d0
  interface_wind = 0  ! 0=use wind routine, 1= use interface wind
!  restart        = 0  ! not used for amuse
!  restart_num    = 360 ! not used                
end subroutine  

subroutine default_physical_parameters()
  begin_time = 0.0
  wind_sigma = -99.0
  
  tau        = 0.05d0
  A_H        = 100.d0
  R_H        = 0.d0
  lambda0    = 0.d0
  lambda1    = 2.0000d-05
  e111       = 0.d0
  phi1z0     = 1.4142135623731d0
  H          = 4000.d0
  rho        = 1000.d0
  beta0      = 1.8616d-11
  Nm         = 1
end subroutine

subroutine initialize_arrays() 
  integer :: i
! properly should have some check whether allocations succeeded and return 
! the non-zero ret value

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

! specialty: indexing of boundary and normal grid the same!
  allocate(boundaries(1)%psi(Nm,0:nbc-1,0:Ny+nbc-1), &
           boundaries(1)%chi(Nm,0:nbc-1,0:Ny+nbc-1) )
  allocate(boundaries(2)%psi(Nm,Nx:Nx+nbc-1,0:Ny+nbc-1), &
           boundaries(2)%chi(Nm,Nx:Nx+nbc-1,0:Ny+nbc-1) )
  allocate(boundaries(3)%psi(Nm,0:Nx+nbc-1,0:nbc-1), &
           boundaries(3)%chi(Nm,0:Nx+nbc-1,0:nbc-1) )
  allocate(boundaries(4)%psi(Nm,0:Nx+nbc-1,Ny:Ny+nbc-1), &
           boundaries(4)%chi(Nm,0:Nx+nbc-1,Ny:Ny+nbc-1) )
  do i=1,4
    boundaries(i)%psi=0.
    boundaries(i)%chi=0.
  enddo  
  allocate(tau_x(Nx,Ny))
  
end subroutine

function get_index_range_inclusive(i1,i2,j1,j2,k1,k2) result(ret)
  integer :: ret, i1,i2, j1,j2,k1,k2

  if(Nx==0.OR.Ny==0.OR.Nm==0) then
    ret=1
    return
  endif
  i1=1
  i2=Nx
  j1=1
  j2=Ny
  k1=1
  k2=Nm
  ret=0
end function

function get_wind_field_index_range_inclusive(i1,i2,j1,j2) result(ret)
  integer :: ret, i1,i2, j1,j2

  if(Nx==0.OR.Ny==0.OR.Nm==0) then
    ret=1
    return
  endif
  i1=1
  i2=Nx
  j1=1
  j2=Ny
  ret=0
end function

function get_boundary_index_range_inclusive(index_of_boundary,i1,i2,j1,j2,k1,k2) result(ret)
  integer :: index_of_boundary,ret, i1,i2, j1,j2,k1,k2

  if(Nx==0.OR.Ny==0.OR.Nm==0) then
    ret=1
    return
  endif
  i1=lbound(boundaries(index_of_boundary)%psi,2)
  i2=ubound(boundaries(index_of_boundary)%psi,2)
  j1=lbound(boundaries(index_of_boundary)%psi,3)
  j2=ubound(boundaries(index_of_boundary)%psi,3)
  k1=lbound(boundaries(index_of_boundary)%psi,1)
  k2=ubound(boundaries(index_of_boundary)%psi,1)
  ret=0
end function

function commit_parameters() result(ret)
  integer :: ret

  Nx = Lx/dx + 1
  Ny = Ly/dy + 1
  !Nt = T/dt + 1
  Nt = INT(T/dt)/savecounter + 1

  call initialize_arrays()

  t_curr   = begin_time

  ret=0
end function

function initialize_grid() result(ret)
  integer :: ret

  if(timestep_method.EQ."euler") then
    ret=initialize_euler()
    return
  endif
  if(timestep_method.EQ."leapfrog") then
    ret=initialize_leapfrog()
    return
  endif
  if(timestep_method.EQ."adams_bashforth") then
    ret=initialize_ab()
    return
  endif
  ret=-1

end function

function initialize_euler() result(ret)
  integer :: ret,i

  chi_prev=0.
! initialize vis_*_prev from psi_2
  call vis_bot(Nm,Nx,Ny,boundary,psi_2,vis_bot_prev,boundaries)
  call vis_lat(Nm,Nx,Ny,boundary,psi_2,vis_lat_prev,boundaries)
! initialize vis_*_curr from psi_1
  call vis_bot(Nm,Nx,Ny,boundary,psi_1,vis_bot_curr,boundaries)
  call vis_lat(Nm,Nx,Ny,boundary,psi_1,vis_lat_curr,boundaries)

! chi is already calculated here so it becomes available for other codes
  call chi(tau,A_H,R_H,Lx,Ly,lambda0,lambda1,e111,phi1z0, &
       &    Nm,Nx,Ny,dx,H,rho,beta0,err_tol,max_it,relax_coef, &
       &    psi_1,boundary,vis_bot_curr,vis_bot_curr,vis_lat_curr,chi_prev, &
       &    chii,boundaries)

  counter=counter+1
  ret=0
end function

function initialize_ab() result(ret)
  integer :: ret,i

  chi_prev=0.
! initialize vis_*_prev from psi_2
  call vis_bot(Nm,Nx,Ny,boundary,psi_2,vis_bot_prev,boundaries)
  call vis_lat(Nm,Nx,Ny,boundary,psi_2,vis_lat_prev,boundaries)
! initialize vis_*_curr from psi_1
  call vis_bot(Nm,Nx,Ny,boundary,psi_1,vis_bot_curr,boundaries)
  call vis_lat(Nm,Nx,Ny,boundary,psi_1,vis_lat_curr,boundaries)

! chi is already calculated here so it becomes available for other codes
  call chi(tau,A_H,R_H,Lx,Ly,lambda0,lambda1,e111,phi1z0, &
       &    Nm,Nx,Ny,dx,H,rho,beta0,err_tol,max_it,relax_coef, &
       &    psi_1,boundary,vis_bot_curr,vis_bot_curr,vis_lat_curr,chi_prev, &
       &    chii,boundaries)

  counter=counter+1
  ret=0
end function

function initialize_leapfrog() result(ret)
  integer :: ret,i

! naive initialization of psi_2
  psi_2=psi_1
  chi_prev=0.
! note in principle the only way to restart consistently is by writing and
! reading psi_2 and chi_prev! 
! note chi_prev does not actually do anything atm
! psi_2 may be initialized by a warm up step...

! actually maybe not a good idea (depends on timestepping)
!  do i=1,Nm
!    if(boundary(1)==2) psi_1(1,1,1:Ny)=boundaries(1)%psi(1,1,1:Ny)
!    if(boundary(2)==2) psi_1(1,Nx,1:Ny)=boundaries(2)%psi(1,Nx,1:Ny)
!    if(boundary(3)==2) psi_1(1,1:Nx,1)=boundaries(3)%psi(1,1:Nx,1)
!    if(boundary(4)==2) psi_1(1,1:Nx,Ny)=boundaries(4)%psi(1,1:Nx,Ny)
!  enddo

! initialize vis_*_prev from psi_2
  call vis_bot(Nm,Nx,Ny,boundary,psi_2,vis_bot_prev,boundaries)
  call vis_lat(Nm,Nx,Ny,boundary,psi_2,vis_lat_prev,boundaries)
! initialize vis_*_curr from psi_1
  call vis_bot(Nm,Nx,Ny,boundary,psi_1,vis_bot_curr,boundaries)
  call vis_lat(Nm,Nx,Ny,boundary,psi_1,vis_lat_curr,boundaries)

! chi is already calculated here so it becomes available for other codes
  call chi(tau,A_H,R_H,Lx,Ly,lambda0,lambda1,e111,phi1z0, &
       &    Nm,Nx,Ny,dx,H,rho,beta0,err_tol,max_it,relax_coef, &
       &    psi_1,boundary,vis_bot_curr,vis_bot_prev,vis_lat_prev,chi_prev, &
       &    chii,boundaries)
!  psi_2=psi_1-dt*chii ! initialization of psi_2 !

  counter=counter+1
  ret=0
end function

function evolve_model(tend) result(ret)
  integer :: ret
  real(8) :: tend
  
  if(timestep_method.EQ."euler") then
    ret=evolve_euler(tend)
    return
  endif
  if(timestep_method.EQ."leapfrog") then
    ret=evolve_leapfrog(tend)
    return
  endif
  if(timestep_method.EQ."adams_bashforth") then
    ret=evolve_ab(tend)
    return
  endif
  
  ret=-1
end function

function evolve_euler(tend) result(ret)
  integer :: ret
  real(8) :: tend
  
  do while ( t_curr < tend)
      
   t_curr = t_curr + 2*dt
       
   psi_2 = psi_1
   psi_1 = psi_1 + 2.d0*dt*chii
      
   chi_prev = chii
   vis_bot_prev  = vis_bot_curr
   call vis_bot(Nm,Nx,Ny,boundary,psi_1,vis_bot_curr,boundaries)
   vis_lat_prev  = vis_lat_curr
   call vis_lat(Nm,Nx,Ny,boundary,psi_1,vis_lat_curr,boundaries)
   call chi(tau,A_H,R_H,Lx,Ly,lambda0,lambda1,e111,phi1z0, &
       &    Nm,Nx,Ny,dx,H,rho,beta0,err_tol,max_it,relax_coef, &
       &    psi_1,boundary,vis_bot_curr,vis_bot_curr,vis_lat_curr,chi_prev, &
       &    chii,boundaries)
   counter = counter + 1 
  enddo
  ret=0
end function

function evolve_ab(tend) result(ret)
  integer :: ret
  real(8) :: tend
  
  do while ( t_curr < tend)
      
   t_curr = t_curr + 2*dt
       
   psi_2 = psi_1
   psi_1 = psi_1 + dt*(3*chii-chi_prev)
      
   chi_prev = chii
   vis_bot_prev  = vis_bot_curr
   call vis_bot(Nm,Nx,Ny,boundary,psi_1,vis_bot_curr,boundaries)
   vis_lat_prev  = vis_lat_curr
   call vis_lat(Nm,Nx,Ny,boundary,psi_1,vis_lat_curr,boundaries)
   call chi(tau,A_H,R_H,Lx,Ly,lambda0,lambda1,e111,phi1z0, &
       &    Nm,Nx,Ny,dx,H,rho,beta0,err_tol,max_it,relax_coef, &
       &    psi_1,boundary,vis_bot_curr,vis_bot_curr,vis_lat_curr,chi_prev, &
       &    chii,boundaries)
   counter = counter + 1 
  enddo
  ret=0
end function


function evolve_leapfrog(tend) result(ret)
  integer :: ret
  integer :: save_num
  real(8) :: tend
  
  do while ( t_curr < tend)
      
   t_curr = t_curr + dt
       
  !Robert-Asselin filter + williams
   inter = psi_2 + 2.d0*dt*chii
   psi_1 = psi_1 + ra_alpha*(inter-2.d0*psi_1+psi_2)*raw_alpha
   psi_2 = inter - ra_alpha*(inter-2.d0*psi_1+psi_2)*(1.d0-raw_alpha)
  !!!
  ! do exactly the same thing again with opposite psi arrays
  ! (specials: leap-frog time stepping, and using previous viscosity)
   t_curr = t_curr + dt
   
   chi_prev = chii
   vis_bot_prev  = vis_bot_curr
   call vis_bot(Nm,Nx,Ny,boundary,psi_2,vis_bot_curr,boundaries)
   vis_lat_prev  = vis_lat_curr
   call vis_lat(Nm,Nx,Ny,boundary,psi_2,vis_lat_curr,boundaries)
   call chi(tau,A_H,R_H,Lx,Ly,lambda0,lambda1,e111,phi1z0, &
       &    Nm,Nx,Ny,dx,H,rho,beta0,err_tol,max_it,relax_coef, &
       &    psi_2,boundary,vis_bot_curr,vis_bot_prev,vis_lat_prev,chi_prev, &
       &    chii,boundaries)
   counter = counter + 1

  !Robert-Asselin filter
   inter = psi_1 + 2.d0*dt*chii
   psi_2 = psi_2 + ra_alpha*(inter-2.d0*psi_2+psi_1)*raw_alpha
   psi_1 = inter - ra_alpha*(inter-2.d0*psi_2+psi_1)*(1.d0-raw_alpha)
  !!!

   chi_prev = chii
  ! update viscosity
   vis_bot_prev  = vis_bot_curr
   call vis_bot(Nm,Nx,Ny,boundary,psi_1,vis_bot_curr,boundaries)
   vis_lat_prev  = vis_lat_curr
   call vis_lat(Nm,Nx,Ny,boundary,psi_1,vis_lat_curr,boundaries)
  !find chi, and take a step
   call chi(tau,A_H,R_H,Lx,Ly,lambda0,lambda1,e111,phi1z0, &
       &    Nm,Nx,Ny,dx,H,rho,beta0,err_tol,max_it,relax_coef, &
       &    psi_1,boundary,vis_bot_curr,vis_bot_prev,vis_lat_prev,chi_prev, &
       &    chii,boundaries)
   counter = counter + 1
 
  enddo
  ret=0
end function

function get_counter(c) result (ret)
  integer :: ret,c
  c=counter
  ret=0
end function

function get_time(tnow) result (ret)
  integer :: ret
  real(8) :: tnow
  tnow=t_curr
  ret=0
end function

function set_time(tnow) result (ret)
  integer :: ret
  real(8) :: tnow
  t_curr=tnow
  ret=0
end function

function get_begin_time(t) result (ret)
  integer :: ret
  real(8) :: t
  t=begin_time
  ret=0
end function

function set_begin_time(t) result (ret)
  integer :: ret
  real(8) :: t
  begin_time=t
  ret=0
end function

function get_wind_sigma(t) result (ret)
  integer :: ret
  real(8) :: t
  t=wind_sigma
  ret=0
end function

function set_wind_sigma(t) result (ret)
  integer :: ret
  real(8) :: t
  wind_sigma=t
  ret=0
end function

function get_ra_alpha(t) result (ret)
  integer :: ret
  real(8) :: t
  t=ra_alpha
  ret=0
end function

function set_ra_alpha(t) result (ret)
  integer :: ret
  real(8) :: t
  ra_alpha=t
  ret=0
end function

function get_raw_alpha(t) result (ret)
  integer :: ret
  real(8) :: t
  t=raw_alpha
  ret=0
end function

function set_raw_alpha(t) result (ret)
  integer :: ret
  real(8) :: t
  raw_alpha=t
  ret=0
end function

function get_timestep_method(t) result (ret)
  integer :: ret
  character(len=15) :: t
  t=timestep_method
  ret=0
end function

function set_timestep_method(t) result (ret)
  integer :: ret
  character(len=15) :: t
  timestep_method=t
  ret=0
end function

function set_boundary_conditions(lowx,highx,lowy,highy) result(ret)
  integer :: ret
  character(len=10) :: lowx,highx,lowy,highy

  ret=0
  if(lowx=="no_slip") then
    boundary(1)=0
  elseif (lowx=="free_slip") then
    boundary(1)=1
  elseif (lowx=="interface") then
    boundary(1)=2
  else
    ret=ret+1
  endif
  
  if(highx=="no_slip") then
    boundary(2)=0
  elseif (highx=="free_slip") then
    boundary(2)=1
  elseif (highx=="interface") then
    boundary(2)=2
  else
    ret=ret+2
  endif
  
  if(lowy=="no_slip") then
    boundary(3)=0
  elseif (lowy=="free_slip") then
    boundary(3)=1
  elseif (lowy=="interface") then
    boundary(3)=2
  else
    ret=ret+4
  endif
  
  if(highy=="no_slip") then
    boundary(4)=0
  elseif (highy=="free_slip") then
    boundary(4)=1
  elseif (highy=="interface") then
    boundary(4)=2
  else
    ret=ret+8
  endif  
  
end function

function get_boundary_conditions(lowx,highx,lowy,highy) result(ret)
  integer :: ret
  character(len=10) :: lowx,highx,lowy,highy

  if(boundary(1)==0) lowx="no_slip"
  if(boundary(1)==1) lowx="free_slip"
  if(boundary(1)==2) lowx="interface"
  if(boundary(2)==0) highx="no_slip"
  if(boundary(2)==1) highx="free_slip"
  if(boundary(2)==2) highx="interface"
  if(boundary(3)==0) lowy="no_slip"
  if(boundary(3)==1) lowy="free_slip"
  if(boundary(3)==2) lowy="interface"
  if(boundary(4)==0) highy="no_slip"
  if(boundary(4)==1) highy="free_slip"
  if(boundary(4)==2) highy="interface"
  
  ret=0
end function


function cleanup_code() result(ret)
  integer :: ret
  if(allocated(psi_1)) then
    deallocate (psi_1)
    deallocate (psi_2)
    deallocate (inter)
    deallocate (chii)
    deallocate (chi_prev)
    deallocate (vis_bot_prev)
    deallocate (vis_bot_curr)
    deallocate (vis_lat_prev)
    deallocate (vis_lat_curr)
    deallocate (max_psi)
  
    deallocate(boundaries(1)%psi, &
             boundaries(1)%chi)
    deallocate(boundaries(2)%psi, &
             boundaries(2)%chi )
    deallocate(boundaries(3)%psi, &
             boundaries(3)%chi )
    deallocate(boundaries(4)%psi, &
             boundaries(4)%chi)
    deallocate(tau_x)
  endif
  ret=0
end function

function get_dpsi_dt(i,j,k,dpsi,n) result(ret)
  integer :: ret,n
  integer :: ii,i(n),j(n),k(n)
  real(8) :: dpsi(n)
  
  do ii=1,n
    dpsi(ii)=chii(k(ii),i(ii),j(ii))
  enddo
  ret=0
end function

function get_psi1_state(i,j,k,psi1,n) result(ret)
  integer :: ret,n
  integer :: ii,i(n),j(n),k(n)
  real(8) :: psi1(n)
  
  do ii=1,n
    psi1(ii)=psi_1(k(ii),i(ii),j(ii))
  enddo
  ret=0
end function

function get_psi2_state(i,j,k,psi2,n) result(ret)
  integer :: ret,n
  integer :: ii,i(n),j(n),k(n)
  real(8) :: psi2(n)
  
 do ii=1,n
    psi2(ii)=psi_2(k(ii),i(ii),j(ii))
  enddo
  ret=0
end function

function get_wind_field(i,j,t,n) result(ret)
  integer :: ret,n
  integer :: ii,i(n),j(n)
  real(8) :: t(n)
  
  do ii=1,n
    t(ii)=tau_x(i(ii),j(ii))
  enddo
  ret=0
end function

function set_wind_field(i,j,t,n) result(ret)
  integer :: ret,n
  integer :: ii,i(n),j(n)
  real(8) :: t(n)
  
  do ii=1,n
    tau_x(i(ii),j(ii))=t(ii)
  enddo
  ret=0
end function

function set_psi1_state(i,j,k,psi1,n) result(ret)
  integer :: ret,n
  integer :: ii,i(n),j(n),k(n)
  real(8) :: psi1(n)
  
 do ii=1,n
   psi_1(k(ii),i(ii),j(ii))= psi1(ii)
  enddo
  ret=0
  ret=0
end function

function set_psi2_state(i,j,k,psi2,n) result(ret)
   integer :: ret,n
  integer :: ii,i(n),j(n),k(n)
  real(8) :: psi2(n)

  do ii=1,n
    psi_2(k(ii),i(ii),j(ii))=psi2(ii)
  enddo
  ret=0
end function

function get_boundary_state(i,j,k,index_of_boundary,psi,dpsi,n) result(ret)
  integer :: ret,n
  integer :: ii,i(n),j(n),k(n),index_of_boundary(n)
  real(8) :: psi(n),dpsi(n)

  ret=n
  do ii=1,n
    if(.NOT.(index_of_boundary(ii).GE.1.AND.index_of_boundary(ii).LE.4)) cycle
    if(k(ii).LT.lbound(boundaries(index_of_boundary(ii))%psi,1).OR. &
       k(ii).GT.ubound(boundaries(index_of_boundary(ii))%psi,1)) cycle
    if(i(ii).LT.lbound(boundaries(index_of_boundary(ii))%psi,2).OR. &
       i(ii).GT.ubound(boundaries(index_of_boundary(ii))%psi,2)) cycle
    if(j(ii).LT.lbound(boundaries(index_of_boundary(ii))%psi,3).OR. &
       j(ii).GT.ubound(boundaries(index_of_boundary(ii))%psi,3)) cycle
    psi(ii)=boundaries(index_of_boundary(ii))%psi(k(ii),i(ii),j(ii))
    dpsi(ii)=boundaries(index_of_boundary(ii))%chi(k(ii),i(ii),j(ii))
    ret=ret-1
  enddo
end function

function set_boundary_state(i,j,k,psi,dpsi,index_of_boundary,n) result(ret)
  integer :: ret,n
  integer :: ii,i(n),j(n),k(n),index_of_boundary(n)
  real(8) :: psi(n),dpsi(n)

  ret=n
  do ii=1,n
    if(.NOT.(index_of_boundary(ii).GE.1.AND.index_of_boundary(ii).LE.4)) cycle
    if(k(ii).LT.lbound(boundaries(index_of_boundary(ii))%psi,1).OR. &
       k(ii).GT.ubound(boundaries(index_of_boundary(ii))%psi,1)) cycle
    if(i(ii).LT.lbound(boundaries(index_of_boundary(ii))%psi,2).OR. &
       i(ii).GT.ubound(boundaries(index_of_boundary(ii))%psi,2)) cycle
    if(j(ii).LT.lbound(boundaries(index_of_boundary(ii))%psi,3).OR. &
       j(ii).GT.ubound(boundaries(index_of_boundary(ii))%psi,3)) cycle
    boundaries(index_of_boundary(ii))%psi(k(ii),i(ii),j(ii))=psi(ii)
    boundaries(index_of_boundary(ii))%chi(k(ii),i(ii),j(ii))=dpsi(ii)
    ret=ret-1
  enddo
end function

function get_position_of_index(i,j,k,x,y,n) result(ret)
  integer :: ret,n
  integer :: ii,i(n),j(n),k(n)
  real(8) :: x(n),y(n)
  
  do ii=1,n
    x(ii)=(i(ii)-1)*dx
    y(ii)=(j(ii)-1)*dy
  enddo
  ret=0
end function

function get_wind_field_position_of_index(i,j,x,y,n) result(ret)
  integer :: ret,n
  integer :: ii,i(n),j(n)
  real(8) :: x(n),y(n)
  
  do ii=1,n
    x(ii)=(i(ii)-1)*dx
    y(ii)=(j(ii)-1)*dy
  enddo
  ret=0
end function


function get_index_of_position(x,y,i,j,n) result(ret)
  integer :: ret,n
  integer :: ii,i(n),j(n)
  real(8) :: x(n),y(n)
  
  do ii=1,n
    i(ii)=floor(x(ii)/dx+0.5)+1
    j(ii)=floor(y(ii)/dy+0.5)+1
  enddo
  ret=0
end function

function get_boundary_position_of_index(i,j,k,index_of_boundary,x,y,n) result(ret)
  integer :: ret,n
  integer :: ii,i(n),j(n),k(n),index_of_boundary(n)
  real(8) :: x(n),y(n)
! its the same as for the normal grid!  
  do ii=1,n
    x(ii)=(i(ii)-1)*dx
    y(ii)=(j(ii)-1)*dy
  enddo
  ret=0
end function

function area_averaging_grid_sample(d,x,y,k,psi,dpsi) result(ret)
  integer i,j,ret,k
  integer imin,imax,jmin,jmax,ii,jj
  real(8) :: d,x,y,psi,dpsi,cell_volume,overlap
  real(8) :: xc1,xc2,yc1,yc2
  real(8) :: xg1,xg2,yg1,yg2
  cell_volume=d*d
  
  xc1=x-d/2
  xc2=x+d/2
  yc1=y-d/2
  yc2=y+d/2
  imin=floor(xc1/dx+0.5)+1
  imax=floor(xc2/dx+0.5)+1
  jmin=floor(yc1/dy+0.5)+1
  jmax=floor(yc2/dy+0.5)+1
  psi=0
  dpsi=0  
  do i=imin,imax
    do j=jmin,jmax
      xg1=(i-1)*dx-dx/2
      xg2=i*dx-dx/2
      yg1=(j-1)*dy-dy/2
      yg2=j*dy-dy/2
      overlap=max(0.,min(xc2,xg2)-max(xc1,xg1))*max(0.,min(yc2,yg2)-max(yc1,yg1))
      ii=min(max(i,1),Nx)
      jj=min(max(j,1),Ny)
!      print*,"overlap",overlap,ii,jj,psi_1(k,ii,jj)
      psi=psi+overlap*psi_1(k,ii,jj)
      dpsi=dpsi+overlap*chii(k,ii,jj)
    enddo
  enddo
  psi=psi/cell_volume
  dpsi=dpsi/cell_volume
  
  ret=0
end function

function get_psi_state_at_point(d,x,y,k,psi,dpsi,n) result(ret)
  integer :: ret,n,i
  integer :: k(n)
  real(8) :: d(n),x(n),y(n),psi(n),dpsi(n),p,dp

  if(Nx==0.OR.Ny==0.OR.Nm==0) then
    ret=-1
    return
  endif
  
  ret=0
  do i=1,n
    ret=ret+area_averaging_grid_sample(d(i),x(i),y(i),k(i),p,dp)
    psi(i)=p
    dpsi(i)=dp
  enddo
  
end function

function recommit_parameters() result(ret)
  integer :: ret
  ret=0
end function

! here come the parameter getters and setters

function get_Lx(x) result (ret)
  integer :: ret
  real(8) :: x
  x=Lx
  ret=0
end function

function set_Lx(x) result (ret)
  integer :: ret
  real(8) :: x
  Lx=x
  ret=0
end function

function get_Ly(x) result (ret)
  integer :: ret
  real(8) :: x
  x=Ly
  ret=0
end function

function set_Ly(x) result (ret)
  integer :: ret
  real(8) :: x
  Ly=x
  ret=0
end function

function get_dy(x) result (ret)
  integer :: ret
  real(8) :: x
  x=dy
  ret=0
end function

function set_dy(x) result (ret)
  integer :: ret
  real(8) :: x
  dy=x
  ret=0
end function

function get_dx(x) result (ret)
  integer :: ret
  real(8) :: x
  x=dx
  ret=0
end function

function set_dx(x) result (ret)
  integer :: ret
  real(8) :: x
  dx=x
  ret=0
end function

function get_dt(x) result (ret)
  integer :: ret
  real(8) :: x
  x=2*dt
  ret=0
end function

function set_dt(x) result (ret)
  integer :: ret
  real(8) :: x
  dt=x/2
  ret=0
end function

function get_T(x) result (ret)
  integer :: ret
  real(8) :: x
  x=T
  ret=0
end function

function set_T(x) result (ret)
  integer :: ret
  real(8) :: x
  T=x
  ret=0
end function

function get_A_H(x) result (ret)
  integer :: ret
  real(8) :: x
  x=A_H
  ret=0
end function

function set_A_H(x) result (ret)
  integer :: ret
  real(8) :: x
  A_H=x
  ret=0
end function

function get_R_H(x) result (ret)
  integer :: ret
  real(8) :: x
  x=R_H
  ret=0
end function

function set_R_H(x) result (ret)
  integer :: ret
  real(8) :: x
  R_H=x
  ret=0
end function

function get_lambda0(x) result (ret)
  integer :: ret
  real(8) :: x
  x=lambda0
  ret=0
end function

function set_lambda0(x) result (ret)
  integer :: ret
  real(8) :: x
  lambda0=x
  ret=0
end function

function get_lambda1(x) result (ret)
  integer :: ret
  real(8) :: x
  x=lambda1
  ret=0
end function

function set_lambda1(x) result (ret)
  integer :: ret
  real(8) :: x
  lambda1=x
  ret=0
end function

function get_e111(x) result (ret)
  integer :: ret
  real(8) :: x
  x=e111
  ret=0
end function

function set_e111(x) result (ret)
  integer :: ret
  real(8) :: x
  e111=x
  ret=0
end function

function get_phi1z0(x) result (ret)
  integer :: ret
  real(8) :: x
  x=phi1z0
  ret=0
end function

function set_phi1z0(x) result (ret)
  integer :: ret
  real(8) :: x
  phi1z0=x
  ret=0
end function

function get_H(x) result (ret)
  integer :: ret
  real(8) :: x
  x=H
  ret=0
end function

function set_H(x) result (ret)
  integer :: ret
  real(8) :: x
  H=x
  ret=0
end function

function get_rho(x) result (ret)
  integer :: ret
  real(8) :: x
  x=rho
  ret=0
end function

function set_rho(x) result (ret)
  integer :: ret
  real(8) :: x
  rho=x
  ret=0
end function

function get_beta0(x) result (ret)
  integer :: ret
  real(8) :: x
  x=beta0
  ret=0
end function

function set_beta0(x) result (ret)
  integer :: ret
  real(8) :: x
  beta0=x
  ret=0
end function

function get_savecounter(x) result (ret)
  integer :: ret,x
  x=savecounter
  ret=0
end function

function set_savecounter(x) result (ret)
  integer :: ret,x
  savecounter=x
  ret=0
end function

function get_tau(x) result (ret)
  integer :: ret
  real(8) :: x
  x=tau
  ret=0
end function

function set_tau(x) result (ret)
  integer :: ret
  real(8) :: x
  tau=x
  ret=0
end function

function get_err_tol(x) result (ret)
  integer :: ret
  real(8) :: x
  x=err_tol
  ret=0
end function

function set_err_tol(x) result (ret)
  integer :: ret
  real(8) :: x
  err_tol=x
  ret=0
end function

function get_max_it(x) result (ret)
  integer :: ret,x
  x=max_it
  ret=0
end function

function set_max_it(x) result (ret)
  integer :: ret,x
  max_it=x
  ret=0
end function

function get_relax_coef(x) result (ret)
  integer :: ret
  real(8) :: x
  x=relax_coef
  ret=0
end function

function set_relax_coef(x) result (ret)
  integer :: ret
  real(8) :: x
  relax_coef=x
  ret=0
end function

function get_Nx(x) result (ret)
  integer :: ret
  integer :: x
  x=Nx
  ret=0
end function
function get_Ny(x) result (ret)
  integer :: ret
  integer :: x
  x=Ny
  ret=0
end function
function get_Nm(x) result (ret)
  integer :: ret
  integer :: x
  x=Nm
  ret=0
end function
function set_Nm(x) result (ret)
  integer :: ret,x
  Nm=x
  ret=0
end function
function get_Nt(x) result (ret)
  integer :: ret
  integer :: x
  x=Nt
  ret=0
end function

function set_interface_wind(x) result (ret)
  integer :: ret,x
  interface_wind=x
  ret=0
end function
function get_interface_wind(x) result (ret)
  integer :: ret
  integer :: x
  x=interface_wind
  ret=0
end function


end module




! move wind to here in order to...
subroutine wind(Nm,Nx,Ny,windy)
 use qgmodel, only: wind_sigma, interface_wind,tau_x,tau
implicit none
integer, intent(in) :: Nx,Ny,Nm
real(8), dimension(Nm,Nx,Ny), intent(out) :: windy

integer :: i,j,m
real(8) :: pi = 3.14159265358979d0
!real(8), dimension(Nx,Ny) :: tau

windy(:,:,:) = 0.d0

if(.not.allocated(tau_x)) stop "tau_x not allocated"


if(interface_wind.EQ.0) then
  tau_x(:,:)     = 0.d0
  do i=1,Nx
    do j=1,Ny
      if(wind_sigma<-1.) then
  ! jan's:
    tau_x(i,j) = cos(2.*pi*((j-1.)/(Ny-1.)-0.5))+2.*sin(pi*((j-1.)/(Ny-1.)-0.5))
      else
  ! dijkstra:
    tau_x(i,j)= - ( wind_sigma*cos(pi*(j-1.)/(Ny-1.))+(1-wind_sigma)*cos(2.*pi*((j-1.)/(Ny-1.))) )
      endif
    end do
  end do
endif

do m=1,Nm
  do i=2,Nx-1
    do j=2,Ny-1
    windy(m,i,j) = -tau_x(i,j+1)+tau_x(i,j-1)    ! include the whole curl, i.e., also x-derivative for generality ?????
    end do
  end do
end do

! fix wind field normalization
if(interface_wind.EQ.0) then
  tau_x=tau*tau_x
else
  if(tau.NE.0) windy=windy/tau
endif

! boundary wind values not needed atm

end


!      real(8), allocatable, dimension (:,:) :: wind_term
!      allocate (wind_term(Nx,Ny))
!      call wind(Nx,Ny,wind_term)



