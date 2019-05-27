module qgcm_interface
  use omuse_qgcm
  implicit none

  real*8 :: begin_time=0.

  logical :: use_interface_forcings=.FALSE.

contains

include "getter_setters.f90"

subroutine parameter_defaults

include "default_param.f90"  

end subroutine

function get_model_time(t) result(ret)
  integer :: ret
  real*8 :: t
  t=tday
  ret=0
end function

function get_ocean_time(t) result(ret)
  integer :: ret
  real*8 :: t
  t=tocean
  ret=0
end function
function get_ocean_prev_time(t) result(ret)
  integer :: ret
  real*8 :: t
  t=tocean-(tos-tosm)
  ret=0
end function

function initialize_code() result(ret)
  integer :: ret
  
  call omuse_initialize

  call parameter_defaults
  
  ret=0
end function

function initialize_grids() result(ret)
  integer :: ret
  
!~   call init_topo_rad_state
  
  ret=0
end function

function commit_grids() result(ret)
  integer :: ret
   
  call init_topo_rad_state

  ! needs to be improved for exact restarts, i.e.
  ! - "replay" timestepping to tini to get current dto, dta
  ! - initialize pom, pam etc from this
  tos=begin_time
  tosm=begin_time-dto

  tas=begin_time
  tasm=begin_time-dta

  if(oceanonly.NE.1) then 
    pam=pa-dpa_dt*(tas-tasm)
    astm=ast-dast_dt*(tas-tasm)
    hmixam=hmixa-dhmixa_dt*(tas-tasm)
  endif
  if(atmosonly.NE.1) then
    pom=po-dpo_dt*(tos-tosm)
    sstm=sst-dsst_dt*(tos-tosm)
  endif
  
  tini=begin_time/secsyr

  call init_grids1
  if(.not.use_interface_forcings) call read_forcings
  call init_grids2
  
  ret=0
end function

function evolve_model(tend) result(ret)
  integer :: ret
  real*8 :: tend,tnow

  tnow=tday
  if(tnow.LT.tend) then
    ! calc number of timesteps, ensuring always a ocean step
    nsteps=nsteps0+nstr*ceiling( secday*(tend-tnow)/(nstr*dta) )
    print*, nt, nsteps0,nsteps
    call mainloop(nsteps0, nsteps)
    nsteps0=nsteps    
    print*, nt, nsteps0,nsteps,tday,ntdone
    print*,"timestep done"
    ret=print_diag()
  endif

  if(oceanonly.NE.1) then 
    dpa_dt=(pa-pam)/dta
    dast_dt=(ast-astm)/dta
    dhmixa_dt=(hmixa-hmixam)/dta
  endif
  if(atmosonly.NE.1) then
    dpo_dt=(po-pom)/dto
    dsst_dt=(sst-sstm)/dto
  endif
    
  ret=0  
end function

function commit_parameters() result(ret)
  integer :: ret
  
  call check_parameters()
    
  ret=0
end function

function recommit_parameters() result(ret)
  integer :: ret
    
  ret=-2
end function

function cleanup_code() result(ret)
  integer :: ret
  
  call finalize
  
  ret=0
end function


end module
