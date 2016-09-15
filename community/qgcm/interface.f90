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
    
  if(oceanonly.NE.1) then 
    pam=pa-dpa_dt*dta
    astm=ast-dast_dt*dta
    hmixam=hmixa-dhmixa_dt*dta
  endif
  if(atmosonly.NE.1) then
    pom=po-dpo_dt*dto
    sstm=sst-dsst_dt*dto
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
