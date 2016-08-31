module qgcm_interface
  use amuse_qgcm

  double precision :: begin_time=0.

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
  
  call amuse_initialize

  call parameter_defaults
  
  ret=0
end function

function commit_grids() result(ret)
  integer :: ret
  
  call init_grids()
  
  ret=0
end function

function evolve_model(tend) result(ret)
  integer :: ret
  real*8 :: tend
  
  call mainloop(nsteps0, nsteps)
  
  ret=0  
end function

function commit_parameters() result(ret)
  integer :: ret
  
  call check_parameters()
  
  tnow=nsteps0
  
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
