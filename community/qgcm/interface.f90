module qgcm_interface
  use amuse_qgcm

contains

include "getter_setters.f90"

subroutine parameter_defaults

include "default_param.f90"  

end subroutine

function initialize_code() result(ret)
  integer :: ret
  
  call amuse_initialize

  call parameter_defaults
  
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
