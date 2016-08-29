module qgcm_interface
  use amuse_qgcm

contains

include "getter_setters.f90"

function initialize_code() result(ret)
  integer :: ret
  ret=0
end function
function commit_parameters() result(ret)
  integer :: ret
  ret=-2
end function
function recommit_parameters() result(ret)
  integer :: ret
  ret=-2
end function
function cleanup_code() result(ret)
  integer :: ret
  ret=-2
end function


end module
