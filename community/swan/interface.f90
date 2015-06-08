module swan_interface
  use amuse_swan
  implicit none

contains

function initialize_code() result(ret)
  integer :: ret
  PROJID='AMUSE'
  PROJNR='1'
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
