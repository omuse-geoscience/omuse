module swan_interface
  use amuse_swan
  implicit none

contains

function initialize_code() result(ret)
  integer :: ret
  ret=swan_entry()
  if(ret.EQ.0) ret=swan_init()
  PROJID='AMUSE'
  PROJNR='1'
  
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
  ret=swan_cleanup()
end function

end module
