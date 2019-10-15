
function initialize() result(ret)
  use dflowfm_omuse_lib
  integer :: result
  character(len=256) :: somefile = "westerscheldt.mdu"

  ret=initialize_dflowfm(somefile)

end function

