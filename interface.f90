module adcirc_interface
  use amuse_adcirc

contains

function test_amuse() result(ret)
  integer :: ret
  ret=test()
end function
  
end module
