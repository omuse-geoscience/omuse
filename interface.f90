module adcirc_interface
  use amuse_adcirc

contains

function test_amuse() result(ret)
  integer :: ret
  ret=test()
end function
  
function get_node_state(ind,eta2_,uu2_,vv2_) result(ret)
  integer :: ind,ret
  real*8 :: eta2_,uu2_,vv2_
  
  eta2_=ETA2(ind)
  uu2_=UU2(ind)
  vv2_=VV2(ind)
  ret=0
end function


end module
