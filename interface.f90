module adcirc_interface
  use amuse_adcirc

contains
  
function get_node_state(ind,eta2_,uu2_,vv2_) result(ret)
  integer :: ind,ret
  real*8 :: eta2_,uu2_,vv2_
  
  eta2_=ETA2(ind)
  uu2_=UU2(ind)
  vv2_=VV2(ind)
  ret=0
end function

function get_node_position(ind,x_,y_) result(ret)
  integer :: ind,ret
  real*8 :: x_,y_
  
  x_=X(ind)
  y_=Y(ind)
  ret=0
end function

function get_element_nodes(ind,n1,n2,n3) result(ret)
  integer :: ind,n1,n2,n3,ret
  if(ind.LT.1.OR.ind.GT.NE) then
    ret=-1
    return
  endif
  n1=NM(ind,1)
  n2=NM(ind,2)
  n3=NM(ind,3)
  ret=0
end function
  
end module
