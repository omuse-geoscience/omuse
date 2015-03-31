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
  if(ind.LT.1.OR.ind.GT.NP) then
    ret=-1
    return
  endif
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
  
function get_element_position(ind,x_,y_) result(ret)
  integer :: ind,n1,n2,n3,ret
  real*8 :: x_,y_

  if(ind.LT.1.OR.ind.GT.NE) then
    ret=-1
    return
  endif
  n1=NM(ind,1)
  n2=NM(ind,2)
  n3=NM(ind,3)  
  x_=(X(n1)+X(n2)+X(n3))/3
  y_=(Y(n1)+Y(n2)+Y(n3))/3
  ret=0
end function
  
function get_number_of_nodes(nout) result(ret)
  integer :: ret,nout
  nout=NP
  ret=0  
end function
  
function get_number_of_elements(nout) result(ret)
  integer :: ret,nout
  nout=NE
  ret=0  
end function

function get_number_of_elevation_boundary_segments(nout) result(ret)
  integer :: ret, nout
  nout=NOPE
  ret=0
end function

function get_number_of_nodes_in_elevation_boundary_segment(iseg,nout) result(ret)
  integer :: ret, iseg, nout
  nout=NVDLL(iseg)
  ret=0
end function
  
function get_elevation_boundary_node(inode,iseg,index_node) result(ret)
  integer :: ret, inode,iseg,index_node  
  index_node=NBDV(iseg,inode)
  ret=0
end function

function get_number_of_flow_boundary_segments(nout) result(ret)
  integer :: ret, nout
  nout=NBOU
  ret=0
end function

function get_number_of_nodes_in_flow_boundary_segment(iseg,nout) result(ret)
  integer :: ret, iseg, nout
  nout=NVELL(iseg)
  ret=0
end function
  
function get_flow_boundary_node(inode,iseg,index_node) result(ret)
  integer :: ret, inode,iseg,index_node  
  index_node=NBVV(iseg,inode)
  ret=0
end function

function get_flow_boundary_type(inode,iseg,type_) result(ret)
  integer :: ret, inode,iseg,type_  
  type_=IBTYPE(iseg)
  ret=0
end function





end module


