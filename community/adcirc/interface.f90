module adcirc_interface
  use amuse_adcirc
  implicit none

contains

function set_rootdir(rootdir_) result(ret)
  integer :: ret
  character(256) :: rootdir_
  ROOTDIR=rootdir_
  ret=0
end function
function get_rootdir(rootdir_) result(ret)
  integer :: ret
  character(256) :: rootdir_
  rootdir_=ROOTDIR
  ret=0
end function

function initialize_code() result(ret)
  integer :: ret
  ROOTDIR='.'
  ret=0
end function

function evolve_model(tend) result(ret)
  integer :: ret
  real(8) :: tend
  do while(ITIME_BGN*DTDP+STATIM*86400.D0<tend-DTDP/2)
    call ADCIRC_Run(1)
  enddo
  if(.not.use_interface_elevation_boundary) ESBIN(1:NETA)=ETA2(NBD(1:NETA))
  if(.not.use_interface_met_forcing) then
    WSX(1:NP)=WSX2(1:NP)
    WSY(1:NP)=WSY2(1:NP)  
  endif
  ret=0
end function

function cleanup_code() result(ret)
  integer :: ret
  call ADCIRC_Final
  ret=0  
end function

function commit_parameters() result(ret)
  integer :: ret
  CALL ADCIRC_Init(ROOTD=ROOTDIR)
  allocate( ESBIN(MNETA) )
  ESBIN(:)=0.
  allocate( WSX(MNP),WSY(MNP) )
  WSX(:)=0.
  WSY(:)=0.
  ret=0
end function

function recommit_parameters() result(ret)
 integer :: ret
 ret=-2
end function

function get_model_time(tout) result(ret)
  integer :: ret
  real(8) :: tout  
  tout = STATIM*86400.D0+DTDP*ITIME_BGN
  ret=0
end function

function get_timestep(dtout) result(ret)
  integer :: ret
  real(8) :: dtout  
  dtout = DTDP  
  ret=0
end function
  
function get_use_interface_elevation_boundary(flag) result(ret)
  integer :: ret
  logical flag
  flag=use_interface_elevation_boundary
  ret=0
end function  
function set_use_interface_elevation_boundary(flag) result(ret)
  integer :: ret
  logical flag
  use_interface_elevation_boundary=flag
  ret=0
end function  

function get_use_interface_met_forcing(flag) result(ret)
  integer :: ret
  logical flag
  flag=use_interface_met_forcing
  ret=0
end function  
function set_use_interface_met_forcing(flag) result(ret)
  integer :: ret
  logical flag
  use_interface_met_forcing=flag
  ret=0
end function  

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

function get_node_coriolis_f(ind,corif_) result(ret)
  integer :: ind,ret
  real*8 :: corif_
  corif_=CORIF(ind)
  ret=0
end function
function set_node_coriolis_f(ind,corif_) result(ret)
  integer :: ind,ret
  real*8 :: corif_
  CORIF(ind)=corif_
  ret=0
end function

function get_node_wind_stress(ind,wsx_,wsy_) result(ret)
  integer :: ind,ret
  real*8 :: wsx_,wsy_
  wsx_=WSX(ind)*RhoWat0
  wsy_=WSY(ind)*RhoWat0
  ret=0
end function
function set_node_wind_stress(ind,wsx_,wsy_) result(ret)
  integer :: ind,ret
  real*8 :: wsx_,wsy_
  WSX(ind)=wsx_/RhoWat0
  WSY(ind)=wsy_/RhoWat0
  ret=0
end function

function get_node_sigma(ind,indz,s_) result(ret)
  integer :: ind,indz,ret
  real*8 :: s_
  if(ind.LT.1.OR.ind.GT.NP) then
    ret=-1
    return
  endif
  if(indz.LT.1.OR.indz.GT.NFEN) then
    ret=-1
    return
  endif
  s_=SIGMA(indz)
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

function get_elevation_boundary_eta(inode,iseg,eta) result(ret)
  integer :: ret, inode,iseg,index_node,i  
  real*8 :: eta
  i=sum(NVDLL(1:iseg-1))+inode
  if(NBD(i).NE.NBDV(iseg,inode)) then
    ret=-1
    return
  endif
  eta=ESBIN(i)
  ret=0
end function

function set_elevation_boundary_eta(inode,eta,iseg) result(ret)
  integer :: ret, inode,iseg,index_node,i  
  real*8 :: eta
  i=sum(NVDLL(1:iseg-1))+inode
  if(NBD(i).NE.NBDV(iseg,inode)) then
    ret=-1
    return
  endif
  ESBIN(i)=eta
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

function get_number_of_vertical_nodes(nout) result(ret)
  integer :: ret,nout
  nout=NFEN
  ret=0  
end function


end module


