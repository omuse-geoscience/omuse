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
  integer :: i,ret
  real(8) :: tend
  do while((ITIME_BGN-1)*DTDP+STATIM*86400.D0<tend-DTDP/2)
    call ADCIRC_Run(1)
  enddo
  DETA_DT(1:NP)=(ETA2(1:NP)-ETA1(1:NP))/DTDP
  if(.not.use_interface_elevation_boundary) ESBIN(1:NETA)=ETA2(NBD(1:NETA))
  if(.not.use_interface_flow_boundary) then
    FBQX(1:NVEL)=VV2(NBV(1:NVEL))*(DP(NBV(1:NVEL))+IFNLFA*ETA2(NBV(1:NVEL)))
    FBQY(1:NVEL)=UU2(NBV(1:NVEL))*(DP(NBV(1:NVEL))+IFNLFA*ETA2(NBV(1:NVEL)))
  else
    do i=1, NVEL
      if(LBCODEI(i).NE.2.AND.LBCODEI(i).NE.12.AND.LBCODEI(i).NE.22) then
        FBQX(i)=VV2(NBV(i))*(DP(NBV(i))+IFNLFA*ETA2(NBV(i)))
        FBQY(i)=UU2(NBV(i))*(DP(NBV(i))+IFNLFA*ETA2(NBV(i)))
      endif
    enddo
  endif
  if(.not.use_interface_met_forcing) then
    WSX(1:NP)=WSX2(1:NP)
    WSY(1:NP)=WSY2(1:NP)  
  endif
  if(.not.use_interface_wave_forcing) then
    RSNX(1:NP)=RSNX2(1:NP)
    RSNY(1:NP)=RSNY2(1:NP)  
  endif
  if(.not.use_interface_tidal_forcing.AND.CTIP) TIP(1:NP)=TIP2(1:NP) 
  ret=0
end function

function cleanup_code() result(ret)
  integer :: ret
  call ADCIRC_Final
  ret=0  
end function

subroutine ADCIRC_allocate_interface_storage()
  if(allocated(ESBIN)) deallocate(ESBIN)
  if(allocated(FBQX)) deallocate(FBQX)
  if(allocated(FBQY)) deallocate(FBQY)
  if(allocated(WSX)) deallocate(WSX)
  if(allocated(WSY)) deallocate(WSY)
  if(allocated(RSNX)) deallocate(RSNX)
  if(allocated(RSNY)) deallocate(RSNY)
  if(allocated(TIP)) deallocate(TIP)
  if(allocated(DETA_DT)) deallocate(DETA_DT)
  if(allocated(PR)) deallocate(PR)
  allocate( ESBIN(MNETA) )
  allocate( FBQX(MNVEL),FBQY(MNVEL) )
  allocate( WSX(MNP),WSY(MNP), DETA_DT(MNP), PR(MNP) )
  allocate( RSNX(MNP),RSNY(MNP) )
  allocate( TIP(MNP) )
  ESBIN(:)=0.
  FBQX(:)=0.
  FBQY(:)=0.
  WSX(:)=0.
  WSY(:)=0.
  RSNX(:)=0.
  RSNY(:)=0.
  DETA_DT(:)=0.
  PR(:)=0.
  TIP(:)=0.

  if(allocated(AMUSE_LNM)) deallocate(AMUSE_LNM)
  if(allocated(AMUSE_RESSAL)) deallocate(AMUSE_RESSAL)
  if(allocated(AMUSE_RESTEMP)) deallocate(AMUSE_RESTEMP)
  if(allocated(AMUSE_HFLUX)) deallocate(AMUSE_HFLUX)
  allocate( AMUSE_LNM(MNETA) )
  allocate( AMUSE_RESSAL(MNETA, NFEN) )
  allocate( AMUSE_RESTEMP(MNETA, NFEN) )
  allocate( AMUSE_HFLUX(MNP) )
  AMUSE_LNM=0.
  AMUSE_RESSAL=0.
  AMUSE_RESTEMP=0.
  AMUSE_HFLUX=0.

end subroutine

function commit_parameters() result(ret)
  integer :: ret
  CALL ADCIRC_Init(ROOTD=ROOTDIR)
  call ADCIRC_allocate_interface_storage()
  ret=0
end function

function commit_grid() result(ret)
  integer :: ret
  integer :: i
  REAL(SZ) H2

  ETA1(1:NP)=ETA2(1:NP)-DTDP*DETA_DT(1:NP)

  DO I=1, NP
     ETAS(I)=ETA2(I)-ETA1(I)
     H2=DP(I)+IFNLFA*ETA2(I)
     QX2(I)=UU2(I)*H2
     QY2(I)=VV2(I)*H2
  END DO

  ret=0
end function


function recommit_parameters() result(ret)
 integer :: ret
 ret=-2
end function

function get_model_time(tout) result(ret)
  integer :: ret
  real(8) :: tout  
  tout = STATIM*86400.D0+DTDP*(ITIME_BGN-1)
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

function get_use_interface_flow_boundary(flag) result(ret)
  integer :: ret
  logical flag
  flag=use_interface_flow_boundary
  ret=0
end function  
function set_use_interface_flow_boundary(flag) result(ret)
  integer :: ret
  logical flag
  use_interface_flow_boundary=flag
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

function get_use_interface_wave_forcing(flag) result(ret)
  integer :: ret
  logical flag
  flag=use_interface_wave_forcing
  ret=0
end function  
function set_use_interface_wave_forcing(flag) result(ret)
  integer :: ret
  logical flag
  use_interface_wave_forcing=flag
  ret=0
end function  

function get_use_interface_tidal_forcing(flag) result(ret)
  integer :: ret
  logical flag
  flag=use_interface_tidal_forcing
  ret=0
end function  
function set_use_interface_tidal_forcing(flag) result(ret)
  integer :: ret
  logical flag
  use_interface_tidal_forcing=flag
  ret=0
end function

function get_use_interface_surface_heat_forcing(flag) result(ret)
  integer :: ret
  logical flag
  flag=use_interface_surface_heat_forcing
  ret=0
end function  
function set_use_interface_surface_heat_forcing(flag) result(ret)
  integer :: ret
  logical flag
  use_interface_surface_heat_forcing=flag
  ret=0
end function

function get_use_interface_lnm_boundary(flag) result(ret)
  integer :: ret
  logical flag
  flag=use_interface_lnm_boundary
  ret=0
end function  
function set_use_interface_lnm_boundary(flag) result(ret)
  integer :: ret
  logical flag
  use_interface_lnm_boundary=flag
  ret=0
end function

function get_use_interface_salinity_boundary(flag) result(ret)
  integer :: ret
  logical flag
  flag=use_interface_salinity_boundary
  ret=0
end function  
function set_use_interface_salinity_boundary(flag) result(ret)
  integer :: ret
  logical flag
  use_interface_salinity_boundary=flag
  ret=0
end function

function get_use_interface_temperature_boundary(flag) result(ret)
  integer :: ret
  logical flag
  flag=use_interface_temperature_boundary
  ret=0
end function  
function set_use_interface_temperature_boundary(flag) result(ret)
  integer :: ret
  logical flag
  use_interface_temperature_boundary=flag
  ret=0
end function

function get_node_eta(ind,eta2_) result(ret)
  integer :: ind,ret
  real*8 :: eta2_
  eta2_=ETA2(ind)
  if(NNODECODE(ind).EQ.0) then
    eta2_=-DP(ind)
  endif
  ret=0
end function

function get_node_vx(ind,uu2_) result(ret)
  integer :: ind,ret
  real*8 :: uu2_
  uu2_=UU2(ind)
  ret=0
end function

function get_node_vy(ind,vv2_) result(ret)
  integer :: ind,ret
  real*8 :: vv2_
  vv2_=VV2(ind)
  ret=0
end function

function get_node_eta_prev(ind,eta_) result(ret)
  integer :: ind,ret
  real*8 :: eta_  
  eta_=ETA1(ind)
  ret=0
end function

function get_node_deta_dt(ind,x) result(ret)
  integer :: ind,ret
  real*8 :: x
  x=DETA_DT(ind)
  ret=0
end function

function get_node_status(ind,x) result(ret)
  integer :: ind,ret
  character(len=*) :: x
  if(NNODECODE(ind).EQ.0) then
    x="dry"
  else
    x="wet"
  endif
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
function get_node_coordinates(ind,lon_,lat_) result(ret)
  integer :: ind,ret
  real*8 :: lon_,lat_
  if(ind.LT.1.OR.ind.GT.NP) then
    ret=-1
    return
  endif
  lon_=SLAM(ind)
  lat_=SFEA(ind)
  ret=0
end function

function set_node_eta(ind,eta2_) result(ret)
  integer :: ind,ret
  real*8 :: eta2_
  
  ETA2(ind)=eta2_
  ret=0
end function
function set_node_vx(ind,uu2_) result(ret)
  integer :: ind,ret
  real*8 :: uu2_
  
  UU2(ind)=uu2_
  ret=0
end function
function set_node_vy(ind,vv2_) result(ret)
  integer :: ind,ret
  real*8 :: vv2_
  
  VV2(ind)=vv2_
  ret=0
end function
function set_node_eta_prev(ind,v1) result(ret)
  integer :: ind,ret
  real*8 :: v1  
  ETA1(ind)=v1
  ret=0
end function
function set_node_deta_dt(ind,x) result(ret)
  integer :: ind,ret
  real*8 :: x
  DETA_DT(ind)=x
  ret=0
end function
function set_node_status(ind,v1) result(ret)
  integer :: ind,ret
  character(len=*) :: v1
  if(v1.EQ."wet") NNODECODE(ind)=1
  if(v1.EQ."dry") NNODECODE(ind)=0
  ret=0
end function


function get_node_depth(ind,dp_) result(ret)
  integer :: ind,ret
  real*8 :: dp_
  if(ind.LT.1.OR.ind.GT.NP) then
    ret=-1
    return
  endif
  dp_=DP(ind)
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

function get_node_wave_stress(ind,wsx_,wsy_) result(ret)
  integer :: ind,ret
  real*8 :: wsx_,wsy_
  wsx_=RSNX(ind)*RhoWat0
  wsy_=RSNY(ind)*RhoWat0
  ret=0
end function
function set_node_wave_stress(ind,wsx_,wsy_) result(ret)
  integer :: ind,ret
  real*8 :: wsx_,wsy_
  RSNX(ind)=wsx_/RhoWat0
  RSNY(ind)=wsy_/RhoWat0
  ret=0
end function

function get_node_tidal_potential(ind,x_) result(ret)
  integer :: ind,ret
  real*8 :: x_
  x_=G*TIP(ind)
  ret=0
end function
function set_node_tidal_potential(ind,x_) result(ret)
  integer :: ind,ret
  real*8 :: x_
  TIP(ind)=x_/G
  ret=0
end function

function get_node_atmospheric_pressure(ind,x) result(ret)
  integer :: ind,ret
  real*8 :: x
  x=PR(ind)+reference_pressure
  ret=0
end function
function set_node_atmospheric_pressure(ind,x) result(ret)
  integer :: ind,ret
  real*8 :: x
  PR(ind)=x-reference_pressure
  ret=0
end function

function get_node_surface_heat_flux(ind,hflux_) result(ret)
  integer :: ind,ret
  real*8 :: hflux_
  hflux_=AMUSE_HFLUX(ind)
  ret=0
end function
function set_node_surface_heat_flux(ind,hflux_) result(ret)
  integer :: ind,ret
  real*8 :: hflux_
  AMUSE_HFLUX(ind)=hflux_
  ret=0
end function

function get_node_sigma(ind,indz,s_,z_) result(ret)
  integer :: ind,indz,ret
  real*8 :: s_,z_
  real*8, parameter :: a=1,b=-1
  if(ind.LT.1.OR.ind.GT.NP) then
    ret=-1
    return
  endif
  if(indz.LT.1.OR.indz.GT.NFEN) then
    ret=-1
    return
  endif
  s_=SIGMA(indz)
  z_=(s_-a)/(a-b)*(DP(ind)+ETA2(ind))+ETA2(ind)
  ret=0
end function

function get_node_velocities_3d(ind,indz,vx_,vy_,vz_) result(ret)
  integer :: ind,indz,ret
  real*8 :: vx_,vy_,vz_
  if(ind.LT.1.OR.ind.GT.NP) then
    ret=-1
    return
  endif
  if(indz.LT.1.OR.indz.GT.NFEN) then
    ret=-1
    return
  endif
  vx_=REAL(Q(ind,indz))
  vy_=IMAG(Q(ind,indz))
  vz_=WZ(ind,indz)
  ret=0
end function

function get_node_temperature_3d(ind,indz,temp_) result(ret)
  integer :: ind,indz,ret
  real*8 :: temp_
  if(ind.LT.1.OR.ind.GT.NP) then
    ret=-1
    return
  endif
  if(indz.LT.1.OR.indz.GT.NFEN) then
    ret=-1
    return
  endif
  temp_=TEMP(ind,indz)
  ret=0
end function

function get_element_nodes(ind,n1,n2,n3) result(ret)
  integer :: ind,n1,n2,n3,ret
  if(ind.LT.1.OR.ind.GT.NE) then
    ret=-1
    return
  endif
  n1=NM(ind,1)-1
  n2=NM(ind,2)-1
  n3=NM(ind,3)-1
  ret=0
end function

function get_element_status(ind,x) result(ret)
  integer :: ind,ret
  character(len=*) :: x
  if(NOFF(ind).EQ.0) then
    x="dry"
  else
    x="wet"
  endif
  ret=0
end function

function set_element_status(ind,v1) result(ret)
  integer :: ind,ret
  character(len=*) :: v1
  if(v1.EQ."wet") NOFF(ind)=1
  if(v1.EQ."dry") NOFF(ind)=0
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
function get_element_coordinates(ind,lon_,lat_) result(ret)
  integer :: ind,n1,n2,n3,ret
  real*8 :: lon_,lat_

  if(ind.LT.1.OR.ind.GT.NE) then
    ret=-1
    return
  endif
  n1=NM(ind,1)
  n2=NM(ind,2)
  n3=NM(ind,3)  
  lon_=(SLAM(n1)+SLAM(n2)+SLAM(n3))/3
  lat_=(SFEA(n1)+SFEA(n2)+SFEA(n3))/3
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
  index_node=NBDV(iseg,inode)-1
  ret=0
end function

function get_elevation_boundary_eta(inode,iseg,eta) result(ret)
  integer :: ret, inode,iseg,i  
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
  integer :: ret, inode,iseg,i  
  real*8 :: eta
  i=sum(NVDLL(1:iseg-1))+inode
  if(NBD(i).NE.NBDV(iseg,inode)) then
    ret=-1
    return
  endif
  ESBIN(i)=eta
  ret=0
end function

function get_flow_boundary_fluxx(inode,iseg,vx) result(ret)
  integer :: ret, inode,iseg,i  
  real*8 :: vx
  i=sum(NVELL(1:iseg-1))+inode
  if(NBV(i).NE.NBVV(iseg,inode)) then
    ret=-1
    return
  endif
  vx=FBQX(i)
  ret=0
end function

function set_flow_boundary_fluxx(inode,vx,iseg) result(ret)
  integer :: ret, inode,iseg,i  
  real*8 :: vx
  i=sum(NVELL(1:iseg-1))+inode
  if(NBV(i).NE.NBVV(iseg,inode)) then
    ret=-1
    return
  endif
  if(LBCODEI(i).NE.2.AND.LBCODEI(i).NE.12.AND.LBCODEI(i).NE.22) then
    ret=-2
    return
  endif
  FBQX(i)=vx
  ret=0
end function

function get_flow_boundary_fluxy(inode,iseg,vy) result(ret)
  integer :: ret, inode,iseg,i  
  real*8 :: vy
  i=sum(NVELL(1:iseg-1))+inode
  if(NBV(i).NE.NBVV(iseg,inode)) then
    ret=-1
    return
  endif
  vy=FBQY(i)
  ret=0
end function

function set_flow_boundary_fluxy(inode,vy,iseg) result(ret)
  integer :: ret, inode,iseg,i
  real*8 :: vy
  i=sum(NVELL(1:iseg-1))+inode
  if(NBV(i).NE.NBVV(iseg,inode)) then
    ret=-1
    return
  endif
  if(LBCODEI(i).NE.2.AND.LBCODEI(i).NE.12.AND.LBCODEI(i).NE.22) then
    ret=-2
    return
  endif
  FBQY(i)=vy
  ret=0
end function

function get_boundary_lnm(inode,iseg,x_) result(ret)
  integer :: ret, inode,iseg,i  
  real*8 :: x_
  i=sum(NVDLL(1:iseg-1))+inode
  if(NBD(i).NE.NBDV(iseg,inode)) then
    ret=-1
    return
  endif
  x_=AMUSE_LNM(i)
  ret=0
end function

function set_boundary_lnm(inode,x_,iseg) result(ret)
  integer :: ret, inode,iseg,i  
  real*8 :: x_
  i=sum(NVDLL(1:iseg-1))+inode
  if(NBD(i).NE.NBDV(iseg,inode)) then
    ret=-1
    return
  endif
 AMUSE_LNM(i)=x_
  ret=0
end function

function get_boundary_salinity(inode,zind,iseg,x_) result(ret)
  integer :: ret, inode,iseg,i,zind  
  real*8 :: x_
  i=sum(NVDLL(1:iseg-1))+inode
  if(NBD(i).NE.NBDV(iseg,inode)) then
    ret=-1
    return
  endif
  if(zind.LT.1.OR.zind.GT.NFEN) then
    ret=-2
    return
  endif
  x_=AMUSE_RESSAL(i,zind)
  ret=0
end function

function set_boundary_salinity(inode,zind,x_,iseg) result(ret)
  integer :: ret, inode,iseg,i,zind  
  real*8 :: x_
  i=sum(NVDLL(1:iseg-1))+inode
  if(NBD(i).NE.NBDV(iseg,inode)) then
    ret=-1
    return
  endif
  if(zind.LT.1.OR.zind.GT.NFEN) then
    ret=-2
    return
  endif
  AMUSE_RESSAL(i,zind)=x_
  ret=0
end function

function get_boundary_temperature(inode,zind,iseg,x_) result(ret)
  integer :: ret, inode,iseg,i,zind  
  real*8 :: x_
  i=sum(NVDLL(1:iseg-1))+inode
  if(NBD(i).NE.NBDV(iseg,inode)) then
    ret=-1
    return
  endif
  if(zind.LT.1.OR.zind.GT.NFEN) then
    ret=-2
    return
  endif
  x_=AMUSE_RESTEMP(i,zind)
  ret=0
end function

function set_boundary_temperature(inode,zind,x_,iseg) result(ret)
  integer :: ret, inode,iseg,i,zind  
  real*8 :: x_
  i=sum(NVDLL(1:iseg-1))+inode
  if(NBD(i).NE.NBDV(iseg,inode)) then
    ret=-1
    return
  endif
  if(zind.LT.1.OR.zind.GT.NFEN) then
    ret=-2
    return
  endif
  AMUSE_RESTEMP(i,zind)=x_
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
  index_node=NBVV(iseg,inode)-1
  ret=0
end function

function get_flow_boundary_type(inode,iseg,type_) result(ret)
  integer :: ret, inode,iseg,type_  
  if(IBTYPE(iseg).NE.LBCODEI(NBVV(iseg,inode))) then
   ret=-1
   return
  endif
  type_=IBTYPE(iseg)
  ret=0
end function

function get_number_of_vertical_nodes(nout) result(ret)
  integer :: ret,nout
  nout=NFEN
  ret=0  
end function

function get_reference_pressure(x) result(ret)
  integer :: ret
  real*8 :: x
  x=reference_pressure
  ret=0
end function

end module


