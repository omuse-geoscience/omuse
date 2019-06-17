module amuse_adcirc
  use sizes
  use GLOBAL
  use ADCIRC_Mod
  use BOUNDARIES, only: NETA, NBD, NVDLL, NBDV, NOPE, NBOU, NVELL, NBVV, IBTYPE, &
    NBV, NVEL, CSII,SIII, LBCODEI
  use MESH, only: NP, NE, X, Y,DP,SLAM,SFEA
  use GLOBAL_3DVS, only: NFEN,SIGMA,Q,WZ, RESSAL, RESTEMP, qsurfkp1, qsurf, SigT, Temp, Sal
  implicit none

  logical :: use_interface_elevation_boundary=.FALSE.
  logical :: use_interface_flow_boundary=.FALSE.
  logical :: use_interface_met_forcing=.FALSE.
  logical :: use_interface_wave_forcing=.FALSE.
  logical :: use_interface_tidal_forcing=.FALSE.
  
  logical :: use_interface_lnm_boundary=.FALSE.
  logical :: use_interface_salinity_boundary=.FALSE.
  logical :: use_interface_temperature_boundary=.FALSE.
  logical :: use_interface_surface_heat_forcing=.FALSE.
  
  real(SZ),allocatable ::   ESBIN(:), WSX(:),WSY(:), DETA_DT(:), PR(:), &
                            RSNX(:), RSNY(:), TIP(:), FBQX(:), FBQY(:)
  
  real(SZ),allocatable ::   AMUSE_LNM(:), AMUSE_RESSAL(:,:), AMUSE_RESTEMP(:,:), &
                            AMUSE_HFLUX(:)
  
  real(SZ),parameter :: reference_pressure=1013.0d0

contains

subroutine update_elevation_boundary(rampfactor)
  real(SZ) :: rampfactor
  ETA2(NBD(1:NETA))=rampfactor*ESBIN(1:NETA)
end subroutine

subroutine update_flow_boundary(rampfactor)
  real(SZ) :: rampfactor
! note the minus 
  QN2(NBV(1:NVEL))=-rampfactor* (CSII(1:NVEL)*FBQX(1:NVEL)+SIII(1:NVEL)*FBQY(1:NVEL))
end subroutine

subroutine update_met_forcing(rampfactor)
  real(SZ) :: rampfactor
  WSX2(1:NP)=rampfactor*WSX(1:NP)
  WSY2(1:NP)=rampfactor*WSY(1:NP)
  PR2(1:NP)=100*(reference_pressure+rampfactor*PR(1:NP))/(RHOWAT0*G)
  NWS=OR(NWS,2**30)
end subroutine

subroutine update_wave_forcing(rampfactor)
  real(SZ) :: rampfactor  
  WSX2(1:NP)=WSX2(1:NP)+rampfactor*RSNX(1:NP)
  WSY2(1:NP)=WSY2(1:NP)+rampfactor*RSNY(1:NP)
end subroutine

subroutine update_tidal_forcing()
  TIP2(1:NP)=TIP2(1:NP)+TIP(1:NP)
end subroutine

subroutine update_lnm_boundary()
  LNM_BC(NBD(1:NETA))=AMUSE_LNM(1:NETA)
end subroutine

subroutine update_salt_boundary()
  RESSAL(1:NETA,1:NFEN)=AMUSE_RESSAL(1:NETA,1:NFEN)
end subroutine 

subroutine update_temp_boundary()
  RESTEMP(1:NETA,1:NFEN)=AMUSE_RESTEMP(1:NETA,1:NFEN)
end subroutine 

subroutine update_surface_heat_forcing()
  qsurf(1:NP)=qsurfkp1(1:NP) ! store old value..?
  HFLUX(1:NP)=AMUSE_HFLUX(1:NP)
end subroutine

end module

subroutine AMUSE_elevation_boundary(rampfactor) 
  use amuse_adcirc
  REAL(SZ) :: rampfactor
  if(use_interface_elevation_boundary) call update_elevation_boundary(rampfactor)
end subroutine

subroutine AMUSE_flow_boundary(rampfactor) 
  use amuse_adcirc
  REAL(SZ) :: rampfactor
  if(use_interface_flow_boundary) call update_flow_boundary(rampfactor)
end subroutine

subroutine AMUSE_met_forcing(rampfactor)
  use amuse_adcirc
  REAL(SZ) :: rampfactor
  if(use_interface_met_forcing) call update_met_forcing(rampfactor)
  if(use_interface_wave_forcing) call update_wave_forcing(rampfactor)
end subroutine

subroutine AMUSE_tidal_forcing() 
  use amuse_adcirc
  if(use_interface_tidal_forcing) call update_tidal_forcing()
end subroutine

subroutine AMUSE_lnm_boundary()
  use amuse_adcirc
  if(use_interface_lnm_boundary) call update_lnm_boundary()
end subroutine

subroutine AMUSE_salt_temp_boundary()
  use amuse_adcirc
  if(use_interface_salinity_boundary) call update_salt_boundary()
  if(use_interface_temperature_boundary) call update_temp_boundary()
end subroutine

subroutine AMUSE_surface_heat_forcing()
  use amuse_adcirc
  if(use_interface_surface_heat_forcing) call update_surface_heat_forcing()
end subroutine
