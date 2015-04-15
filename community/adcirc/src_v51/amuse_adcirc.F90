module amuse_adcirc
  use sizes
  use GLOBAL
  use ADCIRC_Mod
  use BOUNDARIES, only: NETA, NBD, NVDLL, NBDV, NOPE, NBOU, NVELL, NBVV, IBTYPE
  use MESH, only: NP, NE, X, Y
  implicit none

  logical :: use_interface_elevation_boundary=.FALSE.
  logical :: use_interface_met_forcing=.FALSE.
  
  real(SZ),allocatable ::   ESBIN(:), WSX(:),WSY(:)
  
#ifdef CSWAN
#error SWAN coupling TBD 
#endif

contains

subroutine update_elevation_boundary()
  ETA2(NBD(1:NETA))=ESBIN(1:NETA)
end subroutine

subroutine update_met_forcing()
  WSX2(1:NP)=WSX(1:NP)
  WSY2(1:NP)=WSY(1:NP)
  NWS=OR(NWS,2**30)
end subroutine

end module


subroutine AMUSE_elevation_boundary() 
  use amuse_adcirc
  if(use_interface_elevation_boundary) call update_elevation_boundary()
end subroutine

subroutine AMUSE_met_forcing()
  use amuse_adcirc
  if(use_interface_met_forcing) call update_met_forcing()
end subroutine
