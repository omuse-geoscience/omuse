module amuse_adcirc
  use sizes
  use GLOBAL
  use ADCIRC_Mod

  logical :: use_interface_elevation_boundary=.FALSE.
  
  real(SZ),allocatable ::   ESBIN(:)
  
#ifdef CSWAN
#error SWAN coupling TBD 
#endif

contains

subroutine update_elevation_boundary()
  ETA2(NBD(1:NETA))=ESBIN(1:NETA)
end subroutine

end module


subroutine AMUSE_elevation_boundary() 
  use amuse_adcirc
  if(use_interface_elevation_boundary) call update_elevation_boundary()
end subroutine

