module amuse_adcirc
  use sizes
  use GLOBAL

  logical :: use_interface_elevation_boundary=.FALSE.
  
  real(SZ),allocatable ::   ESBIN(:)
  
#ifdef CSWAN
#error SWAN coupling TBD 
#endif

contains

function initialize_code() result(ret)
  USE ADCIRC_Mod, ONLY : ADCIRC_Init
  integer :: ret

  CALL ADCIRC_Init()
  
  allocate( ESBIN(MNETA) )
  ESBIN(:)=0.
  
  ret=0
end function

function get_model_time(tout) result(ret)
  USE ADCIRC_Mod

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


function evolve_model(tend) result(ret)
  USE ADCIRC_Mod
  integer :: ret
  real(8) :: tend

  do while(ITIME_BGN*DTDP+STATIM*86400.D0<tend-DTDP/2)
    call ADCIRC_Run(1)
  enddo

  if(.not.use_interface_elevation_boundary) ESBIN(1:NETA)=ETA2(NBD(1:NETA))

  ret=0
end function

function cleanup_code() result(ret)
  USE ADCIRC_Mod, ONLY : ADCIRC_Final

  integer :: ret

  call ADCIRC_Final
  ret=0
  
end function

function commit_parameters() result(ret)
 integer :: ret
 ret=-2
end function

function recommit_parameters() result(ret)
 integer :: ret
 ret=-2
end function

subroutine update_elevation_boundary()
  ETA2(NBD(1:NETA))=ESBIN(1:NETA)
end subroutine

end module


subroutine AMUSE_elevation_boundary() 
  use amuse_adcirc
  if(use_interface_elevation_boundary) call update_elevation_boundary()
end subroutine

