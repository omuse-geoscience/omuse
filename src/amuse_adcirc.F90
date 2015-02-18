module amuse_adcirc
  use sizes
  use GLOBAL

#ifdef CSWAN
#error SWAN coupling TBD 
#endif

contains

function initialize_code() result(ret)
  USE ADCIRC_Mod, ONLY : ADCIRC_Init
  integer :: ret

  CALL ADCIRC_Init
  
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

  ret=0
end function

function cleanup_code() result(ret)
  USE ADCIRC_Mod, ONLY : ADCIRC_Final

  integer :: ret

  call ADCIRC_Final

  ret=0
  
end function

function test() result(ret)
 integer :: ret

 print*," G= ",g
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



end module
