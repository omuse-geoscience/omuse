module amuse_adcirc
  use GLOBAL

#ifdef CSWAN
#err SWAN coupling TBD 
#endif

#err test this

contains

function initialize_code() result(ret)
  USE ADCIRC_Mod, ONLY : ADCIRC_Init
  integer :: ret

  CALL ADCIRC_Init
  
  ret=0
end function

function evolve_model(tend) result(ret)
  USE ADCIRC_Mod, ONLY : ADCIRC_Run
  integer :: ret

  call ADCIRC_Run

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
