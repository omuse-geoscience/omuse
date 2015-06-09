module amuse_swan
  USE TIMECOMM                                                        
  USE OCPCOMM1                                                        
  USE OCPCOMM2                                                        
  USE OCPCOMM3                                                        
  USE OCPCOMM4                                                        
  USE SWCOMM1                                                         
  USE SWCOMM2                                                         
  USE SWCOMM3                                                         
  USE SWCOMM4                                                         
  USE OUTP_DATA                                                       
  USE M_SNL4                                                          
  USE M_GENARR                                                        
  USE M_OBSTA                                                         
  USE M_PARALL                                                        
  USE SwanGriddata                                                    

contains

function swan_entry() result(ret)
  integer :: ret
  logical STPNOW 
  ret=0
  call SWINITMPI
  CALL RDINIT
  IF (STPNOW()) ret=-1
!PUN      CALL MSG_INIT()                                             

end function

function swan_cleanup() result(ret)
  integer :: ret
  CALL SWEXITMPI
!PUN      CALL MSG_FINI()
end function

function swan_init() result(ret)
  integer :: ret
  logical STPNOW 
  ret=0
  call SWINIT(INERR)
!PUN      IF ( MNPROC>1 ) THEN
!PUN         CALL SwanReadfort18
!PUN!NADC         CALL MSG_TABLE()
!PUN!NADC         CALL MSG_START()
!PUN      ENDIF
!TIMG      CALL SWTSTO(2)
  if(INERR.GT.0) ret=-1
  if(STPNOW()) ret=-2
end function

end module
