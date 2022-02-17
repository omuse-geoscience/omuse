module partitioner_omuse
  use dflowfm_omuse_lib
!~    USE M_GRID
!~    USE M_POLYGON
!~    USE M_LANDBOUNDARY
!~    USE M_BOAT
!~    use m_netw
!~    use unstruc_startup
!~    use unstruc_model
!~    use netcdf
!~    use properties
!~    use m_observations
!~    use unstruc_netcdf
!~    use unstruc_messages
!~    USE UNSTRUC_DISPLAY
!~    USE M_WEARELT
!~    use m_flowparameters
  implicit none

   integer :: MODE,NFLD, KEY
   integer :: JQN
   integer :: JDEMO

   COMMON /MODENOW/ MODE,NFLD
   COMMON /QNRGF/ JQN
   COMMON /DEMO/ JDEMO





  character(len=255) :: fnam             !< filename

contains

  include "partitioner_getter_setters.f90"

  function initialize_code() result(ret)
    integer :: ret
    fnam="omuse"
    md_Ndomains=3
    md_jacontiguous=0
    md_icgsolver=6
    md_pmethod=0
    md_dryptsfile=""
    md_encfile=""
    md_genpolygon=1
    
    japartdomain=1

    ret=0
  end function

  function cleanup_code() result(ret)
    integer :: ret
    ret=0
  end function

  function commit_parameters() result(ret)
    integer :: ret

    ret=initialize_dflowfm_core(.TRUE.)

   JDEMO        = 0
   JQN          = 2

   MMAX   = 0
   NMAX   = MMAX
   KMAX   = MMAX*NMAX
   KNX    = 8
   MXB    = 10
   LMAX   = (MMAX-1)*NMAX + (NMAX-1)*MMAX + 2*(MMAX-1)*(NMAX-1)
   MAXLAN = 500
   MAXPOL = MAXLAN
   MAXBOAT = MAXLAN

   call resetFullFlowModel()
!~    call INIDAT()

  end function

  function recommit_parameters() result(ret)
    integer :: ret
    ret=-1
  end function

  function do_partition() result(ret)
    integer :: ret, istat
    
    call loadNetwork(fnam, istat,1)
    
    call partition_from_commandline(fnam, md_Ndomains, md_jacontiguous, md_icgsolver, &
        md_pmethod, md_dryptsfile, md_encfile, md_genpolygon)

    ret=0
  
  end function

end module
