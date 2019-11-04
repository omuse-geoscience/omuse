#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

module dflowfm_omuse_lib
  use unstruc_api
  use unstruc_display, only: jaGUI ! this should be removed when jaGUI = 0 by default

  use m_partitioninfo
  use m_flow
  use m_waves
  use m_ship
  use m_wind
  use m_flowgeom
  use network_data
  use m_flowexternalforcings
  use m_partitioninfo
  use m_sobekdfm
  use m_transport
  use m_xbeach_data
  use dfm_error
  use m_vegetation
  use m_sediment
  use m_integralstats
  use gridoperations
  use unstruc_model
  use unstruc_files
#ifdef HAVE_MPI
  use mpi
#endif

  implicit none

contains

  function initialize_dflowfm(config_file) result(ret)    
     character(len=*) :: config_file
  
     ! Extra local variables
     integer :: inerr, ret  ! number of the initialisation error
     logical :: mpi_initd
  
  !~    integer(c_int), target, allocatable, save :: x(:,:)
  !~    type(c_ptr) :: xptr
     integer :: i,j,k
      
     ret = 0 
#ifdef HAVE_MPI
     call mpi_initialized(mpi_initd, inerr)
     if (.not. mpi_initd) then
        ja_mpi_init_by_fm = 1
        call mpi_init(inerr)
     else
        ja_mpi_init_by_fm = 0
     end if
  
     call mpi_comm_rank(DFM_COMM_DFMWORLD,my_rank,inerr)
     call mpi_comm_size(DFM_COMM_DFMWORLD,numranks,inerr)
  
     if ( numranks.le.1 ) then
        jampi = 0
     end if
  
     !   make domain number string as soon as possible
     write(sdmn, '(I4.4)') my_rank
  
#else
     numranks=1
#endif
  
     ! do this until default has changed
     jaGUI = 0
  
     ! TODO: check why these are needed to avoid a segfault
     KNX    = 8
     MXB    = 10
     MAXLAN = 500
     MAXPOL = MAXLAN
  
     call init_core()
  
     CALL INIDAT()
     call api_loadmodel(config_file)
  
     !PETSC must be called AFTER reading the mdu file, so the icgsolver option is known to startpetsc 
#ifdef HAVE_PETSC
     call startpetsc()
#endif
  
     ret = flowinit()
  
     time_user = tstart_user
  
     ! Just terminate if we get an error....
     if (inerr > 0) ret=-1
  
  end function 
  
  function evolve_model(tend) result(ret)
!~      use messagehandling
     double precision, intent(in) :: tend
     double precision :: dt
  
     integer ::  ierr, ret

     ! bmi comments:
     ! The time loop seems to be located in unstruc->flow
     ! It is important that we can simulate up to a time set from the outside
     ! We might have to set time_user or dt_user
     
     dt=tend-time1
  
     call flow_run_sometimesteps(dt, ierr)
  
     ret = ierr
  end function 
  
  function get_current_time(t) result(ret)
     integer :: ret
     double precision, intent(out) :: t
  
     t = time1
     ret = 0
  end function 

end module
