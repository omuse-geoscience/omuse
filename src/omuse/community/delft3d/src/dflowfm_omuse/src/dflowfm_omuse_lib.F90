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

  use hashMod

  implicit none

  type(hash_type) :: hash
  integer :: ndx_loc, ndxi_loc, ndxi_glob, ndx_glob ! not available in dflowfm !?
  integer, allocatable :: gid(:) ! gloabl id for local cells only

contains

  function initialize_index_map() result(ierr)
     integer :: ierr,i
  
     if(ndxi.NE.ndx2d.OR.ndx1Db.NE.ndxi) then
       print*, "err:", ndxi,ndx2d,ndx,ndx1Db
       ierr=-1 ! assumption broken
       return
     endif

     allocate(gid(ndx))

     if(jampi.EQ.0) then
       if(allocated(iglobal_s)) then
         ierr=-2 ! unexpectedly allocated
         return
       endif

       do i=1,ndx
         gid(i)=i
       enddo
       ndxi_loc=ndxi
       ndxi_glob=ndxi
       ndx_loc=ndx
       ndx_glob=ndx
       
     else
 
       gid=-1

       ndx_loc=0
       ndxi_loc=0
       do i=1,ndxi
         if(idomain(i).EQ.my_rank) then
           ndxi_loc=ndxi_loc+1
         endif           
       enddo
       do i=1,ndx
         if(idomain(i).EQ.my_rank) then
           ndx_loc=ndx_loc+1
           gid(i)=iglobal_s(i)
         endif
       enddo

       call reduce_int_sum(ndxi_loc, ndxi_glob)
       call reduce_int_sum(ndx_loc, ndx_glob)

     endif
  
     call initHash(ndx/2+1, ndx, gid, hash)
    
     ierr=0
     
  end function

  function find_index(global_index) result(local_index)
      integer :: local_index, global_index
      
      local_index=find(global_index, gid, hash)
      
  end function 


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
  
!~      print*, "info:"
!~      print*, ndx, ndx2d, ndkx
!~      print*, ndxi,ndx1Db
!~      print*, ndx2dr
!~      print*, Nglobal_s, nump1d2d
!~      print*, idomain(:5)
!~      print*, iglobal_s(:5),size(iglobal_s)
!~      print*, ndomains
!~      print*, numk
     
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

  function get_x_position_(i, x, n) result (ret)
    integer :: ret,i(n),n,i_, i__
    double precision :: x(n)
    ret=0
    do i_=1,n
      i__=find_index(i(i_))
      if(i__.GT.0) then
        if(gid(i__).NE.i(i_)) ret=-1 
        x(i_)=xz(i__)
      else
        x(i_)=1
      endif
    enddo
    ret=0
#ifdef HAVE_MPI
    call mpi_allreduce(mpi_in_place,x,n,mpi_double_precision,mpi_sum,DFM_COMM_DFMWORLD,ret)
#endif
  end function

  function get_y_position_(i, x,n) result (ret)
    integer :: ret,i(n),n,i_, i__
    double precision :: x(n)
    do i_=1,n
      i__=find_index(i(i_))
      if(i__.NE.0) then
        if(gid(i__).NE.i(i_)) ret=-1 
        x(i_)=yz(i__)
      else
        x(i_)=0
      endif
    enddo
    ret=0
#ifdef HAVE_MPI
    call mpi_allreduce(mpi_in_place,x,n,mpi_double_precision,mpi_sum,DFM_COMM_DFMWORLD,ret)
#endif
  end function

  function get_water_level_(i, x,n) result (ret)
    integer :: ret,i(n),n,i_, i__
    double precision :: x(n)
    do i_=1,n
      i__=find_index(i(i_))
      if(i__.NE.0) then
        if(gid(i__).NE.i(i_)) ret=-1 
        x(i_)=s1(i__)
      else
        x(i_)=0
      endif
    enddo
    ret=0
#ifdef HAVE_MPI
    call mpi_allreduce(mpi_in_place,x,n,mpi_double_precision,mpi_sum,DFM_COMM_DFMWORLD,ret)
#endif
  end function

end module
