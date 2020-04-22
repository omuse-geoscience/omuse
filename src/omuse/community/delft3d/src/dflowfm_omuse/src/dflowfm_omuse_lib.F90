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
  integer, allocatable :: gid(:) ! global id for local cells only (=iglobal_s)

  type(hash_type) :: lhash
  integer :: lnxi_loc, lnx_loc, lnxi_glob, lnx_glob ! not available in dflowfm !?
  integer, allocatable :: lgid(:) ! global id for local links only

contains

  ! this is mirrors get_global_numbers, except no exchange of ghost info
  subroutine get_global_numbering(n, iglobnum, offset)
    integer, intent(in) :: n !< size of iglobnum
    integer, dimension(:), intent(inout) :: iglobnum  !< input: array with 0 if not count or 1 if count output: array with global numbers. 
    integer, optional, intent(in) :: offset

    integer                                     :: i, num
    integer                                     :: ierror

    integer, dimension(:), allocatable :: numglobcells
   
    ierror = 1

    allocate(numglobcells(0:max(ndomains-1,0)))
            
    if ( jampi.eq.1 ) then
         
      num = count(iglobnum(1:n).eq.1)
      
#ifdef HAVE_MPI
      call mpi_allgather(num, 1, MPI_INTEGER, numglobcells, 1, MPI_INTEGER, DFM_COMM_DFMWORLD, ierror)
#else 
      numglobcells(0)=num
#endif

      num = 0
      if ( my_rank.gt.0 ) then
        num = sum(numglobcells(0:my_rank-1))
      end if
    else  ! jampi.eq.0
      numglobcells(0) = count(iglobnum(1:n).eq.1)
      num = 0
    end if 

    if(present(offset)) num=num+offset

    ! renumber
    do i=1,n
      if ( iglobnum(i).ne.0 ) then
        num = num+1
        iglobnum(i) = num
      end if
    end do

    ierror = 0
  1234 continue

    if ( allocated(numglobcells) ) deallocate(numglobcells)
      
  end subroutine get_global_numbering


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
         endif
       gid(i)=iglobal_s(i)
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

  function initialize_link_map() result(ierr)
     integer :: ierr,i, j, l1,l2, jaghost, idmn
     double precision, dimension(:), allocatable :: dum, dum2
      
     if(lnx1D.NE.0.OR.lnx1db.NE.lnxi) then
       print*, "err:", lnx1d,lnxi,lnx1db,lnx
       ierr=-1 ! assumption broken
       return
     endif

     allocate(lgid(lnx))

     if(jampi.EQ.0) then
       do i=1,lnx
         lgid(i)=i
       enddo
       lnxi_loc=lnxi
       lnxi_glob=lnxi
       lnx_loc=lnx
       lnx_glob=lnx
       
     else
 
       lgid=0

       lnxi_loc=0
!~        lnx_loc=0

!~        do i=1,ndx
!~          if(idomain(i).NE.my_rank) cycle
!~          do j=1,nd(i)%lnx
!~            if(nd(i)%ln(j).GT.0) then
!~              lnx_loc=lnx_loc+1
!~              if(nd(i)%ln(j).LE.lnxi) lnxi_loc=lnxi_loc+1
!~            endif
!~          enddo
!~        enddo

!~        do i=1,lnxi
!~          l1=lne(1,i)
!~          l2=lne(2,i)
!~          if(idomain(l1).EQ.my_rank.OR.idomain(l2).EQ.my_rank) then
!~            if(idomain(l1).GE.idomain(l2)) then
!~              lnxi_loc=lnxi_loc+1
!~            endif
!~          endif           
!~        enddo
!~        lnx_loc=lnxi_loc
!~        do i=lnxi+1,lnx
!~          l1=lne(1,i)
!~          l2=lne(2,i)
!~          if(idomain(l1).EQ.my_rank.OR.idomain(l2).EQ.my_rank) then
!~            if(idomain(l1).GE.idomain(l2)) then
!~              lnx_loc=lnx_loc+1
!~            endif
!~          endif           
!~        enddo

! this seems to work
       do i=1,lnxi
         l1=ln(1,i)
         l2=ln(2,i)
         call link_ghostdata(my_rank,idomain(l1),idomain(l2), jaghost, idmn) 
         if(jaghost.EQ.0.AND.idmn.EQ.my_rank) then
             lnxi_loc=lnxi_loc+1
             lgid(i)=1
         endif           
       enddo
       lnx_loc=lnxi_loc
       do i=lnxi+1,lnx
         l1=ln(1,i)
         l2=ln(2,i)
         call link_ghostdata(my_rank,idomain(l1),idomain(l2), jaghost, idmn) 
         if(jaghost.EQ.0.AND.idmn.EQ.my_rank) then
             lnx_loc=lnx_loc+1
             lgid(i)=1
         endif           
       enddo

       call reduce_int_sum(lnxi_loc, lnxi_glob)
       call reduce_int_sum(lnx_loc, lnx_glob)
       
       print*, my_rank, lnxi_loc, lnxi_glob,lnx_loc, lnx_glob

       call get_global_numbering(lnxi, lgid)
       if(lnx_glob.GT.lnxi_glob) call get_global_numbering(lnx-lnxi, lgid(lnxi+1:lnx), lnxi_glob)
      
       allocate(dum(lnx))
      
       dum=dble(lgid)

!~        allocate(dum2(lnx))
!~        where(lgid.NE.0) dum2=1.

       call update_ghosts(ITYPE_U, 1, lnx, dum, ierr)
       lgid=int(dum)

!~        call update_ghosts(ITYPE_U, 1, lnx, dum2, ierr)
!~        print*, ">",my_rank, minval(dum2), maxval(dum2)


       deallocate(dum)
!~        deallocate(dum2)

     endif
  
     call initHash(lnx/2+1, lnx, lgid, lhash)
    
     ierr=0


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
  
     print*, "info:"
     print*, my_rank, numk, numL1D, numL
     print*, my_rank, nump, ndx2d, ndxi, ndx1Db, ndx
     print*, my_rank, lnx1D, lnxi, lnx1D, lnx

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
    integer :: ret_,ret,i(n),n,i_, i__
    double precision :: x(n)
    ret=0
    ret_=0
    do i_=1,n
      x(i_)=0
      i__=find_index(i(i_))
      if(i__.GT.0) then
        if(gid(i__).NE.i(i_)) ret=-1 
        if(idomain(i__).EQ.my_rank) x(i_)=xz(i__)
      endif
    enddo
#ifdef HAVE_MPI
    call mpi_allreduce(mpi_in_place,x,n,mpi_double_precision,mpi_sum,DFM_COMM_DFMWORLD,ret)
    call reduce_int_sum(ret, ret_)
    ret=ret_
#endif
  end function

  function get_y_position_(i, x,n) result (ret)
    integer :: ret_,ret,i(n),n,i_, i__
    double precision :: x(n)
    ret=0
    ret_=0
    do i_=1,n
      x(i_)=0
      i__=find_index(i(i_))
      if(i__.GT.0) then
        if(gid(i__).NE.i(i_)) ret=-1 
        if(idomain(i__).EQ.my_rank) x(i_)=yz(i__)
      endif
    enddo
#ifdef HAVE_MPI
    call mpi_allreduce(mpi_in_place,x,n,mpi_double_precision,mpi_sum,DFM_COMM_DFMWORLD,ret)
    call reduce_int_sum(ret, ret_)
    ret=ret_
#endif
  end function

  function get_water_level_(i, x,n) result (ret)
    integer :: ret_,ret,i(n),n,i_, i__
    double precision :: x(n)
    ret=0
    ret_=0
    do i_=1,n
      x(i_)=0
      i__=find_index(i(i_))
      if(i__.GT.0) then
        if(gid(i__).NE.i(i_)) ret=-1 
        if(idomain(i__).EQ.my_rank) x(i_)=s1(i__)
      endif
    enddo
#ifdef HAVE_MPI
    call mpi_allreduce(mpi_in_place,x,n,mpi_double_precision,mpi_sum,DFM_COMM_DFMWORLD,ret)
    call reduce_int_sum(ret, ret_)
    ret=ret_
#endif
  end function

end module
