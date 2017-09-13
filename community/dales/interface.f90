! This is the omuse-specific native interface to the 
! Dutch Large Eddy Simulation code DALES.
! Written by G. van den Oord, 2016

module dales_interface

    use daleslib, only: initialize,step,finalize,&
                       allocate_z_axis,allocate_2d,allocate_3d,&
                       my_task,master_task,&
                       FIELDID_U,FIELDID_V,FIELDID_W,FIELDID_THL,FIELDID_QT, &
                       u_tend, v_tend, thl_tend, qt_tend, ps_tend, &
                       gatherlayeravg, gathervol, localindex, gatherLWP, gatherCloudFrac
    use modfields, only: u0,v0,w0,thl0,qt0,ql0,e120,tmp0,sv0,um,vm,wm,thlm,qtm
    use modglobal, only: i1,j1,k1,itot,jtot,kmax,lmoist
    use modsurfdata, only: ps, qts
    use modsurface, only: qtsurf
    use modmicrodata, only: iqr
    
    !TODO: Expose everything so this module only depends on daleslib
    use modglobal, only: rtimee,rdt,fname_options,timeleft,tres,timee,rk3step,dt_lim
    use mpi, only: MPI_COMM_WORLD

    implicit none

    real(8),allocatable::tmp1Dz(:)
    real(8),allocatable::tmp2Dxy(:,:)
    real(8),allocatable::tmp3D(:,:,:)

    integer :: MPI_local_comm = MPI_COMM_WORLD
    
contains

    function set_input_file(ifile) result(ret)
        integer::                      ret
        character(256),intent(in)::    ifile

        fname_options=ifile
        ret=0
    end function set_input_file

    function get_input_file(ifile) result(ret)
        integer::                       ret
        character(256),intent(out)::    ifile

        ifile=fname_options
        ret=0
    end function get_input_file

      ! set the working directory of the process. 
      ! returns 0 on success
    function set_workdir(directory) result(ret)
#if defined (__INTEL_COMPILER)
      USE IFPORT   ! for intel chdir function.
#endif

      integer::                      ret
      character(256),intent(in)::   directory
      ret = chdir(directory)  
      write(*,*) "Dales worker changing directory to", directory, "status:", ret
    end function set_workdir
    
    
    function initialize_code() result(ret)
        integer:: ret
        
        fname_options="namoptions.001"

        call initialize(fname_options,MPI_local_comm)
        ret=0
    end function

    function commit_parameters() result(ret)
        integer:: ret

        ret=0
        if(my_task==master_task) then
            ret=allocate_z_axis(tmp1Dz)
            if(ret/=0) return
            ret=allocate_2d(tmp2dxy)
            if(ret/=0) return
            ret=allocate_3d(tmp3D)
        endif
    end function

    function recommit_parameters() result(ret)
        integer:: ret
        ret=0
    end function

    ! evolve the model to time tend
    ! returns 1 if the end of the simulation was reached (timeleft == 0)
    !
    ! EXPERIMENTAL: if exactEnd is nonzero,
    !               set dt_lim so that the timestepping finishes exactly at tend
    !               then the end time from namoptions is not used.
    !
    !  .or. rk3step < 3) SHOULD BE USED TO FINISH A FULL TIME STEP BEFORE EXITING ??
    !                    but it breaks the current result test, which was stored without this.
    !                    see program.f90: main loop - while(timeleft>0 .or. rk3step < 3)
    
    function evolve_model(tend, exactEnd) result(ret)
        integer             :: ret
        real(8),intent(in)  :: tend
        integer,intent(in)  :: exactEnd
        integer             :: i

        ! print *, 'evolve_model'
        do while((rtimee < tend .and. timeleft > 0)) ! .or. rk3step < 3)

           if (exactEnd /= 0) then
              ! set dt_lim to finish at tend (or before)
              dt_lim = min(dt_lim, (int(tend/tres) - timee))
           endif
        
           call step
!           print *, '  ', 'rtimee=',rtimee, 'timeleft=',timeleft, 'rk3step=',rk3step, 'rdt=',rdt           
        enddo
        ret=0
        if(timeleft <= 0) then
          ret=1
        endif
    end function

    function cleanup_code() result(ret)
        integer :: ret

        call finalize
        ret=0
    end function

    function get_timestep(dt) result(ret)
        real(8),intent(out):: dt
        integer            :: ret

        dt=rdt
        ret=0
    end function

    function get_model_time(t) result(ret)
        real(8),intent(out):: t
        integer::             ret

        t=rtimee
        ret=0
    end function get_model_time

    function set_surface_pressure(p) result(ret)
      real(8),intent(in)::  p
      integer::             ret

      ps = p

      ! modtimedep calls qtsurf to update surface values after changing ps, so we'll do the same
      if (lmoist) then
         call qtsurf
      else
         qts = 0.
      endif
      
      ret = 0
    end function set_surface_pressure

    function get_surface_pressure(p) result(ret)
      real(8),intent(out)::  p
      integer::             ret

      p = ps
      ret = 0
    end function get_surface_pressure
    

    function commit_grid() result(ret)
        integer::             ret
        ret=0
    end function


!!!!! get functions for vertical profiles / slab averages
    !   g_k is a dummy array
    !   a is the returned result
    !   n is the array length    
    function get_profile_U_(g_k, a, n) result(ret)
      integer, intent(in)                 :: n
      integer, dimension(n), intent(in)   :: g_k
      real, dimension(n), intent(out)     :: a
      integer                             :: ret

      ret = gatherlayeravg(u0(2:i1, 2:j1, :), a)
    end function get_profile_U_

    function get_profile_V_(g_k, a, n) result(ret)
      integer, intent(in)                 :: n
      integer, dimension(n), intent(in)   :: g_k
      real, dimension(n), intent(out)     :: a
      integer                             :: ret

      ret = gatherlayeravg(v0(2:i1, 2:j1, :), a)
    end function get_profile_V_
 
    function get_profile_W_(g_k, a, n) result(ret)
      integer, intent(in)                 :: n
      integer, dimension(n), intent(in)   :: g_k
      real, dimension(n), intent(out)     :: a
      integer                             :: ret

      ret = gatherlayeravg(w0(2:i1, 2:j1, :), a)
    end function get_profile_W_

    function get_profile_THL_(g_k, a, n) result(ret)
      integer, intent(in)                 :: n
      integer, dimension(n), intent(in)   :: g_k
      real, dimension(n), intent(out)     :: a
      integer                             :: ret

      ret = gatherlayeravg(thl0(2:i1, 2:j1, :), a)
    end function get_profile_THL_

    function get_profile_QT_(g_k, a, n) result(ret)
      integer, intent(in)                 :: n
      integer, dimension(n), intent(in)   :: g_k
      real, dimension(n), intent(out)     :: a
      integer                             :: ret

      ret = gatherlayeravg(qt0(2:i1, 2:j1, :), a)
    end function get_profile_QT_

    function get_profile_QL_(g_k, a, n) result(ret)
      integer, intent(in)                 :: n
      integer, dimension(n), intent(in)   :: g_k
      real, dimension(n), intent(out)     :: a
      integer                             :: ret

      ret = gatherlayeravg(ql0(2:i1, 2:j1, :), a)
    end function get_profile_QL_

    function get_profile_E12_(g_k, a, n) result(ret)
      integer, intent(in)                 :: n
      integer, dimension(n), intent(in)   :: g_k
      real, dimension(n), intent(out)     :: a
      integer                             :: ret

      ret = gatherlayeravg(e120(2:i1, 2:j1, :), a)
    end function get_profile_E12_

    function get_profile_T_(g_k, a, n) result(ret)
      integer, intent(in)                 :: n
      integer, dimension(n), intent(in)   :: g_k
      real, dimension(n), intent(out)     :: a
      integer                             :: ret
      
      ret = gatherlayeravg(tmp0(2:i1, 2:j1, :), a)
    end function get_profile_T_
    
    function get_zf_(g_k, a, n) result(ret)
       use modglobal, only: zf
      integer, intent(in)                 :: n
      integer, dimension(n), intent(in)   :: g_k
      real, dimension(n), intent(out)     :: a
      integer                             :: ret

      a = zf
      ret = 0
    end function get_zf_

    function get_zh_(g_k, a, n) result(ret)
      use modglobal, only: zh
      integer, intent(in)                 :: n
      integer, dimension(n), intent(in)   :: g_k
      real, dimension(n), intent(out)     :: a
      integer                             :: ret

      a = zh
      ret = 0
    end function get_zh_

    function get_presf_(g_k, a, n) result(ret)
      use modfields, only: presf
      integer, intent(in)                  :: n
      integer, dimension(n), intent(in)   :: g_k
      real, dimension(n), intent(out)     :: a
      integer                             :: ret

      a = presf
      ret = 0
    end function get_presf_

    function get_presh_(g_k, a, n) result(ret)
      use modfields, only: presh
      integer, intent(in)                  :: n
      integer, dimension(n), intent(in)   :: g_k
      real, dimension(n), intent(out)     :: a
      integer                             :: ret

      a = presh
      ret = 0
    end function get_presh_


    ! Cloud fraction getter
    ! Note: in contrast to the other profile getters,
    ! this one relies on g_k to define slabs
    ! the result a has the same length as g_k
    function get_cloudfraction(g_k, a, n) result(ret)
      integer, intent(in)                 :: n
      integer, dimension(n), intent(in)   :: g_k
      real, dimension(n), intent(out)     :: a
      integer                             :: ret
      
      ret = gatherCloudFrac(ql0(2:i1, 2:j1, :), g_k, a)
    end function get_cloudfraction    



    
!!! end of vertical profile getters

!!!!! set functions for vertical tendency vectors
    !   a is the profile to set
    !   n is the array length - assumed to be the number of layers in the model   
    function set_tendency_U(a, n) result(ret)
      integer, intent(in)                 :: n
      real, dimension(n), intent(in)      :: a
      integer                             :: ret
      u_tend = a
      ret = 0
    end function set_tendency_U
    
    function set_tendency_V(a, n) result(ret)
      integer, intent(in)                 :: n
      real, dimension(n), intent(in)      :: a
      integer                             :: ret
      v_tend = a
      ret = 0
    end function set_tendency_V

    function set_tendency_THL(a, n) result(ret)
      integer, intent(in)                 :: n
      real, dimension(n), intent(in)      :: a
      integer                             :: ret
      thl_tend = a
      ret = 0
    end function set_tendency_THL

    function set_tendency_QT(a, n) result(ret)
      integer, intent(in)                 :: n
      real, dimension(n), intent(in)      :: a
      integer                             :: ret
      qt_tend = a
      ret = 0
    end function set_tendency_QT

    function set_tendency_surface_pressure(p_tend) result(ret)
      real(8),intent(in)::  p_tend
      integer::             ret

      ps_tend = p_tend
      
      ret = 0
    end function set_tendency_surface_pressure
  
!!! end of vertical tendency setters

!!! getter functions for full 3D fields - using index arrays
    function get_field_U(g_i,g_j,g_k,a,n) result(ret)
      integer, intent(in)                 :: n
      integer, dimension(n), intent(in)   :: g_i,g_j,g_k
      real,    dimension(n), intent(out)  :: a
      integer                             :: ret
      ret = gathervol(g_i,g_j,g_k,a,n,u0(2:i1,2:j1,1:kmax))
    end function get_field_U

    function get_field_V(g_i,g_j,g_k,a,n) result(ret)
      integer, intent(in)                 :: n
      integer, dimension(n), intent(in)   :: g_i,g_j,g_k
      real,    dimension(n), intent(out)  :: a
      integer                             :: ret
      ret = gathervol(g_i,g_j,g_k,a,n,v0(2:i1,2:j1,1:kmax))
    end function get_field_V

    function get_field_W(g_i,g_j,g_k,a,n) result(ret)
      integer, intent(in)                 :: n
      integer, dimension(n), intent(in)   :: g_i,g_j,g_k
      real,    dimension(n), intent(out)  :: a
      integer                             :: ret
      ret = gathervol(g_i,g_j,g_k,a,n,w0(2:i1,2:j1,1:kmax))
    end function get_field_W

     function get_field_THL(g_i,g_j,g_k,a,n) result(ret)
       integer, intent(in)                 :: n
       integer, dimension(n), intent(in)   :: g_i,g_j,g_k
       real,    dimension(n), intent(out)  :: a
       integer                             :: ret
       ret = gathervol(g_i,g_j,g_k,a,n,thl0(2:i1,2:j1,1:kmax))
     end function get_field_THL

     function get_field_QT(g_i,g_j,g_k,a,n) result(ret)
       integer, intent(in)                 :: n
       integer, dimension(n), intent(in)   :: g_i,g_j,g_k
       real,    dimension(n), intent(out)  :: a
       integer                             :: ret
       ret = gathervol(g_i,g_j,g_k,a,n,qt0(2:i1,2:j1,1:kmax))
     end function get_field_QT

     function get_field_QL(g_i,g_j,g_k,a,n) result(ret)
       integer, intent(in)                 :: n
       integer, dimension(n), intent(in)   :: g_i,g_j,g_k
       real,    dimension(n), intent(out)  :: a
       integer                             :: ret
       ret = gathervol(g_i,g_j,g_k,a,n,ql0(2:i1,2:j1,1:kmax))
     end function get_field_QL

     function get_field_E12(g_i,g_j,g_k,a,n) result(ret)
       integer, intent(in)                 :: n
       integer, dimension(n), intent(in)   :: g_i,g_j,g_k
       real,    dimension(n), intent(out)  :: a
       integer                             :: ret
       ret = gathervol(g_i,g_j,g_k,a,n,e120(2:i1,2:j1,1:kmax))
     end function get_field_E12
     
     function get_field_T(g_i,g_j,g_k,a,n) result(ret)
       integer, intent(in)                 :: n
       integer, dimension(n), intent(in)   :: g_i,g_j,g_k
       real,    dimension(n), intent(out)  :: a
       integer                             :: ret
       ret = gathervol(g_i,g_j,g_k,a,n,tmp0(2:i1,2:j1,1:kmax))
     end function get_field_T
     !!! end of full 3D field getter functions

     ! getter function for LWP - a 2D field, vertical integral of ql
     function get_field_LWP(g_i,g_j,a,n) result(ret)
       integer, intent(in)                 :: n
       integer, dimension(n), intent(in)   :: g_i,g_j
       real,    dimension(n), intent(out)  :: a
       integer                             :: ret
       ret = gatherLWP(g_i,g_j,a,n,ql0(2:i1,2:j1,1:kmax))
     end function get_field_LWP

     ! getter function for TWP - Total Water Path - 2D field, vertical integral of qt
     function get_field_TWP(g_i,g_j,a,n) result(ret)
       integer, intent(in)                 :: n
       integer, dimension(n), intent(in)   :: g_i,g_j
       real,    dimension(n), intent(out)  :: a
       integer                             :: ret
       ret = gatherLWP(g_i,g_j,a,n,qt0(2:i1,2:j1,1:kmax))
     end function get_field_TWP

     ! Rain water path - like LWP but for rain water
     ! get it from sv0(:,:,:,iqr)
     function get_field_RWP(g_i,g_j,a,n) result(ret)
       integer, intent(in)                 :: n
       integer, dimension(n), intent(in)   :: g_i,g_j
       real,    dimension(n), intent(out)  :: a
       integer                             :: ret
       ret = gatherLWP(g_i,g_j,a,n,sv0(2:i1,2:j1,1:kmax,iqr))
     end function get_field_RWP

     
     !!! setter functions for full 3D fields - using index arrays
     !!! these functions set BOTH the -m and the -0 fields
     function set_field_U(g_i,g_j,g_k,a,n) result(ret)
       integer, intent(in)                  :: n
       integer, dimension(n), intent(in)    :: g_i,g_j,g_k
       real,    dimension(n), intent(in)    :: a
       integer                              :: ret, m, i, j, k
       !print *, 'set_field_U'      
       do m = 1,n
          if (localindex(g_i(m), g_j(m), g_k(m), i, j, k) /= 0) then
             ! print *, ' setting (', i, j, k, ') <- ', a(m)
             um(i,j,k) = a(m)
             u0(i,j,k) = a(m)
          endif
       enddo
       ret = 0
     end function set_field_U

     function set_field_V(g_i,g_j,g_k,a,n) result(ret)
       integer, intent(in)                  :: n
       integer, dimension(n), intent(in)    :: g_i,g_j,g_k
       real,    dimension(n), intent(in)    :: a
       integer                              :: ret, m, i, j, k

       do m = 1,n
          if (localindex(g_i(m), g_j(m), g_k(m), i, j, k) /= 0) then
             vm(i,j,k) = a(m)
             v0(i,j,k) = a(m)
          endif
       enddo
       ret = 0
     end function set_field_V

     function set_field_W(g_i,g_j,g_k,a,n) result(ret)
       integer, intent(in)                  :: n
       integer, dimension(n), intent(in)    :: g_i,g_j,g_k
       real,    dimension(n), intent(in)    :: a
       integer                              :: ret, m, i, j, k

       do m = 1,n
          if (localindex(g_i(m), g_j(m), g_k(m), i, j, k) /= 0) then
             wm(i,j,k) = a(m)
             w0(i,j,k) = a(m)
          endif
       enddo
       ret = 0
     end function set_field_W

     function set_field_THL(g_i,g_j,g_k,a,n) result(ret)
       integer, intent(in)                  :: n
       integer, dimension(n), intent(in)    :: g_i,g_j,g_k
       real,    dimension(n), intent(in)    :: a
       integer                              :: ret, m, i, j, k

       do m = 1,n
          if (localindex(g_i(m), g_j(m), g_k(m), i, j, k) /= 0) then
             thlm(i,j,k) = a(m)
             thl0(i,j,k) = a(m)
          endif
       enddo
       ret = 0
     end function set_field_THL
    

     function set_field_QT(g_i,g_j,g_k,a,n) result(ret)
       integer, intent(in)                  :: n
       integer, dimension(n), intent(in)    :: g_i,g_j,g_k
       real,    dimension(n), intent(in)    :: a
       integer                              :: ret, m, i, j, k

       do m = 1,n
          if (localindex(g_i(m), g_j(m), g_k(m), i, j, k) /= 0) then
             qtm(i,j,k) = a(m)
             qt0(i,j,k) = a(m)
          endif
       enddo
       ret = 0
     end function set_field_QT    
     !!! end of setter functions for 3D fields

    
    
    
!!! get simulation parameters    
    function get_params_grid (i,j,k,X,Y) result(ret)
      use modglobal, only: itot,jtot,kmax,xsize,ysize
      implicit none
      integer, intent(out)                :: i,j,k
      real, intent(out)                   :: X,Y
      integer                             :: ret
      
      i=itot
      j=jtot
      k=kmax
      X=xsize
      Y=ysize
      ret = 0
    end function get_params_grid


    


    
end module dales_interface
