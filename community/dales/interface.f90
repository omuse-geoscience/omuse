! This is the omuse-specific native interface to the 
! Dutch Large Eddy Simulation code DALES.
! Written by G. van den Oord, 2016

module dales_interface

    use daleslib, only: initialize,step,finalize,&
                       allocate_z_axis,allocate_2d,allocate_3d,&
                       my_task,master_task,&
                       FIELDID_U,FIELDID_V,FIELDID_W,FIELDID_THL,FIELDID_QT, &
                       u_tend, v_tend, thl_tend, qt_tend, &
                      gatherlayeravg, gathervol
    use modfields, only: u0,v0,w0,thl0,qt0
    use modglobal, only: i1,j1,k1,itot,jtot,kmax
    
    !TODO: Expose everything so this module only depends on daleslib
    use modglobal, only: rtimee,rdt,fname_options
    use mpi, only: MPI_COMM_WORLD

    implicit none

    real(8),allocatable::tmp1Dz(:)
    real(8),allocatable::tmp2Dxy(:,:)
    real(8),allocatable::tmp3D(:,:,:)

contains

    function set_input_file(ifile) result(ret)
        integer::                      ret
        character(256),intent(in)::    ifile

        fname_options=ifile
        ret=0
    end function

    function get_input_file(ifile) result(ret)
        integer::                       ret
        character(256),intent(out)::    ifile

        ifile=fname_options
        ret=0
    end function

    function initialize_code() result(ret)
        integer:: ret

        fname_options="namoptions.001"
        call initialize(fname_options,MPI_COMM_WORLD)
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

    function evolve_model(tend) result(ret)
        integer             :: ret
        real(8),intent(in)  :: tend
        integer             :: i
        do while(rtimee < tend)
            call step
        enddo
        ret=0
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
    end function

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

      !call gatherlayeravg3(u0(2:i1, 2:j1, :), a)
      call gatherlayeravg(u0, a)
      !ret = get_field_layer_avg(FIELDID_U,a)
    end function get_profile_U_

    function get_profile_V_(g_k, a, n) result(ret)
      integer, intent(in)                 :: n
      integer, dimension(n), intent(in)   :: g_k
      real, dimension(n), intent(out)     :: a
      integer                             :: ret

      !call gatherlayeravg3(v0(2:i1, 2:j1, :), a)
      call gatherlayeravg(v0, a)
      !ret = get_field_layer_avg(FIELDID_V,a)
    end function get_profile_V_
 
    function get_profile_W_(g_k, a, n) result(ret)
      integer, intent(in)                 :: n
      integer, dimension(n), intent(in)   :: g_k
      real, dimension(n), intent(out)     :: a
      integer                             :: ret

      !call gatherlayeravg3(w0(2:i1, 2:j1, :), a)
      call gatherlayeravg(w0, a)
      !ret = get_field_layer_avg(FIELDID_W,a)
    end function get_profile_W_

    function get_profile_THL_(g_k, a, n) result(ret)
      integer, intent(in)                 :: n
      integer, dimension(n), intent(in)   :: g_k
      real, dimension(n), intent(out)     :: a
      integer                             :: ret

      !call gatherlayeravg3(thl0(2:i1, 2:j1, :), a)
      call gatherlayeravg(thl0, a)
      !ret = get_field_layer_avg(FIELDID_THL,a)
    end function get_profile_THL_

    function get_profile_QT_(g_k, a, n) result(ret)
      integer, intent(in)                 :: n
      integer, dimension(n), intent(in)   :: g_k
      real, dimension(n), intent(out)     :: a
      integer                             :: ret

      !call gatherlayeravg3(qt0(2:i1, 2:j1, :), a)
      call gatherlayeravg(qt0, a)
      !ret = get_field_layer_avg(FIELDID_QT,a)
    end function get_profile_QT_

     function get_zf_(g_k, a, n) result(ret)
       use modglobal, only: zf
      integer, intent(in)                 :: n
      integer, dimension(n), intent(in)   :: g_k
      real, dimension(n), intent(out)     :: a
      integer                             :: ret
      integer                             :: i

      do i = 1, n
         a(i) = zf(i)
      end do 
      ret = 0
    end function get_zf_

    function get_zh_(g_k, a, n) result(ret)
      use modglobal, only: zh
      integer, intent(in)                 :: n
      integer, dimension(n), intent(in)   :: g_k
      real, dimension(n), intent(out)     :: a
      integer                             :: ret
      integer                             :: i
      
      do i = 1, n
         a(i) = zh(i)
      end do
      ret = 0
    end function get_zh_

!!! end of vertical profile getters

!!!!! set functions for vertical tendency vectors
    !   a is the profile to set
    !   n is the array length - assumed to be the number of layers in the model   
    function set_tendency_U(a, n) result(ret)
      integer, intent(in)                 :: n
      real, dimension(n), intent(in)      :: a
      integer                             :: ret
      u_tend = a
    end function set_tendency_U
    
    function set_tendency_V(a, n) result(ret)
      integer, intent(in)                 :: n
      real, dimension(n), intent(in)      :: a
      integer                             :: ret
      v_tend = a
    end function set_tendency_V

    function set_tendency_THL(a, n) result(ret)
      integer, intent(in)                 :: n
      real, dimension(n), intent(in)      :: a
      integer                             :: ret
      thl_tend = a
    end function set_tendency_THL

    function set_tendency_QT(a, n) result(ret)
      integer, intent(in)                 :: n
      real, dimension(n), intent(in)      :: a
      integer                             :: ret
      qt_tend = a
    end function set_tendency_QT
!!! end of vertical tendency setters

!!! getter functions for full 3D fields - using index arrays
    function get_field_U(g_i,g_j,g_k,a,n) result(ret)
      integer, intent(in)                 :: n
       integer, dimension(n), intent(in)   :: g_i,g_j,g_k
       real,    dimension(n), intent(out)  :: a
       integer                             :: ret
       ret = gathervol(g_i,g_j,g_k,a,n,u0)
    end function get_field_U

    function get_field_V(g_i,g_j,g_k,a,n) result(ret)
      integer, intent(in)                 :: n
       integer, dimension(n), intent(in)   :: g_i,g_j,g_k
       real,    dimension(n), intent(out)  :: a
       integer                             :: ret
       ret = gathervol(g_i,g_j,g_k,a,n,v0)
    end function get_field_V

    function get_field_W(g_i,g_j,g_k,a,n) result(ret)
      integer, intent(in)                 :: n
       integer, dimension(n), intent(in)   :: g_i,g_j,g_k
       real,    dimension(n), intent(out)  :: a
       integer                             :: ret
       ret = gathervol(g_i,g_j,g_k,a,n,w0)
     end function get_field_W

     function get_field_THL(g_i,g_j,g_k,a,n) result(ret)
      integer, intent(in)                 :: n
       integer, dimension(n), intent(in)   :: g_i,g_j,g_k
       real,    dimension(n), intent(out)  :: a
       integer                             :: ret
       ret = gathervol(g_i,g_j,g_k,a,n,thl0)
     end function get_field_THL

     function get_field_QT(g_i,g_j,g_k,a,n) result(ret)
      integer, intent(in)                 :: n
       integer, dimension(n), intent(in)   :: g_i,g_j,g_k
       real,    dimension(n), intent(out)  :: a
       integer                             :: ret
       ret = gathervol(g_i,g_j,g_k,a,n,qt0)
     end function get_field_QT
!!! end of full 3D field getter functions
    

    
    
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
