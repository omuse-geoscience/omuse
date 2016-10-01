! This is the omuse-specific native interface to the 
! Dutch Large Eddy Simulation code DALES.
! Written by G. van den Oord, 2016

module dales_interface

    use daleslib, only: initialize,step,finalize,get_field_layer_avg,&
                       &allocate_z_axis,allocate_2d,allocate_3d,&
                       &get_field_2d,get_field_3d,my_task,master_task
    !TODO: Expose everything so this module only depends on daleslib
    use modglobal, only: rtimee,rdt,fname_options
    use mpi, only: MPI_COMM_WORLD

    implicit none

    real(8),allocatable::tmp1Dz(:)
    real(8),allocatable::tmp2Dxy(:,:)
    real(8),allocatable::tmp3D(:,:,:)

contains

    function set_input_file(ifile) result(ret)
        implicit none

        integer::                      ret
        character(256),intent(in)::    ifile

        fname_options=ifile
        ret=0
    end function

    function get_input_file(ifile) result(ret)
        implicit none

        integer::                       ret
        character(256),intent(out)::    ifile

        ifile=fname_options
        ret=0
    end function

    function initialize_code() result(ret)
        implicit none

        integer:: ret

        fname_options="namoptions.001"
        call initialize(fname_options,MPI_COMM_WORLD)
        ret=0
    end function

    function commit_parameters() result(ret)
        implicit none

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
        implicit none

        integer:: ret
        ret=0
    end function

    function evolve_model(tend) result(ret)
        implicit none

        integer             :: ret
        real(8),intent(in)  :: tend
        integer             :: i
        do while(rtimee < tend)
            call step
        enddo
        ret=0
    end function

    function cleanup_code() result(ret)
        implicit none

        integer :: ret

        call finalize
        ret=0
    end function

    function get_timestep(dt) result(ret)
        implicit none

        real(8),intent(out):: dt
        integer            :: ret

        dt=rdt
        ret=0
    end function

    function get_model_time(t) result(ret)
        implicit none

        real(8),intent(out):: t
        integer::             ret

        t=rtimee
        ret=0
    end function

    function commit_grid() result(ret)
        implicit none

        integer::             ret

        ret=0
    end function

    function get_profile_field(g_k,a,n) result(ret)
        implicit none

        integer, intent(in)                 :: n
        integer, dimension(n), intent(in)   :: g_k
        real, dimension(n), intent(out)     :: a
        integer                             :: ret,j

        ret=get_field_layer_avg(5,tmp1Dz)
        if(ret/=0) then
            return
        endif
        if(my_task==master_task) then
            do j=1,n
                a(j)=tmp1Dz(g_k(j))
            enddo
        endif
    end function

    function get_layer_field(g_i,g_j,k,a,n) result(ret)
        implicit none
    
        integer, intent(in)                 :: n,k
        integer, dimension(n), intent(in)   :: g_i,g_j
        real,    dimension(n), intent(out)  :: a
        integer                             :: ret,j

        ret=get_field_2d(1,k,tmp2Dxy)
        if(ret/=0) then
            return
        endif
        if(my_task==master_task) then
            do j=1,n
                a(j)=tmp2Dxy(g_i(j),g_j(j))
            enddo
        endif
    end function

    function get_volume_field(g_i,g_j,g_k,a,n) result(ret)
        implicit none
    
        integer, intent(in)                 :: n
        integer, dimension(n), intent(in)   :: g_i,g_j,g_k
        real,    dimension(n), intent(out)  :: a
        integer                             :: ret,j

        ret=get_field_3d(1,tmp3D)
        if(ret/=0) then
            return
        endif
        if(my_task==master_task) then
            do j=1,n
                a(j)=tmp3D(g_i(j),g_j(j),g_k(j))
            enddo
        endif
    end function

end module dales_interface
