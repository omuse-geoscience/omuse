! This is the omuse interface to OpenIFS. Written by G. van den Oord, 2016

module openifs_interface

    use ifslib,     only: static_init, initialize, finalize, get_gp_field,&
                        & get_gp_geom, get_gp_field_columns, get_gp_field,& 
                        & jstep, step, PFULL, PHALF, WIND_U, WIND_V, TEMPERATURE,&
                        & SPEC_HUMIDITY, ICE_WATER, LIQ_WATER, CLOUD_FRACTION,&
                        & OZONE
    use yomct0,     only: nstart,nstop
    use yomct3,     only: nstep
    use yomdyn,     only: tstep
    use yommp0,     only: myproc
    use yomgem,     only: ngptotg
    use yomdimv,    only: nflevg

    implicit none

    contains

        function initialize_code() result(ret)

            integer:: ret

            call static_init()
            call initialize('TEST')

            ret = 0

        end function initialize_code

        function commit_parameters() result(ret)

            integer:: ret

            ret = 0

        end function

        function commit_grid() result(ret)

            integer:: ret

            ret = 0

        end function

        function recommit_parameters() result(ret)

            integer:: ret

            ret = 0

        end function

        function get_timestep(dt) result(ret)

            real(8), intent(out)::   dt
            integer::                ret

            dt = tstep
            ret = 0

        end function

        function get_grid_sizes(nxy,nz) result(ret)

            integer, intent(out)::   nxy,nz
            integer::                ret

            nxy = ngptotg
            nz = nflevg
            ret = 0

        end function

        function get_model_time(t) result(ret)

            real(8), intent(out)::  t
            integer::               ret

            t = jstep*tstep
            ret = 0

        end function

        function evolve_model(tend) result(ret)

            real(8), intent(in)::   tend
            integer::               ret,ii
            logical::               istat

            ret = 0
            if(jstep >= nstop) then
                ret = 1
                return
            endif
            istat = .true.
            do while(jstep*tstep <= tend)
                call step(istat)
                if(.not.istat) then
                    ret = 1
                    exit
                endif
            enddo

        end function evolve_model

        function evolve_model_single_step() result(ret)
          integer::               ret
          logical::               istat

          ret = 0
          call step(istat)
          if(.not.istat) then
             ret = 1
          endif
        end function evolve_model_single_step
          
        function get_gridpoints(g_i,lats,lons,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i
            real(8), dimension(n), intent(out)::    lats(n),lons(n)
            real(8), allocatable::                  b(:),c(:)
            integer::                               ret,i

            ret = 0

            if(myproc == 1) then
                allocate(b(ngptotg),c(ngptotg),stat=ret)
                if(ret/=0) then
                    return
                endif
            endif

            call get_gp_geom(b,c)
            
            if(myproc == 1) then
                do i = 1,n
                    lats(i) = b(g_i(i) + 1)
                    lons(i) = c(g_i(i) + 1)
                enddo
            endif

        end function

        function get_field_Pfull_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_field(g_i,g_k,a,n,PFULL)

        end function get_field_Pfull_

        function get_field_Phalf_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_field(g_i,g_k,a,n,PHALF)

        end function get_field_Phalf_

        function get_field_U_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_field(g_i,g_k,a,n,WIND_U)

        end function get_field_U_

        function get_field_V_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_field(g_i,g_k,a,n,WIND_V)

        end function get_field_V_

        function get_field_T_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_field(g_i,g_k,a,n,TEMPERATURE)

        end function get_field_T_

        function get_field_SH_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_field(g_i,g_k,a,n,SPEC_HUMIDITY)

        end function get_field_SH_

        function get_field_QL_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_field(g_i,g_k,a,n,LIQ_WATER)

        end function get_field_QL_

        function get_field_QI_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_field(g_i,g_k,a,n,ICE_WATER)

        end function get_field_QI_

        function get_field_O3_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_field(g_i,g_k,a,n,OZONE)

        end function get_field_O3_

        function get_field(g_i,g_k,a,n,fldid) result(ret)

            integer, intent(in)::                   n,fldid
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            real(8), allocatable::                  b(:,:)
            integer::                               ret,ii

            ret = 0
            if(myproc == 1) then
                if(fldid == PHALF) then
                    allocate(b(ngptotg,nflevg+1),stat=ret)
                else
                    allocate(b(ngptotg,nflevg),stat=ret)
                endif
                if(ret/=0) then
                    return
                endif
            endif
            
            call get_gp_field(b,fldid)
            
            if(myproc == 1) then
                do ii = 1,n
                    a(ii) = b(g_i(ii) + 1,g_k(ii) + 1)
                enddo
                deallocate(b)
            endif

        end function get_field

        function cleanup_code() result(ret)

            integer:: ret

            ret = 0
            call finalize()

        end function

end module openifs_interface
