! This is the omuse interface to OpenIFS. Written by G. van den Oord, 2016

module openifs_interface

    use ifslib,     only: static_init, initialize, finalize, get_gp_field,&
                        & get_gp_geom, get_gp_field_columns, get_gp_field,& 
                        & TEMPERATURE, jstep, step
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

        function get_gridpoints(g_i,lats,lons,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i
            real(8), dimension(n), intent(out)::    lats(n),lons(n)
            real(8), allocatable::                  b(:),c(:)
            integer::                               ret,i

            if(myproc == 1) then
                allocate(b(ngptotg),c(ngptotg))
            endif

            ret = 0
            call get_gp_geom(b,c)
            
            if(myproc == 1) then
                do i = 1,n
                    lats(i) = b(g_i(i) + 1)
                    lons(i) = c(g_i(i) + 1)
                enddo
            endif

            end function


        function get_field_T_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            real(8), allocatable::                  b(:,:)
            integer::                               ret,ii

            if(myproc == 1) then
                allocate(b(ngptotg,nflevg))
            endif
            
            ret = 0
            call get_gp_field(b,TEMPERATURE)
            
            if(myproc == 1) then
                do ii = 1,n
                    a(ii) = b(g_i(ii) + 1,g_k(ii) + 1)
                enddo
                deallocate(b)
            endif

        end function get_field_T_

        function get_profiles_T_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            real(8), allocatable::                  b(:,:)
            integer::                               ret,i,ii,k

            if(myproc == 1) then
                allocate(b(n,nflevg))
            endif
            
            ret = 0
            call get_gp_field_columns(b,n,g_i,TEMPERATURE)
            
            if(myproc == 1) then
                do i = 1,size(g_i)
                    do k = 1,nflevg
                        a(ii) = b(g_i(ii) + 1,k)
                        ii = ii + 1
                    enddo
                enddo
                deallocate(b)
            endif

        end function get_profiles_T_

        function cleanup_code() result(ret)

            integer:: ret

            ret = 0
            call finalize()

        end function

end module openifs_interface
