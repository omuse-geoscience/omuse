! This is the omuse interface to OpenIFS. Written by G. van den Oord, 2016

module openifs_interface

    use ifslib, only: static_init, initialize, itime, finalize, get_gp_field,&
                     &get_gp_geom, get_gp_field_columns, TEMPERATURE
    use yomct0, only: nstart,nstop
    use yomct3, only: nstep
    use yomdyn, only: tstep

    implicit none

    contains

        function initialize_code() result(ret)

            integer:: ret

            call static_init()
            ret=0

        end function initialize_code

        function commit_parameters() result(ret)

            integer:: ret

            call initialize()
            ret=0

        end function

        function recommit_parameters() result(ret)

            integer:: ret

            ret=0

        end function

        function evolve_model(tend) result(ret)

            real(8), intent(in)::   tend
            integer::               ret,ii
            logical::               istat

            ret = 0
            if(istep >= nstop) then
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

        function get_field_T(g_i,g_k,a,n,k) result(ret)

            integer, intent(in)::                   n,k
            integer, dimension(n), intent(in)::     g_i,g_k
            real, dimension(n), intent(out)::       a
            integer::                               ret

            ret = 0
            call get_gp_field_columns(a,n,g_i,TEMPERATURE)

        end function get_field_T

        function cleanup_code() result(ret)

            integer:: ret

            ret=0
            call finalize()

        end function

end module openifs_interface
