! This is the omuse interface to OpenIFS. Written by G. van den Oord, 2016

module openifs_interface

    use ifslib, only: static_init, initialize, istart, istop, step, finalize

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

        function evolve_model() result(ret)

            integer:: ret,ii
            logical:: istat

            ret=0
            istat=.true.
            do ii=istart,istop
                call step(istat)
                if(.not.istat) then
                    ret=1
                    exit
                endif
            enddo

        end function evolve_model

        function cleanup_code() result(ret)

            integer:: ret

            ret=0
            call finalize()

        end function

end module openifs_interface
