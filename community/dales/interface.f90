! This is the omuse-specific native interface to the 
! Dutch Large Eddy Simulation code DALES.
! Written by G. van den Oord, 2016

module dales_interface

    use daleslib
    use modglobal, only: rtimee,rdt

    implicit none

    character(256),         inputfile   ! Dales initializer input namelist file

contains

    function set_input_file(ifile) result(ret)
        integer,intent(out)::          ret
        character(256),intent(in)::    ifile

        inputfile=ifile
        ret=0
    end function

    function get_input_file(ifile) result(ret)
        integer,intent(out)::           ret
        character(256),intent(out)::    ifile

        ifile=inputfile
        ret=0
    end function

    function initialize_code() result(ret)
        integer,intent(out):: ret

        inputfile="namoptions"
        ret=0
    end function

    function commit_parameters() result(ret)
        integer,intent(out) :: ret

        call initialize(inputfile)
        ret=0
    end function

    function recommit_parameters() result(ret)
        integer,intent(out) :: ret

        ret=-2
    end function

    function evolve_model(tend) result(ret)
        integer,intent(out) :: ret
        real(8),intent(in)  :: tend
        integer             :: i
! TODO: Fix this when Dales runs with adaptive time-stepping
        do while(rtimee < tend - 0.5*rdt)
! TODO: Make daleslib time step signature more useful (returns time step/end
! time/etc.)
            call step
        enddo
        ret=0
    end function

    function cleanup_code() result(ret)
        integer :: ret

        call finalize
        ret=0
    end function

end module dales_interface
