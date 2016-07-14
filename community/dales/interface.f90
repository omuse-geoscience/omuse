! This is the omuse-specific native interface to the 
! Dutch Large Eddy Simulation code DALES.
! Written by G. van den Oord, 2016

module dales_interface

    use daleslib

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
        integer :: ret
        
        inputfile="namoptions"
        ret=0
    end function

    function commit_parameters() result(ret)
        integer :: ret
        
        call initialize(inputfile)
        ret=0
    end function

    function recommit_parameters() result(ret)
        integer :: ret
        ret=-2
    end function

    function do_steps(n) result(ret)
        integer :: i,ret
        real(8) :: tend
        do i=1,n
            call step
        enddo
    end function

    function cleanup_code() result(ret)
        integer :: ret
        call finalize
        ret=0  
    end function

end module dales_interface
