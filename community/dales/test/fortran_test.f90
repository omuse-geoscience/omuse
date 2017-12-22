program fortran_test

    use dales_interface

    implicit none

    integer:: stat
    real(8):: tend
    real(8):: walltime
    stat=0

    write(*,*) "Initializing fortran interface code..."
    stat=initialize_code()
    if(stat/=0) then
        stop "...failed, stopping."
    else
        write(*,*) "..success."
    endif

    write(*,*) "Running for a minute..."
    tend=60
    stat=evolve_model(tend, 0, walltime)
    if(stat/=0) then
        stop "...failed, stopping."
    else
        write(*,*) "...success."
    endif

    write(*,*) "Exiting fortran interface code..."
    stat=cleanup_code()
    if(stat/=0) then
        stop "...failed, stopping."
    else
        write(*,*) "..success."
    endif

end program
    
