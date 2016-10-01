program fortran_test

    use dales_interface

    implicit none

    integer:: stat
    real(8):: tend

    stat=0

    write(*,*) "Initializing fortran interface code..."
    stat=initialize_code()
    if(stat/=0) then
        write(*,*) "...failed, stopping."
        return
    else
        write(*,*) "..success."
    endif

    write(*,*) "Running for a minute..."
    tend=60
    stat=evolve_model(tend)
    if(stat/=0) then
        write(*,*) "...failed, stopping."
        return
    else
        write(*,*) "...success."
    endif

    write(*,*) "Exiting fortran interface code..."
    stat=cleanup_code()
    if(stat/=0) then
        write(*,*) "...failed, stopping."
        return
    else
        write(*,*) "..success."
    endif

end program
    
