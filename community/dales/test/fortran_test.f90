program fortran_test

    use dales_interface

    implicit none

    integer:: stat

    stat=0

    write(*,*) "Initializing fortran interface code..."
    stat=initialize_code()
    if(stat/=0) then
        write(*,*) "...failed, stopping."
        return
    else
        write(*,*) "..success."
    endif
    
    write(*,*) "Initializing dales instance..."
    stat=commit_parameters()
    if(stat/=0) then
        write(*,*) "...failed, stopping."
        return
    else
        write(*,*) "..success."
    endif

end program
    
