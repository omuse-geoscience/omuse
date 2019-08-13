program fortran_test

    use dales_interface

    implicit none

    integer             :: stat, nargs
    real(8)             :: tend, walltime
    character(len=256)  :: case, rundir
    stat=0

    nargs = command_argument_count()
    if(nargs==0) then
        stop "Please provide a case argument (bomex/rico/sp-testcase/...)"
    end if
    call get_command_argument(1, case)
    rundir="../dales-repo/cases/"//case
    stat = change_dir(rundir)
    if(stat/=0) then
        stop "...failed, stopping."
    else
        write(*,*) "..success."
    endif

    write(*,*) "Initializing fortran interface code..."
    stat=initialize_code()
    if(stat/=0) then
        stop "...failed, stopping."
    else
        write(*,*) "..success."
    endif

    write(*,*) "Initializing model..."
    stat=commit_parameters()
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
    
