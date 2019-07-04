! This is the omuse-specific native interface to the
! Dutch Large Eddy Simulation code DALES.
! Written by G. van den Oord, 2016

module dales_interface

    use daleslib, only : initialize, step, finalize, &
            allocate_z_axis, allocate_2d, allocate_3d, &
            my_task, master_task, gatherlayer, &
            FIELDID_U, FIELDID_V, FIELDID_W, FIELDID_THL, FIELDID_QT, &
            u_tend, v_tend, thl_tend, qt_tend, ql_tend, ps_tend, ql_ref, qt_alpha, &
            gatherlayeravg, gathervol, localindex, gatherLWP, gatherCloudFrac, gather_ice, &
            qt_forcing_type, QT_FORCING_GLOBAL, QT_FORCING_LOCAL, QT_FORCING_VARIANCE, &
            u_nudge, v_nudge, thl_nudge, qt_nudge, u_nudge_time, v_nudge_time, thl_nudge_time, qt_nudge_time

    !TODO: Expose everything so this module only depends on daleslib
    use modfields, only : u0, v0, w0, thl0, qt0, ql0, e120, tmp0, sv0, um, vm, wm, thlm, qtm !, surf_rain
    use modglobal, only : i1, j1, k1, itot, jtot, kmax, lmoist
    use modsurfdata, only : ps, ustar, z0m, z0h, tskin, qskin, LE, H, obl, thlflux, qtflux, svflux, dudz, dvdz, &
            dqtdz, dthldz, thls, qts, thvs, svs, ustin, wqsurf, wtsurf, wsvsurf, z0, z0mav, z0hav
    use modraddata, only : swd, swdir, swdif, swu, lwd, lwu, swdca, swuca, lwdca, lwuca
    use modmicrodata, only : iqr
    use modglobal, only : rtimee, rdt, fname_options, timeleft, tres, timee, rk3step, dt_lim
    use modstartup, only : do_writerestartfiles
    use mpi, only : MPI_COMM_WORLD

    implicit none

    real(8), allocatable :: tmp1Dz(:)
    real(8), allocatable :: tmp2Dxy(:, :)
    real(8), allocatable :: tmp3D(:, :, :)

    integer :: MPI_local_comm = MPI_COMM_WORLD
    integer :: startdate = 0
    integer :: starttime = 0

    logical :: exactEnd = .false.

contains

    ! Sets the variance nudging type
    function set_qt_forcing(forcing_type) result(ret)
        integer :: ret
        integer, intent(in) :: forcing_type
        if (forcing_type >= 0 .and. forcing_type <= 2) then
            qt_forcing_type = forcing_type
            ret = 0
        else
            ret = 1
        endif
    end function set_qt_forcing

    ! Sets the dales input file
    function set_input_file(ifile) result(ret)
        integer :: ret
        character(256), intent(in) :: ifile

        fname_options = ifile
        ret = 0
    end function set_input_file

    ! Returns the dales input file
    function get_input_file(ifile) result(ret)
        integer :: ret
        character(256), intent(out) :: ifile

        ifile = fname_options
        ret = 0
    end function get_input_file

    ! set the working directory of the dales worker process
    function change_dir(directory) result(ret)
#if defined (__INTEL_COMPILER)
        USE IFPORT   ! for intel chdir function.
#endif
        integer :: ret
        character(256), intent(in) :: directory
        ret = chdir(directory)
    end function change_dir

    ! Sets the simulation starting date
    function set_start_date(date) result(ret)
        integer :: ret
        integer, intent(in) :: date

        startdate = date
        ret = 0
    end function set_start_date

    ! Sets the simulation starting timestamp
    function set_start_time(time) result(ret)
        integer :: ret
        integer, intent(in) :: time

        starttime = time
        ret = 0
    end function set_start_time

    ! Dales instantiation function
    function initialize_code() result(ret)
        integer :: ret

        fname_options = "namoptions.001"
        !startdate=0
        !starttime=0
        ret = 0
    end function

    ! Dales initializer function
    function commit_parameters() result(ret)
        integer :: ret

        call initialize(fname_options, mpi_comm = MPI_local_comm, date = startdate, time = starttime)

        ret = 0
        if(my_task==master_task) then
            ret = allocate_z_axis(tmp1Dz)
            if(ret/=0) return
            ret = allocate_2d(tmp2dxy)
            if(ret/=0) return
            ret = allocate_3d(tmp3D)
        endif
    end function

    ! Trivial placeholder method
    function recommit_parameters() result(ret)
        integer :: ret
        ret = 0
    end function

    ! Switch from/to adaptive time stepping
    function set_exact_end(t) result(ret)
        integer :: ret
        logical, intent(in) :: t

        exactEnd = t
        ret = 0
    end function set_exact_end

    ! Returns the adaptive time stepping flag
    function get_exact_end(t) result(ret)
        integer :: ret
        logical, intent(out) :: t

        t = exactEnd
        ret = 0
    end function get_exact_end

    ! evolve the model to time tend
    ! returns 1 if the end of the simulation was reached (timeleft == 0)
    !
    ! EXPERIMENTAL: if exactEnd is nonzero,
    !               set dt_lim so that the timestepping finishes exactly at tend
    !               then the end time from namoptions is not used.
    !
    !  .or. rk3step < 3) SHOULD BE USED TO FINISH A FULL TIME STEP BEFORE EXITING ??
    !                    but it breaks the current result test, which was stored without this.
    !                    see program.f90: main loop - while(timeleft>0 .or. rk3step < 3)
    !
    !  Returns the elapsed time in seconds (real) in walltime
    function evolve_model(tend, exactEnd_, walltime) result(ret)
        use mpi, only : MPI_Wtime
        integer :: ret
        real(8), intent(in) :: tend
        integer, intent(in) :: exactEnd_
        real(8), intent(out) :: walltime
        integer :: i
        real(8) :: start_time

        start_time = MPI_Wtime()

        do while((rtimee < tend .and. timeleft > 0)) ! .or. rk3step < 3)

            if (exactEnd_ /= 0 .or. exactEnd) then
                ! set dt_lim to finish at tend (or before)
                dt_lim = min(dt_lim, (int(tend / tres) - timee))
            endif

            ! thrice to ensure consistency of state (rkstep==3)
            call step
            call step
            call step
            !           print *, '  ', 'rtimee=',rtimee, 'timeleft=',timeleft, 'rk3step=',rk3step, 'rdt=',rdt
        enddo

        walltime = MPI_Wtime() - start_time
        ret = 0
        if(timeleft <= 0) then
            ret = 1
        endif

    end function

    ! Dales finalizer method
    function cleanup_code() result(ret)
        integer :: ret

        call finalize
        ret = 0
    end function

    ! Returns the model fixed timestep
    function get_timestep(dt) result(ret)
        real(8), intent(out) :: dt
        integer :: ret

        dt = rdt
        ret = 0
    end function

    ! Returns the model time
    function get_model_time(t) result(ret)
        real(8), intent(out) :: t
        integer :: ret

        t = rtimee
        ret = 0
    end function get_model_time

    ! Sets the (uniform) surface pressure
    function set_surface_pressure(p) result(ret)

        use modsurface, only : qtsurf

        real(8), intent(in) :: p
        integer :: ret

        ps = p

        ! modtimedep calls qtsurf to update surface values after changing ps, so we'll do the same
        if (lmoist) then
            call qtsurf
        else
            qts = 0.
        endif

        ret = 0
    end function set_surface_pressure

    ! Returns the model surface pressure
    function get_surface_pressure(p) result(ret)
        real(8), intent(out) :: p
        integer :: ret

        p = ps
        ret = 0
    end function get_surface_pressure

    ! Trivial placeholder method
    function commit_grid() result(ret)
        integer :: ret
        ret = 0
    end function

    !!!!! get functions for vertical profiles / slab averages
    !   g_k is a dummy array
    !   a is the returned result
    !   n is the array length

    ! Returns eastward wind profile
    function get_profile_U_(g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_k
        real, dimension(n), intent(out) :: a
        integer :: ret

        ret = gatherlayeravg(u0(2:i1, 2:j1, g_k(1:n)), a)
    end function get_profile_U_

    ! Returns northward wind profile
    function get_profile_V_(g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_k
        real, dimension(n), intent(out) :: a
        integer :: ret

        ret = gatherlayeravg(v0(2:i1, 2:j1, g_k(1:n)), a)
    end function get_profile_V_

    ! Returns vertical wind profile
    function get_profile_W_(g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_k
        real, dimension(n), intent(out) :: a
        integer :: ret

        ret = gatherlayeravg(w0(2:i1, 2:j1, g_k(1:n)), a)
    end function get_profile_W_

    ! Returns the vertical liquid pot.temp. profile
    function get_profile_THL_(g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_k
        real, dimension(n), intent(out) :: a
        integer :: ret

        ret = gatherlayeravg(thl0(2:i1, 2:j1, g_k(1:n)), a)
    end function get_profile_THL_

    ! Returns the total humidity profile
    function get_profile_QT_(g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_k
        real, dimension(n), intent(out) :: a
        integer :: ret

        ret = gatherlayeravg(qt0(2:i1, 2:j1, g_k(1:n)), a)
    end function get_profile_QT_

    ! Returns the liquid water content profile
    function get_profile_QL_(g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_k
        real, dimension(n), intent(out) :: a
        integer :: ret

        ret = gatherlayeravg(ql0(2:i1, 2:j1, g_k(1:n)), a)
    end function get_profile_QL_

    ! Returns the ice water content profile
    function get_profile_QL_ice_(g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_k
        real, dimension(n), intent(out) :: a
        real, dimension(kmax) :: a_
        integer :: i, ret

        ret = gather_ice(a_)
        do i = 1, n
            a(i) = a_(g_k(i))
        enddo
    end function get_profile_QL_ice_

    ! Returns the rain water profile
    function get_profile_QR_(g_k, a, n) result(ret)
        use modmicrodata, only : iqr
        use modglobal, only : nsv
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_k
        real, dimension(n), intent(out) :: a
        integer :: ret

        if (nsv >= iqr) then
            ret = gatherlayeravg(sv0(2:i1, 2:j1, g_k(1:n), iqr), a)
        else
            a = 0  ! Dales doesn't have the requested rain field
            ret = 1
        endif

    end function get_profile_QR_

    ! Returns the TKE profile
    function get_profile_E12_(g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_k
        real, dimension(n), intent(out) :: a
        integer :: ret

        ret = gatherlayeravg(e120(2:i1, 2:j1, g_k(1:n)), a)
    end function get_profile_E12_

    ! Returns the real temperature profile
    function get_profile_T_(g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_k
        real, dimension(n), intent(out) :: a
        integer :: ret

        ret = gatherlayeravg(tmp0(2:i1, 2:j1, g_k(1:n)), a)
    end function get_profile_T_

    ! Returns the vertical layer center heights
    function get_zf_(g_k, a, n) result(ret)
        use modglobal, only : zf
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_k
        real, dimension(n), intent(out) :: a
        integer :: ret

        a = zf(g_k(1:n))
        ret = 0
    end function get_zf_

    ! Returns the vertical layer interface heights
    function get_zh_(g_k, a, n) result(ret)
        use modglobal, only : zh
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_k
        real, dimension(n), intent(out) :: a
        integer :: ret

        a = zh(g_k(1:n))
        ret = 0
    end function get_zh_

    ! Returns the vertical layer center pressures
    function get_presf_(g_k, a, n) result(ret)
        use modfields, only : presf
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_k
        real, dimension(n), intent(out) :: a
        integer :: ret

        a = presf(g_k(1:n))
        ret = 0
    end function get_presf_

    ! Returns the vertical layer interface pressures
    function get_presh_(g_k, a, n) result(ret)
        use modfields, only : presh
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_k
        real, dimension(n), intent(out) :: a
        integer :: ret

        a = presh(g_k(1:n))
        ret = 0
    end function get_presh_

    ! Returns the full level density
    function get_rhof_(g_k, a, n) result(ret)
        use modfields, only : rhof
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_k
        real, dimension(n), intent(out) :: a
        integer :: ret

        a = rhof(g_k(1:n))
        ret = 0
    end function get_rhof_

      ! Returns the full level base density
    function get_rhobf_(g_k, a, n) result(ret)
        use modfields, only : rhobf
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_k
        real, dimension(n), intent(out) :: a
        integer :: ret

        a = rhobf(g_k(1:n))
        ret = 0
    end function get_rhobf_
      

    ! Cloud fraction getter
    ! Note: in contrast to the other profile getters,
    ! this one relies on g_k to define slabs
    ! the result a has the same length as g_k
    function get_cloudfraction(g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_k
        real, dimension(n), intent(out) :: a
        integer :: ret

        ret = gatherCloudFrac(ql0(2:i1, 2:j1, :), g_k, a)
    end function get_cloudfraction
    !!! end of vertical profile getters

    ! Total rain water content getter
!    function get_rain(rain) result(ret)
!        use mpi
!        use modmpi, only : comm3d, my_real, mpierr, myid, nprocs
!        use modglobal, only : itot, jtot
!
!        real, intent(out) :: rain
!        integer :: ret
!
!        rain = sum(surf_rain(2:i1, 2:j1))
!        write (*,*) 'local rain sum', rain
!
!        !in-place reduction
!        if (myid == 0) then
!            CALL mpi_reduce(MPI_IN_PLACE, rain, 1, MY_REAL, MPI_SUM, 0, comm3d, ret)
!            rain = rain / (itot * jtot)
!            write (*,*) 'global rain avg', rain
!        else
!            CALL mpi_reduce(rain, rain, 1, MY_REAL, MPI_SUM, 0, comm3d, ret)
!        endif
!        ret = 0
!    end function get_rain

    ! set functions for vertical tendency vectors
    !   a is the profile to set
    !   n is the array length - assumed to be the number of layers in the model
    ! Sets the eastward wind tendency profile
    function set_tendency_U(a, n) result(ret)
        integer, intent(in) :: n
        real, dimension(n), intent(in) :: a
        integer :: ret
        u_tend = a
        ret = 0
    end function set_tendency_U

    ! Sets the northward wind tendency profile
    function set_tendency_V(a, n) result(ret)
        integer, intent(in) :: n
        real, dimension(n), intent(in) :: a
        integer :: ret
        v_tend = a
        ret = 0
    end function set_tendency_V

    ! Sets the liquid water pot. temperature profile
    function set_tendency_THL(a, n) result(ret)
        integer, intent(in) :: n
        real, dimension(n), intent(in) :: a
        integer :: ret
        thl_tend = a
        ret = 0
    end function set_tendency_THL

    ! Sets the total humidity tendency profile
    function set_tendency_QT(a, n) result(ret)
        integer, intent(in) :: n
        real, dimension(n), intent(in) :: a
        integer :: ret
        qt_tend = a
        ret = 0
    end function set_tendency_QT

    ! Sets the liquid water content tendency
    function set_tendency_QL(a, n) result(ret)
        integer, intent(in) :: n
        real, dimension(n), intent(in) :: a
        integer :: ret
        ql_tend = a
        ret = 0
    end function set_tendency_QL

    ! Sets the liquid water reference profile
    function set_ref_profile_QL(a, n) result(ret)
        integer, intent(in) :: n
        real, dimension(n), intent(in) :: a
        integer :: ret
        ql_ref = a
        ret = 0
    end function set_ref_profile_QL

    ! Sets the variance nudging factor
    function set_qt_variability_factor(a, n) result(ret)
        integer, intent(in) :: n
        real, dimension(n), intent(in) :: a
        integer :: ret
        qt_alpha = a
        ret = 0
    end function set_qt_variability_factor

    !!!!! indexed get/set functions for vertical tendency vectors
    ! Sets the eastward velocity tendency profile
    function set_tendency_U_(g_k, a, n) result(ret)
        integer, intent(in) :: n
        real, dimension(n), intent(in) :: a
        integer, dimension(n), intent(in) :: g_k
        integer :: ret
        u_tend(g_k(1:n)) = a(1:n)
        ret = 0
    end function set_tendency_U_

    ! Returns the eastward velocity tendency profile
    function get_tendency_U_(g_k, a, n) result(ret)
        integer, intent(in) :: n
        real, dimension(n), intent(out) :: a
        integer, dimension(n), intent(in) :: g_k
        integer :: ret
        a(1:n) = u_tend(g_k(1:n))
        ret = 0
    end function get_tendency_U_

    ! Sets the northward velocity tendency profile
    function set_tendency_V_(g_k, a, n) result(ret)
        integer, intent(in) :: n
        real, dimension(n), intent(in) :: a
        integer, dimension(n), intent(in) :: g_k
        integer :: ret
        v_tend(g_k(1:n)) = a(1:n)
        ret = 0
    end function set_tendency_V_

    ! Returns the northward velocity tendency profile
    function get_tendency_V_(g_k, a, n) result(ret)
        integer, intent(in) :: n
        real, dimension(n), intent(out) :: a
        integer, dimension(n), intent(in) :: g_k
        integer :: ret
        a(1:n) = V_tend(g_k(1:n))
        ret = 0
    end function get_tendency_V_

    ! Sets the liquid water pot. temp. tendency profile
    function set_tendency_THL_(g_k, a, n) result(ret)
        integer, intent(in) :: n
        real, dimension(n), intent(in) :: a
        integer, dimension(n), intent(in) :: g_k
        integer :: ret
        thl_tend(g_k(1:n)) = a(1:n)
        ret = 0
    end function set_tendency_THL_

    ! Returns the liquid water pot. temp. tendency profile
    function get_tendency_THL_(g_k, a, n) result(ret)
        integer, intent(in) :: n
        real, dimension(n), intent(out) :: a
        integer, dimension(n), intent(in) :: g_k
        integer :: ret
        a(1:n) = thl_tend(g_k(1:n))
        ret = 0
    end function get_tendency_THL_

    ! Sets the total humidity tendency profile
    function set_tendency_QT_(g_k, a, n) result(ret)
        integer, intent(in) :: n
        real, dimension(n), intent(in) :: a
        integer, dimension(n), intent(in) :: g_k
        integer :: ret
        qt_tend(g_k(1:n)) = a(1:n)
        ret = 0
    end function set_tendency_QT_

    ! Returns the total humidity tendency profile
    function get_tendency_QT_(g_k, a, n) result(ret)
        integer, intent(in) :: n
        real, dimension(n), intent(out) :: a
        integer, dimension(n), intent(in) :: g_k
        integer :: ret
        a(1:n) = qt_tend(g_k(1:n))
        ret = 0
    end function get_tendency_QT_

    ! Sets the surface pressure tendency
    function set_tendency_surface_pressure(p_tend) result(ret)
        real(8), intent(in) :: p_tend
        integer :: ret

        ps_tend = p_tend

        ret = 0
    end function set_tendency_surface_pressure

    ! Returns the surface pressure tendency
    function get_tendency_surface_pressure(p_tend) result(ret)
        real(8), intent(out) :: p_tend
        integer :: ret

        p_tend = ps_tend

        ret = 0
    end function get_tendency_surface_pressure
    !!! end of vertical tendency setters

    !!!!! indexed get/set functions for nudge targets
    ! Sets the eastward wind velocity nudging relaxation time
    function set_nudge_time_U(time) result(ret)
        real, intent(in) :: time
        integer :: ret
        u_nudge_time = time
        ret = 0
    end function set_nudge_time_U

    ! Returns the eastward wind velocity nudging relaxation time
    function get_nudge_time_U(time) result(ret)
        real, intent(out) :: time
        integer :: ret
        time = u_nudge_time
        ret = 0
    end function get_nudge_time_U

    ! Sets the northward wind velocity nudging relaxation time
    function set_nudge_time_V(time) result(ret)
        real, intent(in) :: time
        integer :: ret
        v_nudge_time = time
        ret = 0
    end function set_nudge_time_V

    ! Returns the northward wind velocity nudging relaxation time
    function get_nudge_time_V(time) result(ret)
        real, intent(out) :: time
        integer :: ret
        time = v_nudge_time
        ret = 0
    end function get_nudge_time_V

    ! Sets the liquid water pot. temp. nudging relaxation time
    function set_nudge_time_THL(time) result(ret)
        real, intent(in) :: time
        integer :: ret
        thl_nudge_time = time
        ret = 0
    end function set_nudge_time_THL

    ! Returns the liquid water pot. temp. nudging relaxation time
    function get_nudge_time_THL(time) result(ret)
        real, intent(out) :: time
        integer :: ret
        time = thl_nudge_time
        ret = 0
    end function get_nudge_time_THL

    ! Sets the total humidity nudging relaxation time
    function set_nudge_time_QT(time) result(ret)
        real, intent(in) :: time
        integer :: ret
        qt_nudge_time = time
        ret = 0
    end function set_nudge_time_QT

    ! Returns the total humidity nudging relaxation time
    function get_nudge_time_QT(time) result(ret)
        real, intent(out) :: time
        integer :: ret
        time = qt_nudge_time
        ret = 0
    end function get_nudge_time_QT

    ! Sets the eastward velocity nudging profile
    function set_nudge_U(g_k, a, n) result(ret)
        integer, intent(in) :: n
        real, dimension(n), intent(in) :: a
        integer, dimension(n), intent(in) :: g_k
        integer :: ret
        u_nudge(g_k(1:n)) = a(1:n)
        ret = 0
    end function set_nudge_U

    ! Returns the eastward velocity nudging profile
    function get_nudge_U(g_k, a, n) result(ret)
        integer, intent(in) :: n
        real, dimension(n), intent(out) :: a
        integer, dimension(n), intent(in) :: g_k
        integer :: ret
        a(1:n) = u_nudge(g_k(1:n))
        ret = 0
    end function get_nudge_U

    ! Sets the northward velocity nudging profile
    function set_nudge_V(g_k, a, n) result(ret)
        integer, intent(in) :: n
        real, dimension(n), intent(in) :: a
        integer, dimension(n), intent(in) :: g_k
        integer :: ret
        v_nudge(g_k(1:n)) = a(1:n)
        ret = 0
    end function set_nudge_V

    ! Returns the northward velocity nudging profile
    function get_nudge_V(g_k, a, n) result(ret)
        integer, intent(in) :: n
        real, dimension(n), intent(out) :: a
        integer, dimension(n), intent(in) :: g_k
        integer :: ret
        a(1:n) = V_nudge(g_k(1:n))
        ret = 0
    end function get_nudge_V

    ! Sets the liquid water pot. temperature nudging profile
    function set_nudge_THL(g_k, a, n) result(ret)
        integer, intent(in) :: n
        real, dimension(n), intent(in) :: a
        integer, dimension(n), intent(in) :: g_k
        integer :: ret
        thl_nudge(g_k(1:n)) = a(1:n)
        ret = 0
    end function set_nudge_THL

    ! Returns the liquid water pot. temperature nudging profile
    function get_nudge_THL(g_k, a, n) result(ret)
        integer, intent(in) :: n
        real, dimension(n), intent(out) :: a
        integer, dimension(n), intent(in) :: g_k
        integer :: ret
        a(1:n) = thl_nudge(g_k(1:n))
        ret = 0
    end function get_nudge_THL

    ! Sets the total humidity nudging profile
    function set_nudge_QT(g_k, a, n) result(ret)
        integer, intent(in) :: n
        real, dimension(n), intent(in) :: a
        integer, dimension(n), intent(in) :: g_k
        integer :: ret
        qt_nudge(g_k(1:n)) = a(1:n)
        ret = 0
    end function set_nudge_QT

    ! Returns the total humidity nudging profile
    function get_nudge_QT(g_k, a, n) result(ret)
        integer, intent(in) :: n
        real, dimension(n), intent(out) :: a
        integer, dimension(n), intent(in) :: g_k
        integer :: ret
        a(1:n) = qt_nudge(g_k(1:n))
        ret = 0
    end function get_nudge_QT
    !!! end of nudge target getters/setters

    !!! getter functions for full 3D fields - using index arrays
    ! Eastward wind velocity volume field getter
    function get_field_U(g_i, g_j, g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j, g_k
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gathervol(g_i, g_j, g_k, a, n, u0(2:i1, 2:j1, 1:kmax))
    end function get_field_U

    ! Northward wind velocity volume field getter
    function get_field_V(g_i, g_j, g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j, g_k
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gathervol(g_i, g_j, g_k, a, n, v0(2:i1, 2:j1, 1:kmax))
    end function get_field_V

    ! Upward wind velocity volume field getter
    function get_field_W(g_i, g_j, g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j, g_k
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gathervol(g_i, g_j, g_k, a, n, w0(2:i1, 2:j1, 1:kmax))
    end function get_field_W

    ! Liquid water potential temperature volume field getter
    function get_field_THL(g_i, g_j, g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j, g_k
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gathervol(g_i, g_j, g_k, a, n, thl0(2:i1, 2:j1, 1:kmax))
    end function get_field_THL

    ! Total humidity volume field getter
    function get_field_QT(g_i, g_j, g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j, g_k
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gathervol(g_i, g_j, g_k, a, n, qt0(2:i1, 2:j1, 1:kmax))
    end function get_field_QT

    ! Liquid water content volume field getter
    function get_field_QL(g_i, g_j, g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j, g_k
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gathervol(g_i, g_j, g_k, a, n, ql0(2:i1, 2:j1, 1:kmax))
    end function get_field_QL

    ! Saturation humidity volume field getter
    !function get_field_Qsat(g_i, g_j, g_k, a, n) result(ret)
    !    integer, intent(in) :: n
    !    integer, dimension(n), intent(in) :: g_i, g_j, g_k
    !    real, dimension(n), intent(out) :: a
    !    integer :: ret
    !    ret = gathervol(g_i, g_j, g_k, a, n, qsat(2:i1, 2:j1, 1:kmax))
    !end function get_field_Qsat

    ! Turbulent kinetic energy volume field getter
    function get_field_E12(g_i, g_j, g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j, g_k
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gathervol(g_i, g_j, g_k, a, n, e120(2:i1, 2:j1, 1:kmax))
    end function get_field_E12

    ! Temperature volume field getter
    function get_field_T(g_i, g_j, g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j, g_k
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gathervol(g_i, g_j, g_k, a, n, tmp0(2:i1, 2:j1, 1:kmax))
    end function get_field_T

    ! Downwelling shortwave radiative flux getter
    function get_field_rswd(g_i, g_j, g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j, g_k
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gathervol(g_i, g_j, g_k, a, n, swd(2:i1, 2:j1, 1:kmax))
    end function get_field_rswd

    ! Upwelling shortwave radiative flux getter
    function get_field_rswu(g_i, g_j, g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j, g_k
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gathervol(g_i, g_j, g_k, a, n, swu(2:i1, 2:j1, 1:kmax))
    end function get_field_rswu

    ! Downwelling direct shortwave radiative flux getter
    function get_field_rswdir(g_i, g_j, g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j, g_k
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gathervol(g_i, g_j, g_k, a, n, swdir(2:i1, 2:j1, 1:kmax))
    end function get_field_rswdir

    ! Downwelling diffuse shortwave radiative flux getter
    function get_field_rswdif(g_i, g_j, g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j, g_k
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gathervol(g_i, g_j, g_k, a, n, swdif(2:i1, 2:j1, 1:kmax))
    end function get_field_rswdif

    ! Downwelling longwave radiative flux getter
    function get_field_rlwd(g_i, g_j, g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j, g_k
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gathervol(g_i, g_j, g_k, a, n, lwd(2:i1, 2:j1, 1:kmax))
    end function get_field_rlwd

    ! Upwelling longwave radiative flux getter
    function get_field_rlwu(g_i, g_j, g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j, g_k
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gathervol(g_i, g_j, g_k, a, n, lwu(2:i1, 2:j1, 1:kmax))
    end function get_field_rlwu

    ! Clear-sky downwelling shortwave radiative flux getter
    function get_field_rswdcs(g_i, g_j, g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j, g_k
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gathervol(g_i, g_j, g_k, a, n, swdca(2:i1, 2:j1, 1:kmax))
    end function get_field_rswdcs

    ! Clear-sky upwelling shortwave radiative flux getter
    function get_field_rswucs(g_i, g_j, g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j, g_k
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gathervol(g_i, g_j, g_k, a, n, swuca(2:i1, 2:j1, 1:kmax))
    end function get_field_rswucs

    ! Clear-sky downwelling longwave radiative flux getter
    function get_field_rlwdcs(g_i, g_j, g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j, g_k
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gathervol(g_i, g_j, g_k, a, n, lwdca(2:i1, 2:j1, 1:kmax))
    end function get_field_rlwdcs

    ! Clear-sky upwelling longwave radiative flux getter
    function get_field_rlwucs(g_i, g_j, g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j, g_k
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gathervol(g_i, g_j, g_k, a, n, lwuca(2:i1, 2:j1, 1:kmax))
    end function get_field_rlwucs
    !!! end of full 3D field getter functions

    ! getter function for LWP - a 2D field, vertical integral of ql
    function get_field_LWP(g_i, g_j, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gatherLWP(g_i, g_j, a, n, ql0(2:i1, 2:j1, 1:kmax))
    end function get_field_LWP

    ! getter function for TWP - Total Water Path - 2D field, vertical integral of qt
    function get_field_TWP(g_i, g_j, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gatherLWP(g_i, g_j, a, n, qt0(2:i1, 2:j1, 1:kmax))
    end function get_field_TWP

    ! Rain water path - like LWP but for rain water
    ! get it from sv0(:,:,:,iqr)
    function get_field_RWP(g_i, g_j, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gatherLWP(g_i, g_j, a, n, sv0(2:i1, 2:j1, 1:kmax, iqr))
    end function get_field_RWP

    ! Surface friction velocity field getter
    function get_field_ustar(g_i, g_j, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gatherlayer(g_i, g_j, a, n, ustar(2:i1, 2:j1))
    end function get_field_ustar

    ! Surface momentum roughness length field getter
    function get_field_z0m(g_i, g_j, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gatherlayer(g_i, g_j, a, n, z0m(2:i1, 2:j1))
    end function get_field_z0m

    ! Surface heat roughness length field getter
    function get_field_z0h(g_i, g_j, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gatherlayer(g_i, g_j, a, n, z0h(2:i1, 2:j1))
    end function get_field_z0h

    ! Surface skin temperature field getter
    function get_field_tskin(g_i, g_j, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gatherlayer(g_i, g_j, a, n, tskin(2:i1, 2:j1))
    end function get_field_tskin

    ! Surface skin moisture content field getter
    function get_field_qskin(g_i, g_j, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gatherlayer(g_i, g_j, a, n, qskin(2:i1, 2:j1))
    end function get_field_qskin

    ! Surface latent heat flux field getter
    function get_field_LE(g_i, g_j, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j
        real, dimension(n), intent(out) :: a
        integer :: ret
        if(allocated(LE)) then
            ret = gatherlayer(g_i, g_j, a, n, LE(2:i1, 2:j1))
        else
            ret = 1
        end if
    end function get_field_LE

    ! Surface sensible heat flux field getter
    function get_field_H(g_i, g_j, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j
        real, dimension(n), intent(out) :: a
        integer :: ret
        if(allocated(H)) then
            ret = gatherlayer(g_i, g_j, a, n, H(2:i1, 2:j1))
        else
            ret = 1
        end if
    end function get_field_H

    ! Surface Obukhov length field getter
    function get_field_obl(g_i, g_j, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gatherlayer(g_i, g_j, a, n, obl(2:i1, 2:j1))
    end function get_field_obl

    ! Surface theta-l flux field getter
    function get_field_thlflux(g_i, g_j, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gatherlayer(g_i, g_j, a, n, thlflux(2:i1, 2:j1))
    end function get_field_thlflux

    ! Surface moisture flux field getter
    function get_field_qtflux(g_i, g_j, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gatherlayer(g_i, g_j, a, n, qtflux(2:i1, 2:j1))
    end function get_field_qtflux

    ! Surface eastward wind vertical gradient getter
    function get_field_dudz(g_i, g_j, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gatherlayer(g_i, g_j, a, n, dudz(2:i1, 2:j1))
    end function get_field_dudz

    ! Surface northward wind vertical gradient getter
    function get_field_dvdz(g_i, g_j, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gatherlayer(g_i, g_j, a, n, dvdz(2:i1, 2:j1))
    end function get_field_dvdz

    ! Surface moisture vertical gradient getter
    function get_field_dqtdz(g_i, g_j, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gatherlayer(g_i, g_j, a, n, dqtdz(2:i1, 2:j1))
    end function get_field_dqtdz

    ! Surface liquid water temperature vertical gradient getter
    function get_field_dthldz(g_i, g_j, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j
        real, dimension(n), intent(out) :: a
        integer :: ret
        ret = gatherlayer(g_i, g_j, a, n, dthldz(2:i1, 2:j1))
    end function get_field_dthldz

    ! Surface prescribed uniform heat flux getter
    function get_wt_surf(wtflux) result(ret)
        real, intent(out) :: wtflux
        real :: ret

        ret = 0
        wtflux = wtsurf !- (thl - tskin) / ra
    end function

    ! Surface prescribed uniform heat flux setter
    function set_wt_surf(wtflux) result(ret)
        real, intent(in) :: wtflux
        real :: ret

        ret = 0
        wtsurf = wtflux !- (thl - tskin) / ra
    end function

    ! Surface prescribed uniform moisture flux getter
    function get_wq_surf(wqflux) result(ret)
        real, intent(out) :: wqflux
        real :: ret

        ret = 0
        wqflux = wqsurf !- (qt - qskin) / ra
    end function

    ! Surface prescribed uniform moisture flux setter
    function set_wq_surf(wqflux) result(ret)
        real, intent(in) :: wqflux
        real :: ret

        ret = 0
        wqsurf = wqflux !- (qt - qskin) / ra
    end function

    ! Surface uniform momentum roughness length getter
    function get_z0m_surf(z0) result(ret)
        real, intent(out) :: z0
        real :: ret

        ret = 0
        z0 = z0mav
    end function

    ! Surface uniform momentum roughness length setter
    function set_z0m_surf(z0) result(ret)
        real, intent(in) :: z0
        real :: ret

        ret = 0
        z0m = z0
        z0mav = z0
    end function

    ! Surface uniform heat roughness length getter
    function get_z0h_surf(z0) result(ret)
        real, intent(out) :: z0
        real :: ret

        ret = 0
        z0 = z0hav
    end function

    ! Surface uniform heat roughness length setter
    function set_z0h_surf(z0) result(ret)
        real, intent(in) :: z0
        real :: ret

        ret = 0
        z0h = z0
        z0hav = z0
    end function

    !!! setter functions for full 3D fields - using index arrays
    !!! these functions set BOTH the -m and the -0 fields
    ! Sets the eastward wind velocity total volume field
    function set_field_U(g_i, g_j, g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j, g_k
        real, dimension(n), intent(in) :: a
        integer :: ret, m, i, j, k
        do m = 1, n
            if (localindex(g_i(m), g_j(m), g_k(m), i, j, k) /= 0) then
                um(i, j, k) = a(m)
                u0(i, j, k) = a(m)
            endif
        enddo
        ret = 0
    end function set_field_U

    ! Sets the northward wind velocity total volume field
    function set_field_V(g_i, g_j, g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j, g_k
        real, dimension(n), intent(in) :: a
        integer :: ret, m, i, j, k

        do m = 1, n
            if (localindex(g_i(m), g_j(m), g_k(m), i, j, k) /= 0) then
                vm(i, j, k) = a(m)
                v0(i, j, k) = a(m)
            endif
        enddo
        ret = 0
    end function set_field_V

    ! Sets the vertical wind velocity total volume field
    function set_field_W(g_i, g_j, g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j, g_k
        real, dimension(n), intent(in) :: a
        integer :: ret, m, i, j, k

        do m = 1, n
            if (localindex(g_i(m), g_j(m), g_k(m), i, j, k) /= 0) then
                wm(i, j, k) = a(m)
                w0(i, j, k) = a(m)
            endif
        enddo
        ret = 0
    end function set_field_W

    ! Sets the liquid water temperature total volume field
    function set_field_THL(g_i, g_j, g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j, g_k
        real, dimension(n), intent(in) :: a
        integer :: ret, m, i, j, k

        do m = 1, n
            if (localindex(g_i(m), g_j(m), g_k(m), i, j, k) /= 0) then
                thlm(i, j, k) = a(m)
                thl0(i, j, k) = a(m)
            endif
        enddo
        ret = 0
    end function set_field_THL

    ! Sets the total humidity total volume field
    function set_field_QT(g_i, g_j, g_k, a, n) result(ret)
        integer, intent(in) :: n
        integer, dimension(n), intent(in) :: g_i, g_j, g_k
        real, dimension(n), intent(in) :: a
        integer :: ret, m, i, j, k

        do m = 1, n
            if (localindex(g_i(m), g_j(m), g_k(m), i, j, k) /= 0) then
                qtm(i, j, k) = a(m)
                qt0(i, j, k) = a(m)
            endif
        enddo
        ret = 0
    end function set_field_QT
    !!! end of setter functions for 3D fields

    !!! get simulation parameters
    function get_params_grid(i, j, k, X, Y) result(ret)
        use modglobal, only : itot, jtot, kmax, xsize, ysize
        implicit none
        integer, intent(out) :: i, j, k
        real, intent(out) :: X, Y
        integer :: ret

        i = itot
        j = jtot
        k = kmax
        X = xsize
        Y = ysize
        ret = 0
    end function get_params_grid

    ! Write a restart file
    function write_restart() result(ret)
        integer :: ret
        call do_writerestartfiles

        ret = 0
    end function write_restart


end module dales_interface
