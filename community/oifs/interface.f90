! This is the omuse interface to OpenIFS. Written by G. van den Oord, 2016

module openifs_interface

    use ifslib,     only: static_init, initialize, finalize, get_gp_geom,& 
                        & jstep, step, PFULL, PHALF, WIND_U, WIND_V, TEMPERATURE,&
                        & SPEC_HUMIDITY, ICE_WATER, LIQ_WATER, CLOUD_FRACTION,&
                        & OZONE
    use spcrmlib,   only: step_until_cloud_scheme,step_cloud_scheme,step_from_cloud_scheme,&
                        & gettends,get_cur_field,settend,set_sp_mask,spcrminit,spcrmclean,&
                        & tendcml,settend
    use yomct0,     only: nstart,nstop
    use yomct3,     only: nstep
    use yomdyn,     only: tstep
    use yommp0,     only: myproc
    use yomgem,     only: ngptotg
    use yomdimv,    only: nflevg
    USE yomlun,     only: nulout
    
    implicit none

    contains

        function initialize_code() result(ret)

            integer:: ret

            call static_init()
            call initialize('TEST')
            call spcrminit()

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

        function evolve_model_single_step() result(ret)

          integer::               ret
          logical::               istat

          ret = 0
          istat = .true.
          call step(istat)
          if(.not.istat) then
             ret = 1
          endif

        end function evolve_model_single_step

        function evolve_model_until_cloud_scheme() result(ret)

          integer::               ret
          logical::               istat

          ret = 0
          istat = .true.
          call step_until_cloud_scheme(istat)
          if(.not.istat) then
             ret = 1
          endif

        end function evolve_model_until_cloud_scheme

        function evolve_model_from_cloud_scheme() result(ret)

          integer::               ret
          logical::               istat

          ret = 0
          istat = .true.
          call step_from_cloud_scheme(istat)
          if(.not.istat) then
             ret = 1
          endif

        end function evolve_model_from_cloud_scheme

        function evolve_model_cloud_scheme() result(ret)

          integer::               ret
          logical::               istat

          ret = 0
          istat = .true.
          call step_cloud_scheme()
          if(.not.istat) then
             ret = 1
          endif

        end function evolve_model_cloud_scheme
          
        function get_gridpoints(g_i,lats,lons,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i
            real(8), dimension(n), intent(out)::    lats(n),lons(n)
            real(8), allocatable::                  b(:),c(:)
            integer::                               ret,i

            ret = 0

            if(myproc == 1) then
                allocate(b(ngptotg),c(ngptotg),stat=ret)
                if(ret/=0) then
                    return
                endif
            endif

            call get_gp_geom(b,c)
            
            if(myproc == 1) then
                do i = 1,n
                    lats(i) = b(g_i(i) + 1)
                    lons(i) = c(g_i(i) + 1)
                enddo
            endif

        end function

        function get_field_Pfull_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_field(g_i,g_k,a,n,PFULL)

        end function get_field_Pfull_

        function get_field_Phalf_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_field(g_i,g_k,a,n,PHALF)

        end function get_field_Phalf_

        function get_field_U_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_field(g_i,g_k,a,n,WIND_U)

        end function get_field_U_

        function get_field_V_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_field(g_i,g_k,a,n,WIND_V)

        end function get_field_V_

        function get_field_T_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_field(g_i,g_k,a,n,TEMPERATURE)

        end function get_field_T_

        function get_field_SH_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_field(g_i,g_k,a,n,SPEC_HUMIDITY)

        end function get_field_SH_

        function get_field_QL_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_field(g_i,g_k,a,n,LIQ_WATER)

        end function get_field_QL_

        function get_field_QI_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_field(g_i,g_k,a,n,ICE_WATER)

        end function get_field_QI_

        function get_field_A_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_field(g_i,g_k,a,n,CLOUD_FRACTION)

        end function get_field_A_

        function get_field_O3_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_field(g_i,g_k,a,n,OZONE)

        end function get_field_O3_

        function get_field(g_i,g_k,a,n,fldid) result(ret)

            integer, intent(in)::                   n,fldid
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            real(8), allocatable::                  b(:,:)
            integer::                               ret,ii

            ret = 0
            if(myproc == 1) then
                if(fldid == PHALF) then
                    allocate(b(ngptotg,nflevg+1),stat=ret)
                else
                    allocate(b(ngptotg,nflevg),stat=ret)
                endif
                if(ret/=0) then
                    return
                endif
            endif
            
            call get_cur_field(b,fldid)
            
            if(myproc == 1) then
                do ii = 1,n
                    a(ii) = b(g_i(ii) + 1,g_k(ii) + 1)
                enddo
                deallocate(b)
            endif

        end function get_field

        function get_tendency_U_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_tendency(g_i,g_k,a,n,WIND_U)

        end function get_tendency_U_

        function get_tendency_V_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_tendency(g_i,g_k,a,n,WIND_V)

        end function get_tendency_V_

        function get_tendency_T_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_tendency(g_i,g_k,a,n,TEMPERATURE)

        end function get_tendency_T_

        function get_tendency_SH_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_tendency(g_i,g_k,a,n,SPEC_HUMIDITY)

        end function get_tendency_SH_

        function get_tendency_QL_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_tendency(g_i,g_k,a,n,LIQ_WATER)

        end function get_tendency_QL_

        function get_tendency_QI_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_tendency(g_i,g_k,a,n,ICE_WATER)

        end function get_tendency_QI_

        function get_tendency_A_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_tendency(g_i,g_k,a,n,CLOUD_FRACTION)

        end function get_tendency_A_

        function get_tendency_O3_(g_i,g_k,a,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            integer::                               ret

            ret = get_tendency(g_i,g_k,a,n,OZONE)

        end function get_tendency_O3_

        function get_tendency(g_i,g_k,a,n,fldid) result(ret)

            integer, intent(in)::                   n,fldid
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(out)::    a(n)
            real(8), allocatable::                  b(:,:)
            integer::                               ret,ii

            ret = 0
            if(myproc == 1) then
                if(fldid == PHALF) then
                    allocate(b(ngptotg,nflevg+1),stat=ret)
                else
                    allocate(b(ngptotg,nflevg),stat=ret)
                endif
                if(ret/=0) then
                    return
                endif
            endif
            
            call gettends(b,fldid,tendcml)
            
            if(myproc == 1) then
                do ii = 1,n
                    a(ii) = b(g_i(ii) + 1,g_k(ii) + 1)
                enddo
                deallocate(b)
            endif

        end function get_tendency

        function set_mask(g_i,n) result(ret)

            integer, dimension(n), intent(in)::     g_i
            integer, intent(in)::                   n
            integer::                               ret

            ret = 0
            write (NULOUT,*) "calling set_sp_mask", g_i
            call set_sp_mask(g_i+1,n,.TRUE.)
            ! convert to Fortran indexing
        end function set_mask

        function reset_mask(g_i,n) result(ret)

            integer, dimension(n), intent(in)::     g_i
            integer, intent(in)::                   n
            integer::                               ret

            ret = 0
            call set_sp_mask(g_i+1,n,.FALSE.)
            ! convert to Fortran indexing
        end function reset_mask

        function set_tendency_U_(g_i,g_k,vals,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(in)::     vals
            integer::                               ret

            ret = set_tendency(g_i,g_k,vals,n,WIND_U)

        end function set_tendency_U_

        function set_tendency_V_(g_i,g_k,vals,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(in)::     vals
            integer::                               ret

            ret = set_tendency(g_i,g_k,vals,n,WIND_V)

        end function set_tendency_V_

        function set_tendency_T_(g_i,g_k,vals,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(in)::     vals
            integer::                               ret

            ret = set_tendency(g_i,g_k,vals,n,TEMPERATURE)

        end function set_tendency_T_

        function set_tendency_SH_(g_i,g_k,vals,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(in)::     vals
            integer::                               ret

            ret = set_tendency(g_i,g_k,vals,n,SPEC_HUMIDITY)

        end function set_tendency_SH_

        function set_tendency_QL_(g_i,g_k,vals,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(in)::     vals
            integer::                               ret

            ret = set_tendency(g_i,g_k,vals,n,LIQ_WATER)

        end function set_tendency_QL_

        function set_tendency_QI_(g_i,g_k,vals,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(in)::     vals
            integer::                               ret

            ret = set_tendency(g_i,g_k,vals,n,ICE_WATER)

        end function set_tendency_QI_

        function set_tendency_A_(g_i,g_k,vals,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(in)::     vals
            integer::                               ret

            ret = set_tendency(g_i,g_k,vals,n,CLOUD_FRACTION)

        end function set_tendency_A_

        function set_tendency_O3_(g_i,g_k,vals,n) result(ret)

            integer, intent(in)::                   n
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(in)::     vals
            integer::                               ret

            ret = set_tendency(g_i,g_k,vals,n,OZONE)

        end function set_tendency_O3_

        function set_tendency(g_i,g_k,vals,n,fldid) result(ret)

            integer, intent(in)::                   n,fldid
            integer, dimension(n), intent(in)::     g_i,g_k
            real(8), dimension(n), intent(in)::     vals
            integer::                               ret,i,nvals
            integer, allocatable::                  klevs(:)

            ret = 0
            nvals = 0

            allocate(klevs(nflevg))

            ! convert from sequence of point to setting vertical profiles
            ! needs testing with multiple 
            do i = 1,n
                nvals = nvals + 1
                if(nvals > nflevg) then
                    ret = 1
                    return
                endif
                klevs(nvals) = g_k(i) + 1 ! to Fortran indexing (g_k)
                if(i == n) then
                    call settend(g_i(i)+1,klevs(1:nvals),nvals,vals((i - nvals + 1):i),fldid)
                    nvals = 0           ! to fortran indexing (g_i)
                elseif(g_i(i + 1) /= g_i(i)) then
                    call settend(g_i(i)+1,klevs(1:nvals),nvals,vals((i - nvals + 1):i),fldid)
                    nvals = 0
                endif
            enddo

            deallocate(klevs)

        end function set_tendency

        function cleanup_code() result(ret)

            integer:: ret

            ret = 0
            call finalize()
            call spcrmclean()

        end function

end module openifs_interface
