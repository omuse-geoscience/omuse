      module omuse_qgcm
      
      USE parameters
      USE atconst
      USE occonst
      USE radiate
      USE timinfo
      USE monitor
#ifndef ocean_only
      USE athomog
      USE amlsubs
      USE qgasubs
      USE atisubs
#endif
#if defined tau_udiff || !defined ocean_only
      USE atstate
#endif
#ifndef atmos_only
      USE ochomog
      USE ocstate
      USE omlsubs
      USE qgosubs
      USE ocisubs
#endif
      USE intrfac
      USE conhoms
      USE vorsubs
      USE xfosubs, ONLY : xforc
      USE eigmode
      USE omuse_valsubs
      USE radsubs
      USE topsubs
#ifdef use_netcdf
      USE nc_subs
#endif /* use_netcdf */
*     Diagnostics modules
      USE mondiag
#ifdef get_areav
      USE areasubs
#endif
#ifndef atmos_only
#  ifdef qoc_diag
      USE qocdiag
#  endif /* qoc_diag */
#endif /* atmos_only */
#ifdef get_covar
      USE covaria
#endif /* get_covar */
      USE timavge

      IMPLICIT NONE

#ifdef use_netcdf
      INCLUDE 'netcdf.inc'
#endif /* use_netcdf */

*
*     Local parameters
      character (len=*), parameter :: subnam = 'omuse_q-gcm'
*     Define I/O unit numbers:
*     ipunit is for input.params
*     odunit is for output directory
      integer ipunit, odunit
      parameter ( ipunit=44, odunit=45 )
*
      DOUBLE PRECISION PI,TWOPI,PIBY2
      PARAMETER ( PI=3.14159265358979324D0,
     &            TWOPI=6.28318530717958648D0,
     &            PIBY2=1.57079632679489662D0 )
*
*     Arrays passed as arguments
*
*     Local variables
      integer i,j,k,i1,iorest,nsko,nska,nt,lenod,lename,oceanonly,
     &        atmosonly,getcovar,outfloc(7),outflat(7),
     &        cyclicoc,hflxsb,hflxnb,tauudiff
      double precision tspind
      character (len=72) :: name
      character (len=80) :: outdir, inpbuf
      logical resdmp,solnok,avgatm,avgocn
*     Extra variables for topography
      character (len=80) :: topocname,topatname

#ifdef use_netcdf
*     Extra variables for netCDF input/output
      integer status,tempid,varid
#  ifndef atmos_only
      double precision poc4dsz
#  endif /* not atmos_only */
#  ifndef ocean_only
      double precision pat4dsz
#  endif /* not ocean_only */
#else
#  ifdef atmos_only
      character (len=*), parameter :: atforav='./atforce.avg.binary'
#  endif /* atmos_only */
#  ifdef ocean_only
      character (len=*), parameter :: ocforav='./ocforce.avg.binary'
#  endif /* ocean_only */
      logical fexist
#endif

*     Extra variables for OpenMP
!$    integer nprocs,nthmax,numthr,thrnum
!$    logical dynadj,nested

!$    integer, external :: OMP_GET_NUM_PROCS,OMP_GET_MAX_THREADS,
!$   &        OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM
!$    logical, external :: OMP_GET_DYNAMIC,OMP_GET_NESTED


! auxiliary omuse variables
      double precision :: tocean, tos, tosm! ocean state timers
      double precision :: tatm, tas, tasm! atmos. state timers
#ifndef atmos_only
      double precision :: dpo_dt(nxpo,nypo,nlo) ! auxiliary for pom
#endif
#ifndef ocean_only
      double precision :: dpa_dt(nxpa,nypa,nla) ! auxiliary for pam
#endif
      double precision :: dsst_dt(nxto,nyto) ! auxiliary for sstm
      double precision :: dhmixa_dt(nxta,nyta) ! auxiliary for hmixa
      double precision :: dast_dt(nxta,nyta) ! auxiliary for ast

#ifdef ocean_only
      double precision :: pa(1,1,1) ! dummies because interface.f90 expects them 
      double precision :: pam(1,1,1)
      double precision :: qa(1,1,1) ! 
      double precision :: dpa_dt(1,1,1) 

#endif       
#ifdef atmos_only
      double precision :: po(1,1,1) ! 
      double precision :: pom(1,1,1) ! 
      double precision :: qo(1,1,1) ! 
      double precision :: dpo_dt(1,1,1) !
#endif       


      contains

      subroutine omuse_initialize()
      
      dpo_dt=0
      dpa_dt=0
      dsst_dt=0
      dhmixa_dt=0
      dast_dt=0
            
      call initialize()
      
      end subroutine
      
      subroutine initialize()
      
      print *,' '
      write(*,*) ' Computational parameters:'
      write(*,*) ' ========================='
      write(*,*) ' Model version is Q-GCM v1.5.0'
*     Compute typical array sizes, and local arrays for stack allocation
      call memreq

*     Examine OpenMP environment
*     --------------------------
!$    nprocs = OMP_GET_NUM_PROCS()
!$    nthmax = OMP_GET_MAX_THREADS()
!$    dynadj = OMP_GET_DYNAMIC()
!$    nested = OMP_GET_NESTED()
!$    write(*,*) ' '
!$    write(*,*) ' OpenMP parallelism is activated'
!$    write(*,*) ' -------------------------------'
!$    write(*,'(a,i5)') '  No. of processors available = ',nprocs
!$    write(*,'(a,i3)') '  Max. no. of threads available = ',nthmax
!$    write(*,*) ' Dynamic adjustment = ',dynadj
!$    write(*,*) ' Nested parallelism = ',nested
*     Test OpenMP is functioning correctly
!$    write(*,*) ' Outside parallel section:'
!$    numthr = OMP_GET_NUM_THREADS()
!$    write(*,*) ' Number of threads = ',numthr
!$    thrnum = OMP_GET_THREAD_NUM()
!$    write(*,'(a,i4)') '  thrnum = ',thrnum
!$    write(*,*) ' Test trivial parallel loop:'
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (j,thrnum)
!$OMP&         SHARED (numthr)
!$OMP SINGLE
!$    numthr = OMP_GET_NUM_THREADS()
!$    write(*,*) ' Number of threads = ',numthr
!$OMP END SINGLE
!$OMP DO SCHEDULE (STATIC)
!$    do j=1,numthr
!$      thrnum = OMP_GET_THREAD_NUM()
!$      write(*,'(a,i6,i4)') '  j, thrnum = ',j,thrnum
!$    enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL
*     The following loop is the most realistic test of OpenMP //ism
*     Enable (if necessary) by removing the first ! from each line
*     Comment out again once OpenMP is working correctly,
*     to avoid lots of spurious output
!!$    write(*,*) ' Test typical j-range loop (ocean p-points)'
!!$OMP PARALLEL DO DEFAULT (NONE)
!!$OMP&         PRIVATE (j,thrnum)
!!$OMP&         SCHEDULE (STATIC)
!!$    do j=1,nypo
!!$      thrnum = OMP_GET_THREAD_NUM()
!!$      write(*,'(a,i6,i4)') '  j, thrnum = ',j,thrnum
!!$    enddo
!!$OMP END PARALLEL DO

*     Read input parameters via included file
*     skipped for OMUSE
*      INCLUDE './in_param.f'

*     Specify configuration being used
*     --------------------------------
*     Specify as a set of integers so that the config-
*     uration information can be conveniently written
*     to a Matlab-compatible dump via out_param.f
#ifdef ocean_only
      oceanonly = 1
#else
      oceanonly = 0
#endif
#ifdef atmos_only
      atmosonly = 1
#else
      atmosonly = 0
#endif /* atmos_only */
#ifdef get_covar
      getcovar = 1
#else
      getcovar = 0
#endif /* get_covar */
#ifdef cyclic_ocean
      cyclicoc = 1
#else
      cyclicoc = 0
#endif
#ifdef sb_hflux
      hflxsb = 1
#else
      hflxsb = 0
#endif
#ifdef nb_hflux
      hflxnb = 1
#else
      hflxnb = 0
#endif
#ifdef tau_udiff
      tauudiff = 1
#else
      tauudiff = 0
#endif

      end subroutine initialize
      
      subroutine check_parameters

      print *,' '
      print *,' Control parameters:'
      print *,' ==================='
      lenod = index(outdir, '   ') - 1
      lename = index(name, '   ') - 1
      print *,' outdir = ',outdir(1:lenod)
      print *,' name = ',name(1:lename)
      if ( oceanonly.eq.1 .and. atmosonly.eq.1 ) then
        print *,' '
        print *,' Invalid model configuration: ocean_only and'
        print *,' atmos_only options cannot both be selected'
        print *,' Program terminates'
        stop
       else if ( oceanonly.eq.1 ) then
        print *,' Model running in ocean_only configuration'
       else if ( atmosonly.eq.1 ) then
        print *,' Model running in atmos_only configuration'
       else
        print *,' Model running in full coupled configuration'
      endif
      if ( atmosonly.eq.0 ) then
#ifdef cyclic_ocean
        print *,' Model running in cyclic ocean configuration'
#else
        print *,' Model running in finite box ocean configuration'
#endif
        if ( oceanonly.eq.0 ) then
*         tau-udiff option only applies in coupled mode
*         i.e. neither atmos-only nor ocean_only are set
#ifdef tau_udiff
          print *,' Windstress depends on atmos.-ocean vel. diff.'
#else
          print *,' Windstress depends on atmos. velocity only'
#endif
        endif
      endif
#ifdef get_areav
      print *,' Model computes area average statistics'
#else
      print *,' Model does not compute area average stats'
#endif
#ifdef qoc_diag
      print *,' Model computes ocean dq/dt diagnostics '
#else
      print *,' Model does not compute ocean dq/dt diagnostics'
#endif
      if ( getcovar.eq.1 ) then
        print *,' Model accumulates covariance statistics'
        write(*,201) '  Sampling intervals nscvat, nscvoc = ',
     &               nscvat,nscvoc
       else
        print *,' Model does not accumulate covariance stats'
      endif
      if ( hflxsb.eq.1 .and. hflxnb.eq.1 ) then
        print *,' '
        print *,' Invalid model configuration: sb_hflux and'
        print *,' nb_hflux options cannot both be selected'
        print *,' Program terminates'
        stop
       else if ( hflxsb.eq.1 ) then
        print *,' Model running with modified o.m.l. southern b.c.'
       else if ( hflxnb.eq.1 ) then
        print *,' Model running with modified o.m.l. northern b.c.'
       else
        print *,' Model running with no heat flux N & S'
      endif
*     Check some of the grid parameters
      if ( nla.lt.2 .or. nlo.lt.2 ) then
        print *,' '
        print *,' nla, nlo = ',nla,nlo
        print *,' Inadequate nla or nlo, needs to be at least 2'
        print *,' Program terminates'
        stop
      endif
#ifdef cyclic_ocean
*     Check cyclic ocean case is properly set up
      if ( nxta.ne.nxaooc ) then
        print *,' '
        print *,' nxta, nxaooc = ',nxta,nxaooc
        print *,' For cyclic ocean nxta should equal nxaooc'
        print *,' Program terminates'
        stop
      endif
#else
*     Check finite box ocean case is properly set up
      if ( nxta.lt.nxaooc ) then
        print *,' '
        print *,' nxta, nxaooc = ',nxta,nxaooc
        print *,' Inadequate nxta needs to be at least nxaooc'
        print *,' Program terminates'
        stop
      endif
#endif
      if ( nyta.lt.nyaooc ) then
        print *,' '
        print *,' nyta, nyaooc = ',nyta,nyaooc
        print *,' Inadequate nyta needs to be at least nyaooc'
        print *,' Program terminates'
        stop
      endif
#ifdef get_covar
      if ( mod(nxta,nscvat).ne.0 .or. mod(nyta,nscvat).ne.0 ) then
        print *,' '
        print *,' nxta, nyta, nscvat = ',nxta,nyta,nscvat
        print *,' Need nxta, nyta both to be integer multiples of'
        print *,' nscvat for covariance subsampling to work correctly'
        print *,' Program terminates'
        stop
      endif
      if ( mod(nxto,nscvoc).ne.0 .or. mod(nyto,nscvoc).ne.0 ) then
        print *,' '
        print *,' nxto, nyto, nscvoc = ',nxto,nyto,nscvoc
        print *,' Need nxto, nyto both to be integer multiples of'
        print *,' nscvoc for covariance subsampling to work correctly'
        print *,' Program terminates'
        stop
      endif
#endif /* get_covar */
*     Check heat-flux boundary conditions
#ifdef sb_hflux
      if ( fnot.lt.0.0d0 ) then
        print *,' '
        print *,' Southern boundary o.m.l. heat-flux activated'
        print *,' Sign of fnot -> running in southern hemisphere'
        print *,' These are inconsistent choices'
        print *,' Program terminates'
        stop
      endif
#endif
#ifdef nb_hflux
      if ( fnot.gt.0.0d0 ) then
        print *,' '
        print *,' Northern boundary o.m.l. heat-flux activated'
        print *,' Sign of fnot -> running in northern hemisphere'
        print *,' These are inconsistent choices'
        print *,' Program terminates'
        stop
      endif
#endif

*     Derive atmos gridspace and ocean timestep
*     -----------------------------------------
*     Derive larger from smaller to ensure integer ratio
      dxa = ndxr*dxo
      dto = nstr*dta
      write(*,201) '  Atmos/ocean grid ratio ndxr = ',ndxr
      write(*,201) '  Oc/atm. timestep ratio nstr = ',nstr
      write(*,201) '  Atmos. gridcells over ocean = ',nxaooc,nyaooc
      write(*,201) '  Ocn start indices  nx1, ny1 = ',nx1,ny1
      write(*,215) '  Coriolis par. f0 (rad s^-1) = ',fnot
      write(*,215) '  Beta =df/dy (rad s^-1 m^-1) = ',beta

*     Atmospheric grid
*     ----------------
      dya = dxa
      hdxam1 = 0.5d0/dxa
      dxam2 = 1.0d0/(dxa*dxa)
      xla = nxta*dxa
      yla = nyta*dya
      do i=1,nxpa
        xpa(i) = (i-1)*dxa
      enddo
      do i=1,nxta
        xta(i) = xpa(i) + 0.5d0*dxa
      enddo
      do j=1,nypa
        ypa(j) = (j-1)*dya
        yparel(j) = ypa(j) - 0.5d0*yla
      enddo
      do j=1,nyta
        yta(j) = ypa(j) + 0.5d0*dya
        ytarel(j) = yta(j) - 0.5d0*yla
      enddo

*     Oceanic grid
*     ------------
      dyo = dxo
      hdxom1 = 0.5d0/dxo
      dxom2 = 1.0d0/(dxo*dxo)
      xlo = nxto*dxo
      ylo = nyto*dyo
      do i=1,nxpo
        xpo(i) = (i-1)*dxo + (nx1-1)*dxa
      enddo
      do i=1,nxto
        xto(i) = xpo(i) + 0.5d0*dxo
      enddo
      do j=1,nypo
        ypo(j) = (ny1-1)*dya + (j-1)*dyo
        yporel(j) = ypo(j) - 0.5d0*yla
      enddo
      do j=1,nyto
        yto(j) = ypo(j) + 0.5d0*dyo
        ytorel(j) = yto(j) - 0.5d0*yla
      enddo

*     Derive useful constants
*     -----------------------
      rdxaf0 = 1.0d0/(dxa*fnot)
      rdxof0 = 1.0d0/(dxo*fnot)
      rrcpat = 1.0d0/(rhoat*cpat)
      rrcpoc = 1.0d0/(rhooc*cpoc)
      raoro = rhoat/rhooc
      tdto = 2.0d0*dto
      tdta = 2.0d0*dta

*     Derive total thicknesses
*     ------------------------
      hto = 0.0d0
      do k=1,nlo
        hto = hto + hoc(k)
      enddo
      hta = 0.0d0
      do k=1,nla
        hta = hta + hat(k)
      enddo

*     Print out a few interesting numbers
*     -----------------------------------
      call eigmod (nlo, gpoc, hoc, 'Ocean', amatoc,
     &             cphsoc, rdefoc, rdm2oc, ctl2moc, ctm2loc)
      print *,' '
      print *,' Oceanic parameters:'
      print *,' -------------------'
      write(*,201) '  No. of ocean QG layers  nlo = ',nlo
      write(*,201) '  No. of gridcells nxto, nyto = ',nxto,nyto
      write(*,204) '  Gridlength dxo         (km) = ',1.0d-3*dxo
      write(*,203) '  Domain sizes xlo, ylo  (km) = ',
     &             1.0d-3*xlo,1.0d-3*ylo
      write(*,206) '  Rossby number   Beta*ylo/f0 = ',beta*ylo/abs(fnot)
      write(*,215) '  f range S -> N   (rad s^-1) = ',
     &             fnot+beta*yporel(1),fnot+beta*yporel(nypo)
      write(*,215) '  Midlatitude Coriolis param  = ',
     &             fnot+beta*0.5d0*( yporel(1) + yporel(nypo) )
      write(*,205) '  Timestep dto      (minutes) = ',dto/60.0d0
      write(*,201) '  No. of timesteps per day    = ',nint(secday/dto)
      write(*,203) '  Mixed layer thickness   (m) = ',hmoc
      write(*,213) '  Sp. ht. cap. (J kg^-1 K^-1) = ',cpoc
      write(*,213) '  SST grad-2 diff  (m^2 s^-1) = ',st2d
      call diffts (2, nlo, st2d, 1, dxo, rdefoc, dto)
      write(*,213) '  SST grad-4 diff  (m^4 s^-1) = ',st4d
      call diffts (4, nlo, st4d, 1, dxo, rdefoc, dto)
      write(*,205) '  Abs. pot. temp. tabsoc  (K) = ',
     &             (tabsoc(k),k=1,nlo)
      write(*,203) '  Layer thicknesses hoc   (m) = ',(hoc(k),k=1,nlo)
      write(*,203) '  Total thickness   hto   (m) = ',hto
      write(*,207) '  Reduced gravities  (m s^-2) = ',
     &             (gpoc(k),k=1,nlo-1)
      write(*,206) '  Baroclinic wavespeeds (m/s) = ',
     &             (cphsoc(k),k=2,nlo)
      write(*,206) '  Courant number(s)           = ',
     &             ( (dto/dxo)*cphsoc(k),k=2,nlo)
      write(*,205) '  Deformation radii      (km) = ',
     &             (1.0d-3*rdefoc(k),k=2,nlo)
      write(*,206) '  Gridlengths per def. radius = ',
     &             (rdefoc(k)/dxo,k=2,nlo)
      write(*,207) '  Long Rossby wavespeed (m/s) = ',
     &             (-beta*cphsoc(k)**2/fnot**2,k=2,nlo)
      write(*,213) '  Del-sqd coeffts  (m^2 s^-1) = ',(ah2oc(k),k=1,nlo)
      call diffts (2, nlo, ah2oc, nlo, dxo, rdefoc, dto)
      write(*,213) '  Del-4th coeffts  (m^4 s^-1) = ',(ah4oc(k),k=1,nlo)
      call diffts (4, nlo, ah4oc, nlo, dxo, rdefoc, dto)
      write(*,204) '  Munk b.l. width scale  (km) = ',
     &             (1.0d-3*(ah4oc(k)/beta)**0.2d0,k=1,nlo)
      write(*,204) '  Munk b.l. width scale (pts) = ',
     &             (((ah4oc(k)/beta)**0.2d0)/dxo,k=1,nlo)
      write(*,204) '  Bottom Ekm. layer thickness = ',delek
      write(*,213) '  Bottom layer Ekman number   = ',
     &             (delek/hoc(nlo))**2
      if ( delek.lt.0.0d0 ) then
        print *,' Invalid -ve value of delek'
        print *,' Program terminates'
        stop
       else if (delek.eq.0.0d0 ) then
        tspind = 0.0d0
       else
        tspind = 2.0d0*hoc(nlo)/(abs(fnot)*delek)/secday
      endif
      write(*,204) '  Spindown timescale   (days) = ',tspind
      write(*,213) '  Mixed BC coeff. bccooc (nd) = ',bccooc

      call eigmod (nla, gpat, hat, 'Atmosphere', amatat,
     &             cphsat, rdefat, rdm2at, ctl2mat, ctm2lat)
      print *,' '
      print *,' Atmospheric parameters:'
      print *,' -----------------------'
      write(*,201) '  No. of atmos. QG layers nla = ',nla
      write(*,201) '  No. of gridcells nxta, nyta = ',nxta,nyta
      write(*,201) '  At ocean res., no. of cells = ',nxtaor,nytaor
      write(*,204) '  Gridlength dxa         (km) = ',1.0d-3*dxa
      write(*,203) '  Domain sizes xla, yla  (km) = ',
     &             1.0d-3*xla,1.0d-3*yla
      write(*,206) '  Rossby number   Beta*yla/f0 = ',beta*yla/abs(fnot)
      write(*,215) '  f range S -> N   (rad s^-1) = ',
     &             fnot+beta*yparel(1),fnot+beta*yparel(nypa)
      write(*,215) '  Midlatitude Coriolis param  = ',
     &             fnot+beta*0.5d0*( yparel(1) + yparel(nypa) )
      write(*,205) '  Timestep dta      (minutes) = ',dta/60.0d0
      write(*,201) '  No. of timesteps per day    = ',nint(secday/dta)
      write(*,203) '  Mixed layer thickness   (m) = ',hmat
      write(*,203) '  Min. mixed layer thick. (m) = ',hmamin
      write(*,213) '  Quad. drag coefft Cd (nond) = ',cdat
      write(*,213) '  Sp. ht. cap. (J kg^-1 K^-1) = ',cpat
      write(*,213) '  AST grad-2 diff  (m^2 s^-1) = ',at2d
      call diffts (2, nla, at2d, 1, dxa, rdefat, dta)
      write(*,213) '  AST grad-4 diff  (m^4 s^-1) = ',at4d
      call diffts (4, nla, at4d, 1, dxa, rdefat, dta)
      write(*,213) '  hmix diffusivity (m^2 s^-1) = ',ahmd
      call diffts (2, nla, ahmd, 1, dxa, rdefat, dta)
      write(*,213) '  hmix damping coefft  hmadmp = ',hmadmp
      write(*,203) '  Layer thicknesses hat   (m) = ',(hat(k),k=1,nla)
      write(*,203) '  Total thickness   hta   (m) = ',hta
      write(*,205) '  Abs. pot. temp. tabsat  (K) = ',
     &             (tabsat(k),k=1,nla)
      write(*,207) '  Reduced gravities  (m s^-2) = ',
     &             (gpat(k),k=1,nla-1)
      write(*,206) '  Baroclinic wavespeeds (m/s) = ',
     &             (cphsat(k),k=2,nla)
      write(*,206) '  Courant number(s)           = ',
     &             ( (dta/dxa)*cphsat(k),k=2,nla)
      write(*,205) '  Deformation radii      (km) = ',
     &             (1.0d-3*rdefat(k),k=2,nla)
      write(*,206) '  Gridlengths per def. radius = ',
     &             (rdefat(k)/dxa,k=2,nla)
      write(*,213) '  Del-4th coeffts  (m^4 s^-1) = ',(ah4at(k),k=1,nla)
      call diffts (4, nla, ah4at, nla, dxa, rdefat, dta)
      write(*,213) '  Mixed BC coeff. bccoat (nd) = ',bccoat

      print *,' '
      print *,' Coupling parameters:'
      print *,' --------------------'
      write(*,205) '  Coefft. Lambda   (W m^-2/K) = ',xlamda
      write(*,204) '  Ast coupling  coefft  xcexp = ',xcexp
      write(*,204) '  Sst advection coefft  ycexp = ',ycexp

  201 format(a,9i13)
  203 format(a,9f13.3)
  204 format(a,9f13.4)
  205 format(a,9f13.5)
  206 format(a,9f13.6)
  207 format(a,9f13.7)
  213 format(a,1p,9d13.3)
* 214 format(a,1p,9d13.4)
  215 format(a,1p,9d13.5)
* 225 format(a,i2,a,9f13.5)
* 226 format(a,i2,a,9f13.6)

      end subroutine check_parameters

      subroutine init_topo_rad_state

*     Topography section
*     ==================
      call topset (topocname, topatname, outdir, lenod, .true.)

*     Radiation section
*     =================
*     Compute mean state radiative balance and perturbation
*     radiation coefficients A, B, C and D. Also compute atmosphere
*     and ocean mixed layer temperatures that ensure equilibrium
      call radiat

*     Initialise pressure and temperature fields
*     ==========================================
      if ( name.eq.'zero' ) then
        call zeroin
        tini = 0.0d0
       else if ( name.eq.'rbal' ) then
        call rbalin
        tini = 0.0d0
       else
#ifdef use_netcdf
!! Read netCDF restart dump
        status = nf_open(name, NF_NOWRITE, restid)
        if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
        call restart_nc
        status = nf_close(restid)
        if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#else
*       Read unformatted restart dump
        open (11, file=name, form='unformatted', status='old')

c   time stamp
        read(11) tini

c   ocean pressure first
#  ifndef atmos_only
        read(11) po,pom
#  endif /* atmos_only */

c   initial condition for pa
#  ifndef ocean_only
        read(11) pa,pam
#  endif

c   sst initial condition
        read(11) sst,sstm

c   atmospheric surface temperature
        read(11) ast,astm

c   atmospheric mixed layer thickness
        read(11) hmixa,hmixam

        close (11)
#endif
      endif

      end subroutine init_topo_rad_state
      
      subroutine init_grids1

*     Derive run length and interval parameters
*     -----------------------------------------
      print *,' '
      print *,' Run length and interval parameters:'
      print *,' -----------------------------------'
      tend = tini + trun
      write(*,204) '  Start time tini     (years) = ',tini
      write(*,204) '  Run length trun     (years) = ',trun
      write(*,204) '  Final time tend     (years) = ',tend
      nsteps0 = nint(tini*secsyr/dta)
      nsteps = nint(tend*secsyr/dta)
      ntsrun = nsteps - nsteps0
      write(*,201) '  Start no. of (atmos)  steps = ',nsteps0
      write(*,201) '  Final no. of (atmos)  steps = ',nsteps
      write(*,201) '  Total no. of (atmos)  steps = ',ntsrun
*     Do all dumps when atmos. and ocean are in step
      nvalid = nint( valday*secday/dto )*nstr
      noutoc = nint( odiday*secday/dto )*nstr
      noutat = nint( adiday*secday/dto )*nstr
      write(*,204) '  Validity test int.   (days) = ',valday
      write(*,204) '  Ocean data dump int. (days) = ',odiday
      write(*,204) '  Atmos data dump int. (days) = ',adiday
      write(*,204) '  Diagnostics dump in. (days) = ',dgnday
      write(*,204) '  Printout dump inter. (days) = ',prtday
      resdmp = resday.gt.0.0d0
      if ( resdmp ) then
        write(*,204) '  Restart dump interval (day) = ',resday
        noutre = nint( resday*secday/dto )*nstr
        iorest = 12
      endif
      nocmon = nint( dgnday*secday/dto )*nstr
      nocprt = nint( prtday*secday/dto )*nstr
*     Control intervals for time-averaging of forcings
      avgatm = dtavat.gt.0.0d0 .and. oceanonly.eq.0
      if ( avgatm ) then
        write(*,204) '  Atmos. averaging int. (day) = ',dtavat
        ntavat = nint( dtavat*secday/dto )*nstr
        write(*,201) '  No. of steps in avg. inter. = ',ntavat
        if ( mod(ntsrun,ntavat).ne.0 .or. mod(ntavat,2).ne.0 ) then
          print *,' Unsuitable choice of dtavat; program stops'
          stop
        endif
        nmidat = ntavat/2
      endif
      avgocn = dtavoc.gt.0.0d0 .and. atmosonly.eq.0
      if ( avgocn ) then
        write(*,204) '  Ocean  averaging int. (day) = ',dtavoc
        ntavoc = nint( dtavoc*secday/dto )*nstr
        write(*,201) '  No. of steps in avg. inter. = ',ntavoc
        if ( mod(ntsrun,ntavoc).ne.0 .or. mod(ntavoc,2).ne.0 ) then
          print *,' Unsuitable choice of dtavoc; program stops'
          stop
        endif
        nmidoc = ntavoc/2
      endif
#ifdef get_covar
*     Control intervals for taking covariances
      ntcovat = nint( dtcovat*secday/dto )*nstr
      ntcovoc = nint( dtcovoc*secday/dto )*nstr
#endif /* get_covar */

!~ *     Write Matlab-readable copy of input and derived parameters
!~ *     ----------------------------------------------------------
!~       open (10, file=outdir(1:lenod)//'/input_parameters.m',
!~      &      status='unknown')
!~       INCLUDE './out_param.f'
!~       close(10)

*     Compute initial mass and momentum integrals for constraints
*     -----------------------------------------------------------
      call constr

*     Compute initial potential vorticities from pressures
*     ====================================================
#ifndef atmos_only
*     Ocean
*     -----
*     Internal points
      call qcomp (qo, po, amatoc, yporel, dxom2,
     &            nxpo, nypo, nlo, ddynoc, nlo)
      call qcomp (qom,pom,amatoc, yporel, dxom2,
     &            nxpo, nypo, nlo, ddynoc, nlo)
*     Zonal boundaries + meridional if not cyclic (mixed condition)
      call ocqbdy (qo, po )
      call ocqbdy (qom,pom)
#  ifdef cyclic_ocean
*     Meridional boundaries (periodic)
      call merqcy (qo, po,  amatoc, yporel, dxom2,
     &             nxpo, nypo, nlo, ddynoc, nlo)
      call merqcy (qom,pom, amatoc, yporel, dxom2,
     &             nxpo, nypo, nlo, ddynoc, nlo)
#  endif
#endif /* atmos_only */
#ifndef ocean_only
*     Atmosphere
*     ----------
*     Internal points
      call qcomp (qa, pa, amatat, yparel, dxam2,
     &            nxpa, nypa, nla, ddynat, 1)
      call qcomp (qam,pam,amatat, yparel, dxam2,
     &            nxpa, nypa, nla, ddynat, 1)
*     Zonal boundaries (mixed condition)
      call atqzbd (qa, pa )
      call atqzbd (qam,pam)
*     Meridional boundaries (periodic)
      call merqcy (qa, pa,  amatat, yparel, dxam2,
     &             nxpa, nypa, nla, ddynat, 1)
      call merqcy (qam,pam, amatat, yparel, dxam2,
     &             nxpa, nypa, nla, ddynat, 1)
#endif

  201 format(a,9i13)
  203 format(a,9f13.3)
  204 format(a,9f13.4)
  205 format(a,9f13.5)
  206 format(a,9f13.6)
  207 format(a,9f13.7)
  213 format(a,1p,9d13.3)
* 214 format(a,1p,9d13.4)
  215 format(a,1p,9d13.5)
* 225 format(a,i2,a,9f13.5)
* 226 format(a,i2,a,9f13.6)

      end subroutine init_grids1

      subroutine read_forcings
      
!! Read in mean ocean temperature if we're in atmosphere only mode
#ifdef atmos_only
#  ifdef use_netcdf
      print *,' Mean ocean state for atmos_only case read from netCDF'
      status = nf_open('avges.nc', NF_NOWRITE, tempid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_inq_varid(tempid, 'sst', varid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_get_var_double(tempid, varid, sst)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_get_var_double(tempid, varid, sstm)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_close(tempid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#  else
      print *,' Mean forcing for atmos_only case read from binary'
      inquire (file=atforav, exist=fexist)
      if ( fexist ) then
        open (50, file=atforav, form='unformatted', status='old')
        read (50) sst
        do j=1,nyto
          do i=1,nxto
            sstm(i,j) = sst(i,j)
          enddo
        enddo
        close(50)
        print *,' Atmos mean forcing filename = ',atforav
       else
        print *,' Mean atmos forcing unformatted binary file not found'
        print *,' Filename = ',atforav
        print *,' Program terminates'
        stop
      endif
#  endif
#endif /* atmos_only */

!! Read in mean forcings if we're in ocean only mode
*     Need to supply fnetoc, tauxo and tauyo
#ifdef ocean_only
#  ifdef use_netcdf
      print *,' Mean forcing for ocean_only case read from netCDF'
      status = nf_open('avges.nc', NF_NOWRITE, tempid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_inq_varid(tempid, 'fnetoc', varid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_get_var_double(tempid, varid, fnetoc)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_inq_varid(tempid, 'tauxo', varid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_get_var_double(tempid, varid, tauxo)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_inq_varid(tempid, 'tauyo', varid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_get_var_double(tempid, varid, tauyo)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_close(tempid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#  else
      print *,' Mean forcing for ocean_only case read from binary'
      inquire (file=ocforav, exist=fexist)
      if ( fexist ) then
        open (50, file=ocforav, form='unformatted', status='old')
        read (50) fnetoc
        read (50) tauxo
        read (50) tauyo
        close(50)
        print *,' Ocean mean forcing filename = ',ocforav
       else
        print *,' Mean ocean forcing unformatted binary file not found'
        print *,' Filename = ',ocforav
        print *,' Program terminates'
        stop
      endif
#  endif
#endif
      end subroutine read_forcings

      subroutine init_grids2

!! Next section is for the xcexp experiments
!! There are a few options which require adjustments to
!! the comments to make them current

!! Option 1: use astbar as determined by radiative equilibrium:
      do j=1,nyta
        do i=1,nxta
           xc1ast(i,j) = ( 1.0d0 - xcexp )*astbar(j)
        enddo
      enddo

!! Options 2 & 3: both require reading in a file from the
!! local directory which holds averaged ast information:
!      print *,' Reading in mean AST for xcexp experiment'
!      status = nf_open('avges.nc', NF_NOWRITE, tempid)
!      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
!      status = nf_inq_varid(tempid, 'ast', varid)
!      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
!      status = nf_get_var_double(tempid, varid, xc1ast)
!      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)

!! Option 2: Directly use the mean AST fields:
!      do j=1,nyta
!        do i=1,nxta
!           xc1ast(i,j) = ( 1.0d0 - xcexp )*xc1ast(i,j)
!        enddo
!      enddo

!! Option 3: Use the AST from data, but average
!! it in the x-direction, so that there are no 'bumps':
!      do j=1,nyta
!        do i=2,nxta
!           xc1ast(1,j) = xc1ast(1,j) + xc1ast(i,j)
!        enddo
!        xc1ast(1,j) = ( 1.0d0 - xcexp )*xc1ast(1,j)/dble(nxta)
!        do i=2,nxta
!           xc1ast(i,j) = xc1ast(1,j)
!        enddo
!      enddo

*     Initialise forcing arrays
*     =========================
      call xforc

#ifndef atmos_only
*     Ocean forcings
*     --------------
*     Set those quantities which have not already been either
*     set by xforc, or in the ocean only case, read from a file
!$OMP PARALLEL DO DEFAULT (NONE)
!$OMP&         PRIVATE (i,j)
!$OMP&         SHARED  (entoc)
!$OMP&         SCHEDULE (STATIC)
      do j=1,nypo
        do i=1,nxpo
          entoc(i,j) = 0.0d0
        enddo
      enddo
!$OMP END PARALLEL DO
      do k=1,nlo-1
        xon(k) = 0.0d0
      enddo
#  ifdef cyclic_ocean
      do k=1,nlo-1
        enisoc(k) = 0.0d0
        eninoc(k) = 0.0d0
      enddo
*     Also zero ocean boundary Jacobian arrays
*     These will be correctly recomputed in the cyclic
*     case, and are not used in the finite box case
      do k=1,nlo
        ajisoc(k) = 0.0d0
        ajinoc(k) = 0.0d0
      enddo
#  endif
#endif /* atmos_only */

#ifndef ocean_only
*     Atmosphere forcings
*     -------------------
      do j=1,nypa
        do i=1,nxpa
          entat(i,j) = 0.0d0
        enddo
      enddo
      do k=1,nla-1
        enisat(k) = 0.0d0
        eninat(k) = 0.0d0
        xan(k) = 0.0d0
      enddo
#endif

*     Compute tridiagonal coefficients for Helmholtz equations
*     ========================================================
*     Compute Del-sqd part once and for all; add def. rad. later.
*     In all cases the tridiagonal system is diagonally dominant and
*     thus numerically well-behaved. The Del-sqd operator to be
*     inverted is the finite-difference form in both x- and y-directions.
*     Ordering of wavenumber components in bd2oc,
*     bd2at is appropriate to FFTPACK.

#ifndef atmos_only
*     Oceanic
*     -------
      aoc = 1.0d0/( dyo*dyo )
#  ifdef cyclic_ocean
*     Zonally cyclic ocean
      do i=2,nxto/2
        i1 = 2*i - 1
        bd2oc(i1-1) = -2.0d0*aoc
     &               + 2.0d0*dxom2*( cos( (i-1)*TWOPI/nxto ) - 1.0d0 )
        bd2oc( i1 ) = bd2oc(i1-1)
      enddo
*     When using FFTPACK, special values go to 1 and nxto
      bd2oc(  1 ) = -2.0d0*aoc
      bd2oc(nxto) = -2.0d0*aoc - 4.0d0*dxom2
*     Initialise oceanic FFT coefficient array
      call drffti (nxto, oftwrk)
#  else
*     Finite box ocean
      do i=2,nxto
        bd2oc(i-1) = -2.0d0*aoc
     &              + 2.0d0*dxom2*( cos( (i-1)*PI/nxto ) - 1.0d0 )
      enddo
      bd2oc(nxto) = 0.0d0
*     Initialise oceanic FFT coefficient array
      call dsinti (nxto-1, oftwrk)
#  endif
#endif /* atmos_only */

#ifndef ocean_only
*     Atmospheric
*     -----------
      aat = 1.0d0/( dya*dya )
      do i=2,nxta/2
        i1 = 2*i - 1
        bd2at(i1-1) = -2.0d0*aat
     &               + 2.0d0*dxam2*( cos( (i-1)*TWOPI/nxta ) - 1.0d0 )
        bd2at( i1 ) = bd2at(i1-1)
      enddo
*     When using FFTPACK, special values go to 1 and nxta
      bd2at(  1 ) = -2.0d0*aat
      bd2at(nxta) = -2.0d0*aat - 4.0d0*dxam2
*     Initialise atmospheric FFT coefficient array
      call drffti (nxta, aftwrk)
#endif

*     Compute homogeneous solutions for ocean and atmosphere
      call homsol

      noutstepoc = ntsrun/noutoc + 1
      noutstepat = ntsrun/noutat + 1
      numoutsteps = ntsrun/nocmon + 1
      nt = nsteps0
      ntdone = 0
      tday = tini*daysyr
      tyrs = tini
      call monnc_comp

*     Check validity of initial state, before
*     proceeding with opening of output files
*     ---------------------------------------
      solnok = .true.
      call valids (solnok)
      if ( .not.solnok ) then
        print *,' *** valids has detected invalid values ***'
        write(*,'(a,i12,f12.2,f11.4)')
     &        '  problem occurs at nt, tday, tyrs = ',nt,tday,tyrs
        print *,' program will terminate after diagnostic printout'
        call prsamp (nt)
        call cfltry (nt)
        print *,' '
        print *,' program terminates'
        stop
      endif

#ifdef use_netcdf
*     netCDF output initialisation
*     ============================
*     N.B. if the error message "One or more variable sizes
*     violate format constraints" occurs when initialising
*     a netCDF file, try either coarsening spatial
*     and/or temporal sampling, or create that file with
*     "64-bit offset" format; see netCDF manuals for details
      print *,' '

#  ifndef atmos_only
*     Work out the size (in Gbytes) of a 4-D oceanic p-grid variable
*     (believed to be the maximal size variable in the netCDF dumps)
*     This assumes single precision (4 byte) dumps; if using
*     double precision the initial 4 should be changed to 8.
      poc4dsz = 4.0d0*(1 + nxto/nsko)*(1 + nyto/nsko)*nlo*noutstepoc
     &         /dble(1024**3)
      write(*,'(a,a,f7.3)') '  netCDF3: 4-dim oceanic p-grid',
     &                      ' record size (Gbytes) = ',poc4dsz
      if ( poc4dsz.ge.4.0d0 ) then
        print *,' netCDF3 ERROR: netCDF3 limit = 4.0 Gbytes'
        print *,' netCDF3 ERROR: program terminates'
        STOP
      endif
*     Ocean output files
*     ------------------
      status = nf_create (outdir(1:lenod)//'/ocpo.nc',
***  &                    NF_CLOBBER, ocpid)
*     Alternative creation call with "Large File Support"
     &                    IOR(NF_CLOBBER,NF_64BIT_OFFSET), ocpid)
      if ( status.ne.NF_NOERR ) then
        print *,' problem creating file ocpo.nc'
        call handle_err (status, subnam)
       else
        print *,' ocpo.nc file created'
      endif
      status = nf_create (outdir(1:lenod)//'/ocsst.nc',
     &                    NF_CLOBBER, octid)
*     Alternative creation call with "Large File Support"
***  &                    IOR(NF_CLOBBER,NF_64BIT_OFFSET), octid)
      if ( status.ne.NF_NOERR ) then
        print *,' problem creating file ocsst.nc'
        call handle_err (status, subnam)
       else
        print *,' ocsst.nc file created'
      endif
      call ocnc_init (nsko, outfloc)
      print *,' Ocean netCDF files initialised'
      call ocnc_out (nsko, outfloc)
#    ifdef qoc_diag
*     Optional ocean vorticity diagnostics
      status = nf_create (outdir(1:lenod)//'/qocdiag.nc',
***  &                    NF_CLOBBER, qocncid)
*     Alternative creation call with "Large File Support"
     &                    IOR(NF_CLOBBER,NF_64BIT_OFFSET), qocncid)
      if ( status.ne.NF_NOERR ) then
        print *,' problem creating file qocdiag.nc'
        call handle_err (status, subnam)
       else
        print *,' qocdiag.nc file created'
      endif
      call qocdiag_init (nsko)
      print *,' qocdiag netCDF file initialised'
#    endif
#  endif /* atmos_only */

#  ifndef ocean_only
*     Work out the size (in Gbytes) of a 4-D atmospheric p-grid variable
*     (believed to be the maximal size variable in the netCDF dumps)
*     This assumes single precision (4 byte) dumps; if using
*     double precision the initial 4 should be changed to 8.
      pat4dsz = 4.0d0*(1 + nxta/nska)*(1 + nyta/nska)*nla*noutstepat
     &         /dble(1024**3)
      write(*,'(a,a,f7.3)') '  netCDF3: 4-dim atmospheric p-grid',
     &                      ' record size (Gbytes) = ',pat4dsz
      if ( pat4dsz.ge.4.0d0 ) then
        print *,' netCDF3 ERROR: netCDF3 limit = 4.0 Gbytes'
        print *,' netCDF3 ERROR: program terminates'
        STOP
      endif
*     Atmosphere output files
*     -----------------------
      status = nf_create (outdir(1:lenod)//'/atpa.nc',
     &                    NF_CLOBBER, atpid)
*     Alternative creation call with "Large File Support"
***  &                    IOR(NF_CLOBBER,NF_64BIT_OFFSET), atpid)
      if ( status.ne.NF_NOERR ) then
        print *,' problem creating file atpa.nc'
        call handle_err (status, subnam)
       else
        print *,' atpa.nc file created'
      endif
      status = nf_create (outdir(1:lenod)//'/atast.nc',
     &                    NF_CLOBBER, attid)
*     Alternative creation call with "Large File Support"
***  &                    IOR(NF_CLOBBER,NF_64BIT_OFFSET), attid)
      if ( status.ne.NF_NOERR ) then
        print *,' problem creating file atast.nc'
        call handle_err (status, subnam)
       else
        print *,' atast.nc file created'
      endif
      call atnc_init (nska, outflat)
      print *,' Atmos. netCDF files initialised'
      call atnc_out (nska, outflat)
#  endif

*     Additional monitoring files
*     ---------------------------
      status = nf_create (outdir(1:lenod)//'/monit.nc',
     &                    NF_CLOBBER, monncid)
      if ( status.ne.NF_NOERR ) then
        print *,' problem creating file monit.nc'
        call handle_err (status, subnam)
       else
        print *,' monit.nc file created'
      endif
*     Ensure that monitoring and area-averaging are done
*     synchronously as post-processing routines assume this.
      call monnc_init
      print *,' Diagnostics netCDF file initialised'
      call monnc_out
#  ifdef get_areav
*     Area average output file
*     ------------------------
      status = nf_create (outdir(1:lenod)//'/areas.nc',
     &                    NF_CLOBBER, arencid)
*     Alternative creation call with "Large File Support"
***  &                    IOR(NF_CLOBBER,NF_64BIT_OFFSET), arencid)
      if ( status.ne.NF_NOERR ) then
        print *,' problem creating file areas.nc'
        call handle_err (status, subnam)
       else
        print *,' areas.nc file created'
      endif
#  endif
#endif

*     Print sample of initial state
*     -----------------------------
*     (validity of initial state already tested)
      call prsamp (nt)
#ifdef get_areav
      call areavg
#endif

*     Time averaging of forcing fields
      if ( avgatm .or. avgocn ) call tavini

#ifdef get_covar
      call covini
#endif /* get_covar */

  201 format(a,9i13)
  203 format(a,9f13.3)
  204 format(a,9f13.4)
  205 format(a,9f13.5)
  206 format(a,9f13.6)
  207 format(a,9f13.7)
  213 format(a,1p,9d13.3)
* 214 format(a,1p,9d13.4)
  215 format(a,1p,9d13.5)
* 225 format(a,i2,a,9f13.5)
* 226 format(a,i2,a,9f13.6)

      end subroutine init_grids2


      subroutine mainloop(nsteps0, nsteps)
      double precision :: told
      integer nsteps0, nsteps

****************    START OF MAIN TIMESTEPPING LOOP     ****************

      do 1000 nt=nsteps0+1,nsteps

        if ( mod(nt,nstr).eq.1 ) then

*         Compute forcing fields, if on ocean step
#ifndef ocean_only
          call xforc
#endif

#ifndef atmos_only

* a consistent calculation independent of averaging, filtering 
* of dpo_dt, dsst_dt etc can be made here by saving po in dpo_dt 
* and differencing (po-dpo_dt)/tdto after ocqbdy, currently removed

*         Step ocean mixed layer
          call oml
#    ifdef qoc_diag
*         Compute dq/dt diagnostics here, after entoc has
*         been computed, but before po and qo are updated
          if ( mod(ntdone,noutoc).eq.0 ) then
            call qocdiag_out (nsko)
          endif
#    endif

*         Step ocean qg model
          call qgostep

*         step omuse timers
          told=tos
          tos=tosm+2*dto
          tosm=told
          tocean=tocean+(tos-tosm)
          if(tosm.GT.2*dto) then
            tosm=tosm-2*dto
            tos=tos-2*dto
          endif
          
*         Invert pv into oceanic pressures
          call ocinvq

*         Compute pv on ocean boundaries (mixed condition)
          call ocqbdy (qo, po)

#endif /* not atmos_only */

        endif

#ifndef ocean_only
*       Step atmospheric mixed layer
        call aml

*       Step atmospheric qg channel model
        call qgastep

        told=tas
        tas=tasm+2*dta
        tasm=told
        tatm=tatm+(tas-tasm)
        if(tasm.GT.2*dta) then
          tasm=tasm-2*dta
          tas=tas-2*dta
        endif

*       Invert pv into atmospheric pressures
        call atinvq

*       Compute pv on atmosphere zonal boundary (mixed condition)
        call atqzbd (qa, pa)
#endif

*       Timestep done; do checking and diagnostics as necessary
        ntdone=ntdone+1
!        ntdone = nt - nsteps0
        tday = nt*dta/secday
        tyrs = nt*dta/secsyr

*       Periodically check validity of solution
*       ---------------------------------------
        if ( mod(ntdone,nvalid).eq.0 ) then
          call valids (solnok)
          if ( .not.solnok ) then
            print *,' *** valids has detected invalid values ***'
            write(*,'(a,i12,f12.2,f11.4)')
     &            '  problem occurs at nt, tday, tyrs = ',nt,tday,tyrs
            print *,' program will terminate after diagnostic printout'
            call monnc_comp
            call prsamp (nt)
            call cfltry (nt)
#ifdef use_netcdf
            print *,' make netCDF dumps of fields before exit'

            call monnc_out
            status = nf_close(monncid)
            if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#  ifdef get_areav
            status = nf_close(arencid)
            if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#  endif

#  ifndef atmos_only
            call ocnc_out (nsko, outfloc)
            status = nf_close(ocpid)
            if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
            status = nf_close(octid)
            if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#    ifdef qoc_diag
            call qocdiag_out (nsko)
            status = nf_close(qocncid)
            if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#    endif
#  endif /* atmos_only */

#  ifndef ocean_only
            call atnc_out (nska, outflat)
            status = nf_close(atpid)
            if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
            status = nf_close(attid)
            if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#  endif
#endif
            print *,' program terminates'
            stop
          endif
        endif

*       Suppress computational mode of leapfrog scheme
*       ==============================================
#ifndef atmos_only
        if ( (oalpha.EQ.0.).AND.mod(nt-1,25*nstr).eq.0 ) then
*         Average oceanic time levels
*         ---------------------------
          tocean=tocean-(tos-tosm)/2
          tos=(tos+tosm)/2
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k)
!$OMP&         SHARED  (qo,qom,po,pom,sst,sstm)

          do k=1,nlo
!$OMP       DO SCHEDULE (STATIC)
            do j=1,nypo
              do i=1,nxpo
                qo(i,j,k) = 0.5d0*( qo(i,j,k)+qom(i,j,k) )
                po(i,j,k) = 0.5d0*( po(i,j,k)+pom(i,j,k) )
              enddo
            enddo
!$OMP       END DO
          enddo

!$OMP     DO SCHEDULE (STATIC)
          do j=1,nyto
            do i=1,nxto
              sst(i,j) = 0.5d0*( sst(i,j)+sstm(i,j) )
            enddo
          enddo
!$OMP     END DO NOWAIT

!$OMP     END PARALLEL

*         Also average constraint variables
          do k=1,nlo-1
            dpioc(k) = 0.5d0*( dpioc(k) + dpiocp(k) )
          enddo
#  ifdef cyclic_ocean
          do k=1,nlo
            ocncs(k) = 0.5d0*( ocncs(k) + ocncsp(k) )
            ocncn(k) = 0.5d0*( ocncn(k) + ocncnp(k) )
          enddo
#  endif
        endif
#endif /* atmos_only */

#ifndef ocean_only
        if ( mod(nt-1,100).eq.0 ) then
*         Average atmospheric time levels
*         -------------------------------
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k)
!$OMP&         SHARED  (qa,qam,pa,pam,ast,astm,hmixa,hmixam)

          do k=1,nla
!$OMP       DO SCHEDULE (STATIC)
            do j=1,nypa
              do i=1,nxpa
                qa(i,j,k) = 0.5d0*( qa(i,j,k)+qam(i,j,k) )
                pa(i,j,k) = 0.5d0*( pa(i,j,k)+pam(i,j,k) )
              enddo
            enddo
!$OMP       END DO
          enddo

!$OMP     DO SCHEDULE (STATIC)
          do j=1,nyta
            do i=1,nxta
              ast(i,j) = 0.5d0*( ast(i,j)+astm(i,j) )
              hmixa(i,j) = 0.5d0*( hmixa(i,j)+hmixam(i,j) )
            enddo
          enddo
!$OMP     END DO NOWAIT

!$OMP     END PARALLEL

*         Also average constraint variables
          do k=1,nla-1
            dpiat(k) = 0.5d0*( dpiat(k) + dpiatp(k) )
          enddo
          do k=1,nla
            atmcs(k) = 0.5d0*( atmcs(k) + atmcsp(k) )
            atmcn(k) = 0.5d0*( atmcn(k) + atmcnp(k) )
          enddo
        endif
#endif

*       Occasionally dump restart file
*       ------------------------------
        if ( resdmp ) then
          if ( mod(ntdone,noutre).eq.0 ) then
*           First check the solution is valid
            call valids (solnok)
            if ( solnok ) then
              print *,' '
              print *,' Writing restart dump at nt, tyrs = ',nt,tyrs
*             Open and close restart file each time because
*             flush and rewind not yet working properly under Linux.
#ifdef use_netcdf
              status = nf_create (outdir(1:lenod)//'/restart.nc',
     &                            NF_CLOBBER, restid)
              if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
              call resave_nc
              status = nf_close(restid)
              if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#else
              open (iorest, file=outdir(1:lenod)//'/restart',
     &              form='unformatted', status='unknown')
              call resave (iorest)
              close (iorest)
#endif
            endif
          endif
        endif

*       Monitoring of current state if required
*       ---------------------------------------
*       Ensure that monitoring and area-averaging are done
*       synchronously as post-processing routines assume this.
        if ( mod(ntdone,nocmon).eq.0 ) then
          call monnc_comp
#ifdef use_netcdf
          call monnc_out
#endif
#ifdef get_areav
          call areavg
#endif
        endif
        if ( mod(ntdone,nocprt).eq.0 ) then
          call prsamp (nt)
        endif

#ifdef use_netcdf
*       Write full fields to netCDF files
*       ---------------------------------
#  ifndef atmos_only
        if ( mod(ntdone,noutoc).eq.0 ) then
          call ocnc_out (nsko, outfloc)
        endif
#  endif /* atmos_only */
#  ifndef ocean_only
        if ( mod(ntdone,noutat).eq.0 ) then
          call atnc_out (nska, outflat)
        endif
#  endif
#endif

*       Add contributions to running means of forcing fields
*       ----------------------------------------------------
        if ( avgatm ) then
          if ( mod(ntdone,ntavat).eq.nmidat ) call tavatm
        endif
        if ( avgocn ) then
          if ( mod(ntdone,ntavoc).eq.nmidoc ) call tavocn
        endif

#ifdef get_covar
*       Accumulate contributions to covariance matrices
*       -----------------------------------------------
        if ( mod(ntdone,ntcovat).eq.0 ) call covatm
        if ( mod(ntdone,ntcovoc).eq.0 ) call covocn
#endif /* get_covar */

 1000 continue

****************     END OF MAIN TIMESTEPPING LOOP      ****************

      end subroutine mainloop

      subroutine finalize
      
#ifdef use_netcdf
*     Close monitoring, netCDF and optional restart files
      status = nf_close(monncid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#  ifdef get_areav
      status = nf_close(arencid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#  endif

#  ifndef atmos_only
      status = nf_close(ocpid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_close(octid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#    ifdef qoc_diag
*     Compute qoc diagnostics on final state. Ideally oml would also be
*     called here to update entoc, but this would wrongly change sst and sstm
      call qocdiag_out (nsko)
      status = nf_close(qocncid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#    endif
#  endif /* atmos_only */

#  ifndef ocean_only
      status = nf_close(atpid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      status = nf_close(attid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#  endif
#endif

*     dump restart file (either netCDF or unformatted)
#ifdef use_netcdf
      status = nf_create (outdir(1:lenod)//'/lastday.nc',
     &                    NF_CLOBBER, restid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
      call resave_nc
      status = nf_close(restid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#else
      open (14, file=outdir(1:lenod)//'/last.day', form='unformatted',
     &      status='unknown')
      call resave (14)
      close(14)
#endif

*     Compute and output time-averaged fields
#ifdef use_netcdf
      status = nf_create (outdir(1:lenod)//'/avges.nc',
     &                    NF_CLOBBER, tavncid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#endif
      call tavout
#ifdef use_netcdf
      status = nf_close(tavncid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#endif

#ifdef get_covar
*     Write out covariance matrix
#  ifdef use_netcdf
      status = nf_create (outdir(1:lenod)//'/covar.nc',
     &                    NF_CLOBBER, covncid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#  endif
      call covout
#  ifdef use_netcdf
      status = nf_close(covncid)
      if ( status.ne.NF_NOERR ) call handle_err (status, subnam)
#  endif
#endif /* get_covar */

*     Precautionary call to prsamp, if not already done
      if ( mod(ntdone,nocprt).ne.0 ) then
        call prsamp (nt)
      endif

*     Compute final mass and momentum integrals for constraints
*     ---------------------------------------------------------
      call constr

      print *,' '
      write(*,'(a,i12,f11.2,f11.4)')
     &        '  End of run at nt, tday, tyrs = ',nt,tday,tyrs


      end subroutine finalize

      function print_diag() result(ret)
        integer :: ret
        double precision :: qt2dif(3),qt4dif(3),qotjac(3),qotent(3),dqdt(3)
      
        print*, "time:", tday
        print*, "d2 diss", sum(ah2doc)
        print*, "d4 diss", sum(ah4doc)
        print*, "bottom drag", btdgoc
        print*, "wind", utauoc
        print*, "pe_t", sum(ddtpeoc)
        print*, "ke_t", sum(ddtkeoc)
        print*, "ke", sum(kealoc)
        print*, "pe", sum(0.5*rhooc*gpoc*et2moc)
        print*, "hfmloc ", hfmloc
        print*, "hfluxes ", slhfav, oradav
        print*, "hfluxes2 ", arocav, arlaav
        print*, "entrainment ",  pkenoc, centoc
        print*, "dpioc:", dpioc
        print*, "dpiocp:", dpiocp
#ifdef cyclic_ocean 
        print*, "ocncs:", ocncs
        print*, "ocncsp:", ocncsp
        print*, "ocncn:", ocncn
        print*, "ocncnp:", ocncnp
#endif
        print*, "emfroc:", emfroc
        print*, "ermaso:", ermaso
        
        call qt4dif_int(qt2dif,qt4dif,qotjac,qotent,dqdt)
        qt2dif=rhooc*hoc*qt2dif*ocnorm/fnot
        print*, "qt2dif", sum(qt2dif)
        qt4dif=rhooc*hoc*qt4dif*ocnorm/fnot
        print*, "qt4dif", sum(qt4dif)
        qotjac=rhooc*hoc*qotjac*ocnorm/fnot
        print*, "qotjac", sum(qotjac)
        qotent=rhooc*hoc*qotent*ocnorm/fnot
        print*, "qotent", sum(qotent), utauoc-btdgoc
        dqdt=rhooc*hoc*dqdt*ocnorm/fnot
        print*, "dqdt", sum(dqdt)
        print*, "sum:", sum(ddtkeoc)-(sum(ddtpeoc)+utauoc-btdgoc-sum(ah4doc))
        ret=0
      end function
      
      end module

! include diagnostic routines
#include "qocdiag.F"

