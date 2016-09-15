c subroutines from q-gcm.F

#include "omuse_valsubs.F"

c***********************************************************************
c
      SUBROUTINE ipbget (buffer, iounit)
*
*     Reads records from unit iounit until a valid one (i.e.
*     one not marked with a "comment" character) is found,
*     then returns this valid character string for processing.
*     The comment marker in Q-GCM is deemed to be
*     an exclamation mark "!" in the first column.
*
      IMPLICIT NONE
*
*     Subroutine arguments
      character (len=80) :: buffer
      integer, INTENT(IN) :: iounit
*
*     Local variables

  100 continue
      read (iounit, err=200, fmt='(a80)') buffer
      if ( buffer(1:1).eq.'!' ) goto 100
      return

  200 continue
      print *,' Error reading character buffer from iounit = ',iounit
      print *,' Program terminates in ipbget'
      stop

      END SUBROUTINE ipbget
c
c***********************************************************************
c
      SUBROUTINE zeroin
*
*     Set initial state to zero pressure, and radiative
*     equilibrium with unperturbed (background) forcing,
*     i.e. mixed layer temperature anomalies are all zero.

*     Modules
      USE parameters
#ifndef ocean_only
      USE atstate
#endif
#ifndef atmos_only
      USE ocstate
#endif
      USE intrfac, ONLY : sst, sstm, ast, astm, hmixa, hmixam, hmat

      IMPLICIT NONE

*     Subroutine arguments
*
*     Local variables
      integer i,j,k

*     Initialise atmosphere fields (incl. mixed layer)
*     ------------------------------------------------
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k)
!$OMP&         SHARED  (ast,astm,hmixa,hmixam,hmat)
#ifndef ocean_only
!$OMP&         SHARED  (pa,pam)

*     Initialise atmospheric pressure
      do k=1,nla
!$OMP   DO SCHEDULE (STATIC)
        do j=1,nypa
          do i=1,nxpa
            pa(i,j,k) = 0.0d0
            pam(i,j,k) = 0.0d0
          enddo
        enddo
!$OMP   END DO NOWAIT
      enddo
#endif

*     Initialise atmos. mixed layer rel. temp. and thickness
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyta
        do i=1,nxta
          ast(i,j) = 0.0d0
          astm(i,j) = 0.0d0
          hmixa(i,j) = hmat
          hmixam(i,j) = hmat
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

*     Initialise ocean fields (incl. mixed layer)
*     -------------------------------------------
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k)
!$OMP&         SHARED  (sst,sstm)
#ifndef atmos_only
!$OMP&         SHARED  (po,pom)

*     Initialise ocean pressure
      do k=1,nlo
!$OMP   DO SCHEDULE (STATIC)
        do j=1,nypo
          do i=1,nxpo
            po(i,j,k) = 0.0d0
            pom(i,j,k) = 0.0d0
          enddo
        enddo
!$OMP   END DO NOWAIT
      enddo
#endif /* atmos_only */

*     Initialise oceanic mixed layer rel. temp.
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyto
        do i=1,nxto
          sst(i,j) = 0.0d0
          sstm(i,j) = 0.0d0
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL
      print *,' '
      write (*,*) ' Zero initialisation complete'

      END SUBROUTINE zeroin
c
c***********************************************************************
c
      SUBROUTINE rbalin
*
*     Initialise in radiative equilibrium, with
*     radiative forcing perturbation fsprim included.

*     Modules
      USE parameters
      USE atconst
      USE occonst
#ifndef ocean_only
      USE atstate
#endif /* not ocean_only */
#ifndef atmos_only
      USE ocstate
#endif /* not atmos_only */
      USE intrfac, ONLY : sst, sstm, sstbar, ast, astm, astbar,
     &                    hmixa, hmixam, hmat
      USE xfosubs, ONLY : fsprim
      USE radiate

      IMPLICIT NONE

*     Subroutine arguments
*
*     Local parameters
      logical prtval
      parameter ( prtval = .false. )
*
*     Local variables
      integer i,j,k,nxco,nyco,nxca,nyca
      double precision plfac(nla),play
**    double precision sumc,sumcp,p1off

*     Centres of domains if needed for sample printing
      nxco = (nxpo+1)/2
      nyco = (nypo+1)/2
      nxca = (nxpa+1)/2
      nyca = (nypa+1)/2

      print *,' '
      write (*,*) ' Radiative balance initialisation:'
      write (*,*) ' ---------------------------------'

#ifndef ocean_only
*     Derive suitable multiplier of Fs' for each atmos. layer,
*     from the eta coefficients of Fs' derived in radiat.
*     We have nla layers but only nla-1 eta coeffts,
*     and so need an extra constraint.
      plfac(1) = 0.0d0
      do k=2,nla
        plfac(k) = plfac(k-1) - gpat(k-1)*rbetat(k-1)
      enddo

*     Option 1: leave alone. This gives p(1) = 0 everywhere;
*     no pressure gradient and thus no flow in layer 1.

*     Option 2: apply offset so that barotropic p = 0, as in zeroin.
**    sumc = 0.0d0
**    sumcp = 0.0d0
**    do k=1,nla
**      sumc = sumc + ctl2mat(k,1)
**      sumcp = sumcp + ctl2mat(k,1)*plfac(k)
**    enddo
**    p1off = -sumcp/sumc
**    sumcp = 0.0d0
**    do k=1,nla
**      plfac(k) = plfac(k) + p1off
**      sumcp = sumcp + ctl2mat(k,1)*plfac(k)
**    enddo
**    print *,' Barotropic coefft. in rbalin = ',sumcp
      write(*,206) '  Layer coeffts for pa, plfac = ',
     &             (plfac(k),k=1,nla)
  206 format(a,9f13.6)
#endif

      if ( prtval ) then
        print *,' '
        write (*,*) ' Initial relative ast:'
        do j=nyta,1,-1
          if ( astbar(j).gt.tat(1) ) then
            write (*,'(i6,f16.8,a)') j,astbar(j),'  convect'
           else
            write (*,'(i6,f16.8,a)') j,astbar(j)
          endif
        enddo
      endif

      if ( prtval ) then
        print *,' '
        write (*,*) ' Initial relative sst:'
        do j=nyto,1,-1
          if ( sstbar(j).lt.toc(1) ) then
            write (*,'(i6,f16.8,a)') j,sstbar(j),'  convect'
           else
            write (*,'(i6,f16.8,a)') j,sstbar(j)
          endif
        enddo
      endif

*     Initialise atmosphere fields (incl. mixed layer)
*     ------------------------------------------------
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k,play)
!$OMP&         SHARED  (ast,astm,astbar,hmixa,hmixam,hmat)
#ifndef ocean_only
!$OMP&         SHARED  (pa,pam,plfac,yparel)

*     Initialise atmospheric pressure
      do k=1,nla
!$OMP   DO SCHEDULE (STATIC)
        do j=1,nypa
          play = plfac(k)*fsprim( yparel(j) )
          do i=1,nxpa
            pa(i,j,k) = play
            pam(i,j,k) = play
          enddo
        enddo
!$OMP   END DO NOWAIT
      enddo
#endif /* not ocean_only */

*     Initialise atmos. mixed layer rel. temp. and thickness
*     The a.m.l is initialised even in the ocean_only
*     case, to provide an upper boundary condition
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyta
        do i=1,nxta
          ast(i,j) = astbar(j)
          astm(i,j) = astbar(j)
          hmixa(i,j) = hmat
          hmixam(i,j) = hmat
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

#ifndef ocean_only
      if ( prtval ) then
        print *,' '
*       Print meridional sections at the centre of the domain,
*       since the initial values are a function of y only
        write(*,*) ' Meridional sections of initial',
     &             ' values across atmos. centre'
        do k=1,nla
          print *,' '
          write(*,'(a,i2,a)') '  Layer k = ',k,':'
          write(*,'(a,a)') '     j        pa(k)            pam(k)'
          do j=nypa,1,-1
            write(*,240) j,pa(nxca,j,k),pam(nxca,j,k)
          enddo
        enddo
      endif
#endif /* not ocean_only */

*     Initialise ocean fields (incl. mixed layer)
*     -------------------------------------------
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k,play)
!$OMP&         SHARED  (sst,sstm,sstbar)
#ifndef atmos_only
!$OMP&         SHARED  (po,pom)

*     Initialise ocean pressure
      do k=1,nlo
!$OMP   DO SCHEDULE (STATIC)
        do j=1,nypo
          do i=1,nxpo
            po(i,j,k) = 0.0d0
            pom(i,j,k) = 0.0d0
          enddo
        enddo
!$OMP   END DO NOWAIT
      enddo
#endif /* not atmos_only */

*     Initialise oceanic mixed layer rel. temp.
*     The o.m.l is initialised even in the atmos_only
*     case, to provide a lower boundary condition
!$OMP DO SCHEDULE (STATIC)
      do j=1,nyto
        do i=1,nxto
          sst(i,j) = sstbar(j)
          sstm(i,j) = sstbar(j)
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

#ifndef atmos_only
      if ( prtval ) then
        print *,' '
*       Print meridional sections at the centre of the domain,
*       since the initial values are a function of y only
        write(*,*) ' Meridional sections of initial',
     &             ' values across atmos. centre'
        do k=1,nlo
          print *,' '
          write(*,'(a,i2,a)') '  Layer k = ',k,':'
          write(*,'(a,a)') '     j        po(k)            pom(k)'
          do j=nypo,1,-1
            write(*,240) j,po(nxco,j,k),pom(nxco,j,k)
          enddo
        enddo
      endif
#endif /* not atmos_only */

      print *,' '
      write (*,*) ' Rbal initialisation complete'

  240 format(i6,1p,5d18.10)

      END SUBROUTINE rbalin
c
c***********************************************************************
c
      SUBROUTINE prsamp (nt)
*
*     Print out a sample of the current state. Useful for
*     tracking progress, or for testing/debugging purposes

*     Modules
      USE parameters, ONLY : nxpo, nypo, nxto, nyto, nlo,
     &                       nxpa, nypa, nxta, nyta, nla, fnot
#ifndef ocean_only
      USE atconst, ONLY : ddynat
      USE atstate, ONLY : pa, qa, wekpa
#endif
#ifndef atmos_only
      USE occonst, ONLY : ddynoc
      USE ocstate, ONLY : po, qo, wekpo
#endif
      USE intrfac, ONLY : sst, ast, hmixa
      USE timinfo, ONLY : tday, tyrs
      USE monitor, ONLY : pavgoc, qavgoc, kealoc, et2moc, osfmin,
     &                    osfmax, occirc, kealat, et2mat

      IMPLICIT NONE

*     Subroutine arguments
      integer, INTENT(IN) :: nt
*
*     Local parameters
      double precision extrem
      parameter ( extrem=1.0d30 )
*
*     Local variables
      integer i,j,k
#ifndef atmos_only
      integer nxco,nyco
      double precision sstmin,sstmax,poref(nlo)
#endif /* not atmos_only */
#ifndef ocean_only
      integer nxca,nyca
      double precision astmin,astmax,hmxmin,hmxmax,paref(nla)
#endif /* not ocean_only */

*     Centres of domains for sample printing
#ifndef atmos_only
      nxco = (nxpo+1)/2
      nyco = (nypo+1)/2
#endif /* not atmos_only */
#ifndef ocean_only
      nxca = (nxpa+1)/2
      nyca = (nypa+1)/2
#endif /* not ocean_only */

*     Define reference pressures on equatorward side of domain
*     --------------------------------------------------------
#ifndef ocean_only
      if ( fnot.gt.0.0d0 ) then
        do k=1,nla
          paref(k) = pa(1,1,k)
        enddo
       else if ( fnot.lt.0.0d0 ) then
        do k=1,nla
          paref(k) = pa(1,nypa,k)
        enddo
      endif
#endif /* not ocean_only */
#ifndef atmos_only
      if ( fnot.gt.0.0d0 ) then
        do k=1,nlo
          poref(k) = po(1,1,k)
        enddo
       else if ( fnot.lt.0.0d0 ) then
        do k=1,nlo
          poref(k) = po(1,nypo,k)
        enddo
      endif
#endif /* not atmos_only */

      print *,' '
      write(*,'(a,i12,f11.2,f11.4)')
     &      '  Sample output at nt, tday, tyrs = ',nt,tday,tyrs

*     Print some ocean spot values and extrema
*     ----------------------------------------
#ifndef atmos_only
      write(*,217) '  po(k) at centre = ',(po(nxco,nyco,k),k=1,nlo)
      write(*,217) '  qo(k) at centre = ',(qo(nxco,nyco,k),k=1,nlo)
**    write(*,217) '  po(1) on S. bdy = ',(po(i,  1 ,1),i=nxco-1,nxco+1)
**    write(*,217) '  qo(1) on S. bdy = ',(qo(i,  1 ,1),i=nxco-1,nxco+1)
**    write(*,217) '  po(1) on N. bdy = ',(po(i,nypo,1),i=nxco-1,nxco+1)
**    write(*,217) '  qo(1) on N. bdy = ',(qo(i,nypo,1),i=nxco-1,nxco+1)
**    write(*,217) '   sst  on S. bdy = ',(sst(i,  1 ),i=nxco-1,nxco+1)
**    write(*,217) '   sst  on N. bdy = ',(sst(i,nyto),i=nxco-1,nxco+1)
      write(*,217) '  po(k)  average  = ',(pavgoc(k),k=1,nlo)
      write(*,217) '  qo(k)  average  = ',(qavgoc(k),k=1,nlo)
**    write(*,217) '  Ocean kealoc(k) = ',(kealoc(k),k=1,nlo)
**    write(*,217) '  Ocean et2moc(k) = ',(et2moc(k),k=1,nlo-1)
**    write(*,217) '  osfmin(k)  (Sv) = ',(osfmin(k),k=1,nlo)
**    write(*,217) '  osfmax(k)  (Sv) = ',(osfmax(k),k=1,nlo)
**    write(*,217) '  occirc(k)  (Sv) = ',(occirc(k),k=1,nlo)
*     Meridional section at centre
**    write(*,*) ' Meridional section at ocean centre'
**    write(*,'(a,a)') '     j         po(1)             qo(1)       ',
**   &           '      wekpo           porel(1)           ddynoc'
**    do j=nypo,1,-1
**      write(*,240) j,po(nxco,j,1),qo(nxco,j,1),wekpo(nxco,j),
**   &               po(nxco,j,1)-poref(1),ddynoc(nxco,j)
**    enddo
*     Southern boundary values
**    write(*,*) ' Southern boundary values for ocean'
**    write(*,*) '    i        po(1)             qo(1)             ',
**   &           'wekpo           porel(1)'
**    do i=1,nxpo
**      write(*,240) i,po(i,1,1),qo(i,1,1),wekpo(i,1),
**   &               po(i,1,1)-poref(1)
**    enddo
      sstmin =  extrem
      sstmax = -extrem
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j)
!$OMP&         SHARED  (sstmin,sstmax,sst)

!$OMP DO SCHEDULE (STATIC)
!$OMP&   REDUCTION(min:sstmin) REDUCTION(max:sstmax)
      do j=1,nyto
        do i=1,nxto
          sstmin = min( sstmin, sst(i,j) )
          sstmax = max( sstmax, sst(i,j) )
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL
      write(*,217) '  s.s.t: min, max = ',sstmin,sstmax
#endif /* atmos_only */

*     Print some atmos. spot values and extrema
*     -----------------------------------------
#ifndef ocean_only
      write(*,217) '  pa(k) at centre = ',(pa(nxca,nyca,k),k=1,nla)
      write(*,217) '  qa(k) at centre = ',(qa(nxca,nyca,k),k=1,nla)
**    write(*,217) '  pa(1) on S. bdy = ',(pa(i,  1 ,1),i=nxca-1,nxca+1)
**    write(*,217) '  qa(1) on S. bdy = ',(qa(i,  1 ,1),i=nxca-1,nxca+1)
**    write(*,217) '  pa(1) on N. bdy = ',(pa(i,nypa,1),i=nxca-1,nxca+1)
**    write(*,217) '  qa(1) on N. bdy = ',(qa(i,nypa,1),i=nxca-1,nxca+1)
**    write(*,217) '   ast  on S. bdy = ',(ast(i,  1 ),i=nxca-1,nxca+1)
**    write(*,217) '   ast  on N. bdy = ',(ast(i,nyta),i=nxca-1,nxca+1)
**    write(*,217) '  Atmos kealat(k) = ',(kealat(k),k=1,nla)
**    write(*,217) '  Atmos et2mat(k) = ',(et2mat(k),k=1,nla-1)
*     Meridional section at centre
**    write(*,*) ' Meridional section at atmos. centre'
**    write(*,'(a,a)') '     j         pa(1)             qa(1)       ',
**   &           '      wekpa           parel(1)           ddynat'
**    do j=nypa,1,-1
**      write(*,240) j,pa(nxca,j,1),qa(nxca,j,1),wekpa(nxca,j),
**   &               pa(nxca,j,1)-paref(1),ddynat(nxca,j)
**    enddo
      astmin =  extrem
      astmax = -extrem
      hmxmin =  extrem
      hmxmax = -extrem
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j)
!$OMP&         SHARED  (astmin,astmax,ast,hmxmin,hmxmax,hmixa)

!$OMP DO SCHEDULE (STATIC)
!$OMP&   REDUCTION(min:astmin) REDUCTION(max:astmax)
!$OMP&   REDUCTION(min:hmxmin) REDUCTION(max:hmxmax)
      do j=1,nyta
        do i=1,nxta
          astmin = min( astmin, ast(i,j) )
          astmax = max( astmax, ast(i,j) )
          hmxmin = min( hmxmin, hmixa(i,j) )
          hmxmax = max( hmxmax, hmixa(i,j) )
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL
      write(*,217) '  a.s.t: min, max = ',astmin,astmax
      write(*,217) '  hmixa: min, max = ',hmxmin,hmxmax
#endif /* not ocean_only */

  217 format(a,1p,9d15.7)
  240 format(i6,1p,5d18.10)

      END SUBROUTINE prsamp
c
c***********************************************************************
c
      SUBROUTINE cfltry (nt)
*
*     Check how well the current state satisfies the CFL criterion.
*     From version 1.5.0, additionally checks atmos. mixed layer.

*     Modules
      USE parameters, ONLY : nxpo, nypo, nlo, nxpa, nypa, nla, fnot
#ifndef ocean_only
      USE atconst, ONLY : dta, dxa, rdxaf0
      USE atstate, ONLY : pa
      USE intrfac, ONLY : uekat, vekat
#endif /* not ocean_only */
#ifndef atmos_only
      USE occonst, ONLY : dto, dxo, rdxof0
      USE ocstate, ONLY : po
      USE intrfac, ONLY : tauxo, tauyo, hmoc
#endif /* not atmos_only */
      USE timinfo, ONLY : tday, tyrs

      IMPLICIT NONE

*     Subroutine arguments
      integer, INTENT(IN) :: nt
*
*     Local parameters
      double precision extrem,cflcrit
      logical cflocn,cflatm,cflmix
      parameter ( extrem=1.0d30, cflocn=.true., cflatm=.true.,
     &            cflmix= .true., cflcrit=0.8d0 )
*
*     Local variables
      integer i,j,k
      double precision uabs,uabmax,vabs,vabmax
#ifndef atmos_only
      double precision uamaxo(nlo),vamaxo(nlo),rhf0hm
#endif /* not atmos_only */
#ifndef ocean_only
      double precision uamaxa(nla),vamaxa(nla)
#endif /* not ocean_only */

      print *,' '
      write(*,'(a,i12,f11.2,f11.4)')
     &      '  CFL testing   at nt, tday, tyrs = ',nt,tday,tyrs

#ifndef atmos_only

*     Check extreme ocean velocity components and derive Courant numbers
*     ==================================================================
      if ( cflocn ) then

        if ( cflmix ) then
*         Optionally also infer um, vm including the
*         ageostrophic components (Ekman transport)
          rhf0hm = 0.5d0/(fnot*hmoc)
          uabmax = -extrem
          vabmax = -extrem
!$OMP     PARALLEL DEFAULT (NONE)
!$OMP&             PRIVATE (i,j,uabs,vabs)
!$OMP&             SHARED  (rdxof0,rhf0hm,po,tauxo,tauyo,uabmax,vabmax)

!$OMP     DO SCHEDULE (STATIC)
!$OMP&       REDUCTION(max:uabmax)
          do j=1,nypo-1
            do i=1,nxpo
              uabs = abs( -rdxof0*(  po(i,j+1,1) -  po(i,j,1) )
     &                    +rhf0hm*( tauyo(i,j+1) + tauyo(i,j) ) )
              uabmax = max( uabmax, uabs )
            enddo
          enddo
!$OMP     END DO NOWAIT
!$OMP     DO SCHEDULE (STATIC)
!$OMP&       REDUCTION(max:vabmax)
          do j=1,nypo
            do i=1,nxpo-1
              vabs = abs(  rdxof0*(  po(i+1,j,1) -  po(i,j,1) )
     &                    -rhf0hm*( tauxo(i+1,j) + tauxo(i,j) ) )
              vabmax = max( vabmax, vabs )
            enddo
          enddo
!$OMP     END DO NOWAIT

!$OMP     END PARALLEL

          write(*,217) '  Ocn. m.l. |u| max, CFL = ',
     &                 uabmax,uabmax*dto/dxo
          write(*,217) '  Ocn. m.l. |v| max, CFL = ',
     &                 vabmax,vabmax*dto/dxo
*         If have bad CFL values, scan for and print locations
          if ( uabmax*dto/dxo.ge.cflcrit ) then
            do j=1,nypo-1
              do i=1,nxpo
                uabs = (dto/dxo)*
     &                 abs( -rdxof0*(  po(i,j+1,1) -  po(i,j,1) )
     &                      +rhf0hm*( tauyo(i,j+1) + tauyo(i,j) ) )
                if ( uabs.ge.cflcrit ) then
                  write(*,250) '  Bad m.l. |u|; CFL, i, j = ',uabs,i,j
                endif
              enddo
            enddo
          endif
          if ( vabmax*dto/dxo.ge.cflcrit ) then
            do j=1,nypo
              do i=1,nxpo-1
                vabs = (dto/dxo)*
     &                 abs(  rdxof0*(  po(i+1,j,1) -  po(i,j,1) )
     &                      -rhf0hm*( tauxo(i+1,j) + tauxo(i,j) ) )
                if ( vabs.ge.cflcrit ) then
                  write(*,250) '  Bad m.l. |v|; CFL, i, j = ',vabs,i,j
                endif
              enddo
            enddo
          endif
        endif

*       Infer u, v in Q-G layers geostrophically
*       ----------------------------------------
        do k=1,nlo
          uabmax = -extrem
          vabmax = -extrem
!$OMP     PARALLEL DEFAULT (NONE)
!$OMP&             PRIVATE (i,j,uabs,vabs)
!$OMP&             SHARED  (k,rdxof0,po,uabmax,vabmax)

!$OMP     DO SCHEDULE (STATIC)
!$OMP&       REDUCTION(max:uabmax)
          do j=1,nypo-1
            do i=1,nxpo
              uabs = abs( -rdxof0*( po(i,j+1,k) - po(i,j,k) ) )
              uabmax = max( uabmax, uabs )
            enddo
          enddo
!$OMP     END DO NOWAIT
!$OMP     DO SCHEDULE (STATIC)
!$OMP&       REDUCTION(max:vabmax)
          do j=1,nypo
            do i=1,nxpo-1
              vabs = abs( rdxof0*( po(i+1,j,k) - po(i,j,k) ) )
              vabmax = max( vabmax, vabs )
            enddo
          enddo
!$OMP     END DO NOWAIT

!$OMP     END PARALLEL
          uamaxo(k) = uabmax
          vamaxo(k) = vabmax
        enddo
        write(*,217) '  Ocn. max geos |u|(k) = ',(uamaxo(k),k=1,nlo)
        write(*,217) '  Ocn. CFL geos |u|(k) = ',
     &               (uamaxo(k)*dto/dxo,k=1,nlo)
        write(*,217) '  Ocn. max geos |v|(k) = ',(vamaxo(k),k=1,nlo)
        write(*,217) '  Ocn. CFL geos |v|(k) = ',
     &               (vamaxo(k)*dto/dxo,k=1,nlo)
*       If have bad CFL values, scan for and print locations
        do k=1,nlo
          if ( uamaxo(k)*dto/dxo.ge.cflcrit ) then
            do j=1,nypo-1
              do i=1,nxpo
                uabs = (dto/dxo)*
     &                 abs( -rdxof0*( po(i,j+1,k) - po(i,j,k) ) )
                if ( uabs.ge.cflcrit ) then
                  write(*,250) '  Bad |u|; CFL, i, j, k = ',uabs,i,j,k
                endif
              enddo
            enddo
          endif
          if ( vamaxo(k)*dto/dxo.ge.cflcrit ) then
            do j=1,nypo
              do i=1,nxpo-1
                vabs = (dto/dxo)*
     &                 abs( rdxof0*( po(i+1,j,k) - po(i,j,k) ) )
                if ( vabs.ge.cflcrit ) then
                  write(*,250) '  Bad |v|; CFL, i, j, k = ',vabs,i,j,k
                endif
              enddo
            enddo
          endif
        enddo
      endif

#endif

#ifndef ocean_only

*     Check extreme atmos. velocity components and derive Courant numbers
*     ===================================================================
      if ( cflatm ) then
        if ( cflmix ) then
*         Optionally also infer um, vm including the ageostrophic
*         components (Ekman transport), which as uekat,
*         vekat, have already been computed in xfosubs.F
          uabmax = -extrem
          vabmax = -extrem
!$OMP     PARALLEL DEFAULT (NONE)
!$OMP&             PRIVATE (i,j,uabs,vabs)
!$OMP&             SHARED  (rdxaf0,pa,uekat,vekat,uabmax,vabmax)

!$OMP     DO SCHEDULE (STATIC)
!$OMP&       REDUCTION(max:uabmax)
          do j=1,nypa-1
            do i=1,nxpa
              uabs = abs( -rdxaf0*( pa(i,j+1,1) -  pa(i,j,1) )
     &                   + uekat(i,j) )
              uabmax = max( uabmax, uabs )
            enddo
          enddo
!$OMP     END DO NOWAIT
!$OMP     DO SCHEDULE (STATIC)
!$OMP&       REDUCTION(max:vabmax)
          do j=1,nypa
            do i=1,nxpa-1
              vabs = abs(  rdxaf0*( pa(i+1,j,1) -  pa(i,j,1) )
     &                   + vekat(i,j) )
              vabmax = max( vabmax, vabs )
            enddo
          enddo
!$OMP     END DO NOWAIT

!$OMP     END PARALLEL
          write(*,217) '  Atm. m.l. |u| max, CFL = ',
     &                 uabmax,uabmax*dta/dxa
          write(*,217) '  Atm. m.l. |v| max, CFL = ',
     &                 vabmax,vabmax*dta/dxa
*         If have bad CFL values, scan for and print locations
          if ( uabmax*dta/dxa.ge.cflcrit ) then
            do j=1,nypa-1
              do i=1,nxpa
                uabs = (dta/dxa)*
     &                 abs( -rdxaf0*( pa(i,j+1,1) -  pa(i,j,1) )
     &                     + uekat(i,j) )
                if ( uabs.ge.cflcrit ) then
                  write(*,250) '  Bad m.l. |u|; CFL, i, j = ',uabs,i,j
                endif
              enddo
            enddo
          endif
          if ( vabmax*dta/dxa.ge.cflcrit ) then
            do j=1,nypa
              do i=1,nxpa-1
                vabs = (dta/dxa)*
     &                 abs(  rdxaf0*( pa(i+1,j,1) -  pa(i,j,1) )
     &                     + vekat(i,j) )
                if ( vabs.ge.cflcrit ) then
                  write(*,250) '  Bad m.l. |v|; CFL, i, j = ',vabs,i,j
                endif
              enddo
            enddo
          endif
        endif

*       Infer u, v in Q-G layers geostrophically
*       ----------------------------------------
        do k=1,nla
          uabmax = -extrem
          vabmax = -extrem
!$OMP     PARALLEL DEFAULT (NONE)
!$OMP&             PRIVATE (i,j,uabs,vabs)
!$OMP&             SHARED  (k,rdxaf0,pa,uabmax,vabmax)

!$OMP     DO SCHEDULE (STATIC)
!$OMP&       REDUCTION(max:uabmax)
          do j=1,nypa-1
            do i=1,nxpa
              uabs = abs( -rdxaf0*( pa(i,j+1,k) - pa(i,j,k) ) )
              uabmax = max( uabmax, uabs )
            enddo
          enddo
!$OMP     END DO NOWAIT
!$OMP     DO SCHEDULE (STATIC)
!$OMP&       REDUCTION(max:vabmax)
          do j=1,nypa
            do i=1,nxpa-1
              vabs = abs( rdxaf0*( pa(i+1,j,k) - pa(i,j,k) ) )
              vabmax = max( vabmax, vabs )
            enddo
          enddo
!$OMP     END DO NOWAIT

!$OMP     END PARALLEL
          uamaxa(k) = uabmax
          vamaxa(k) = vabmax
        enddo
        write(*,217) '  Atm. max geos |u|(k) = ',(uamaxa(k),k=1,nla)
        write(*,217) '  Atm. CFL geos |u|(k) = ',
     &               (uamaxa(k)*dta/dxa,k=1,nla)
        write(*,217) '  Atm. max geos |v|(k) = ',(vamaxa(k),k=1,nla)
        write(*,217) '  Atm. CFL geos |v|(k) = ',
     &               (vamaxa(k)*dta/dxa,k=1,nla)
*       If have bad CFL values, scan for and print locations
        do k=1,nla
          if ( uamaxa(k)*dta/dxa.ge.cflcrit ) then
            do j=1,nypa-1
              do i=1,nxpa
                uabs = (dta/dxa)*
     &                 abs( -rdxaf0*( pa(i,j+1,k) - pa(i,j,k) ) )
                if ( uabs.ge.cflcrit ) then
                  write(*,250) '  Bad |u|; CFL, i, j, k = ',uabs,i,j,k
                endif
              enddo
            enddo
          endif
          if ( vamaxa(k)*dta/dxa.ge.cflcrit ) then
            do j=1,nypa
              do i=1,nxpa-1
                vabs = (dta/dxa)*
     &                 abs( rdxaf0*( pa(i+1,j,k) - pa(i,j,k) ) )
                if ( vabs.ge.cflcrit ) then
                  write(*,250) '  Bad |v|; CFL, i, j, k = ',vabs,i,j,k
                endif
              enddo
            enddo
          endif
        enddo
      endif

#endif /* not ocean_only */

  217 format(a,1p,9d15.7)
  250 format(a,1p,d15.7,2i8,i6)

      END SUBROUTINE cfltry
c
c***********************************************************************
c
      SUBROUTINE memreq

*     Compute size of ocean and atmosphere arrays, data module sizes,
*     and local memory allocation in subroutines if significant
*     and required (controlled by internal parameter stacks)

*     Modules
      USE parameters
#ifndef ocean_only
      USE atisubs
#endif
#ifndef atmos_only
      USE ocisubs
#endif
*     Use here those *.mod files of program modules
*     which define array-dimensioning parameters

      IMPLICIT NONE

*     Subroutine arguments
*
*     Local parameters
      double precision dim2mb
      logical stacks
      parameter ( dim2mb=1.0d0/131072.0d0, stacks = .false. )
*     dim2mb converts from no. of elements in a double
*     precision array to size in Megabytes (= 8/(1024**2) )
*
*     Local variables
      double precision sizarr,dynsiz,dgnsiz
      double precision sizmin

*     Compute sizes of main arrays
*     ============================
      print *,' '
      write(*,'(a)') '  Array sizes:'
      write(*,'(a)') '  ------------'
      sizarr = nxpo*nypo
      write(*,203) '  Ocean  2D p-array size (MB) = ',dim2mb*sizarr
      sizarr = nxpa*nypa
      write(*,203) '  Atmos. 2D p-array size (MB) = ',dim2mb*sizarr
*     Compute size of full 3D arrays
#ifndef atmos_only
      sizarr = nxpo*nypo*nlo
      write(*,203) '  Ocean  3D p-array size (MB) = ',dim2mb*sizarr
#endif
#ifndef ocean_only
      sizarr = nxpa*nypa*nla
      write(*,203) '  Atmos. 3D p-array size (MB) = ',dim2mb*sizarr
#endif
      sizarr = nxpaor*nypaor
      write(*,203) '  At higher (ocean) resolution,'
      write(*,203) '  atmos. 2D p-array size (MB) = ',dim2mb*sizarr
#ifdef get_covar
*     Compute size of covariance arrays
#  ifndef atmos_only
      sizarr = nmcvoc
      write(*,203) '  Ocean covar array size (MB) = ',dim2mb*sizarr
      sizarr = nvcvoc
      write(*,203) '  Ocean avge  array size (MB) = ',dim2mb*sizarr
#  endif
#  ifndef ocean_only
      sizarr = nmcvat
      write(*,203) '  Atmos covar array size (MB) = ',dim2mb*sizarr
      sizarr = nvcvat
      write(*,203) '  Atmos avge  array size (MB) = ',dim2mb*sizarr
#  endif
      write(*,201) ' '
#  ifndef atmos_only
      write(*,201) '  Ocean covar array dimension = ',nmcvoc
      write(*,201) '  Ocean avge  array dimension = ',nvcvoc
#  endif
#  ifndef ocean_only
      write(*,201) '  Atmos covar array dimension = ',nmcvat
      write(*,201) '  Atmoc avge  array dimension = ',nvcvat
#  endif
#endif /* get_covar */

*     Compute sizes of data storage in modules
*     ========================================
*     Accumulate the number of double precision variables
      dynsiz = 0.0d0
      dgnsiz = 0.0d0
      print *,' '
      write(*,'(a)') '  Data sizes of modules:'
      write(*,'(a)') '  ======================'
      write(*,'(a)') '  Dynamics:'
      write(*,'(a)') '  ---------'
      write(*,'(a)') '  Data modules:'
*     atconst data module
*     -------------------
      sizarr = 7 + nxpa + 2*nypa + nxta + 2*nyta + 2 + 3*nla + 4
     &        + 2*nla + nxta*nyta + 1 + nla*nla + 3*nla + 2*nla*nla
     &        + 2*nxpa*nypa + 1
      dynsiz = dynsiz + sizarr
      write(*,203) '  atconst mod total size (MB) = ',dim2mb*sizarr
*     athomog data module
*     -------------------
#ifndef ocean_only
      sizarr = 2*nypa*(nla-1) + 5*(nla-1) + nypa + 2 + 3*(nla-1)
     &        + 4*nla + 2*(nla-1) + 4*nla + 2
      dynsiz = dynsiz + sizarr
      write(*,203) '  athomog mod total size (MB) = ',dim2mb*sizarr
#endif
*     atstate data module
*     -------------------
      sizarr = 4*nxpa*nypa*nla + 2*nxpa*nypa + nxta*nyta
      dynsiz = dynsiz + sizarr
      write(*,203) '  atstate mod total size (MB) = ',dim2mb*sizarr
*     intrfac data module
*     -------------------
      sizarr = 2*nxto*nyto + nyto + 2*nxta*nyta + nyta + 2*nxpo*nypo
     &        + 2*nxpa*nypa + nxpa*nyta + nxta*nypa + nxpa*(nypa+1)
     &        + (nxpa+1)*nypa + nxto*nyto + 3*nxta*nyta + 13
      dynsiz = dynsiz + sizarr
      write(*,203) '  intrfac mod total size (MB) = ',dim2mb*sizarr
*     occonst data module
*     -------------------
      sizarr = 7 + nxpo + 2*nypo + nxto + 2*nyto + 2 + (nlo-1)
     &        + 2*nlo + 2 + 3*nlo + 5 + nlo*nlo + 3*nlo + 2*nlo*nlo
     &        + 2*nxpo*nypo + 1
      dynsiz = dynsiz + sizarr
      write(*,203) '  occonst mod total size (MB) = ',dim2mb*sizarr
*     ochomog data module
*     -------------------
#ifndef atmos_only
#  ifdef cyclic_ocean
*     Zonally cyclic ocean
      sizarr = 2*nypo*(nlo-1) + 4*(nlo-1) + nypo + (nlo-1) + 2
     &        + 4*nlo + 2*(nlo-1) + 6*nlo + 4
#  else
*     Finite box ocean
      sizarr = nxpo*nypo*(nlo-1) + (nlo-1) + nlo*(nlo-1)
     &        + 2*(nlo-1)*(nlo-1) + (nlo-1)/2
#  endif
      sizarr = sizarr + 3*(nlo-1)
      dynsiz = dynsiz + sizarr
      write(*,203) '  ochomog mod total size (MB) = ',dim2mb*sizarr
#endif
*     ocstate data module
*     -------------------
#ifndef atmos_only
      sizarr = 4*nxpo*nypo*nlo + 2*nxpo*nypo + nxto*nyto
      dynsiz = dynsiz + sizarr
      write(*,203) '  ocstate mod total size (MB) = ',dim2mb*sizarr
#endif
*     radiate data module
*     -------------------
      sizarr = 4 + nla + 4 + 2*nla*(nla-1) + 4 + 2*nla
     &        + 3 + nla + 2*(nla-1) + 3
      dynsiz = dynsiz + sizarr
      write(*,203) '  radiate mod total size (MB) = ',dim2mb*sizarr
*     timinfo data module
*     -------------------
      sizarr = 3 + 14 + (17+1)/2
#ifdef get_covar
      sizarr = sizarr + 1
#endif
      dynsiz = dynsiz + sizarr
      write(*,203) '  timinfo mod total size (MB) = ',dim2mb*sizarr
      write(*,'(a)') '  Program modules:'
*     atisubs prog module
*     -------------------
#ifndef ocean_only
      sizarr = lwftat + 1 + nxta
      dynsiz = dynsiz + sizarr
      write(*,203) '  atisubs mod total size (MB) = ',dim2mb*sizarr
#endif
*     ocisubs prog module
*     -------------------
#ifndef atmos_only
      sizarr = lwftoc + 1 + nxto
      dynsiz = dynsiz + sizarr
      write(*,203) '  ocisubs mod total size (MB) = ',dim2mb*sizarr
#endif
!!    [UP TO HERE]
*     xfosubs prog module
*     -------------------
      sizarr = 5*16*16 + 16*(ndxr+1)*(ndxr+1)
      dynsiz = dynsiz + sizarr
      write(*,203) '  xfosubs mod total size (MB) = ',dim2mb*sizarr

      if ( dim2mb*dynsiz.lt.10000.00 ) then
        write(*,'(a)') '                                  --------'
        write(*,203) '       Dynamics subtotal (MB) = ',dim2mb*dynsiz
        write(*,'(a)') '                                  --------'
       else
        write(*,'(a)') '                                 ---------'
        write(*,203) '       Dynamics subtotal (MB) = ',dim2mb*dynsiz
        write(*,'(a)') '                                 ---------'
      endif
      write(*,'(a)') '  Diagnostics:'
#ifdef get_covar
*     covaria prog module
*     -------------------
      sizarr = 0.5d0
#  ifndef ocean_only
      sizarr = sizarr + 2*( nmcvat + nvcvat + 1 ) + 1
#  endif
#  ifndef atmos_only
      sizarr = sizarr + 2*( nmcvoc + nvcvoc + 1 ) + 1
#  endif
      dgnsiz = dgnsiz + sizarr
      write(*,203) '  covaria mod total size (MB) = ',dim2mb*sizarr
#endif /* get_covar */
*     monitor prog module
*     -------------------
      sizarr = 9 + (nla-1) + 1 + (nla-1) + (nlo-1) + (nla-1) + 17
     &        + 2*(nla-1) + 1 + 2*(nlo-1) + 3*nlo + 3*nla + 2
     &        + 2*nlo + nla + (nlo-1)+ 2*(nla-1) + 5
     &        + (nla-1) + (nlo-1) + 3*nlo + nla + ( nla + nlo )/2
      dgnsiz = dgnsiz + sizarr
      write(*,203) '  monitor mod total size (MB) = ',dim2mb*sizarr
*     nc_subs prog module
*     -------------------
*     nc_subs: non-stack storage size is negligible
*     timavge prog module
*     -------------------
      sizarr = 0.5d0
#ifndef ocean_only
      sizarr =  sizarr + 2*nxpa*nypa + 3*nxta*nyta
     &        + 2*nxpa*nypa*nla + 3*nxpa*nyta + 3*nxta*nypa
#endif
#ifndef atmos_only
      sizarr =  sizarr + 3*nxpo*nypo + 3*nxto*nyto
     &        + 2*nxpo*nypo*nlo + 3*nxpo*nyto + 3*nxto*nypo + 2
#endif
      dgnsiz = dgnsiz + sizarr
      write(*,203) '  timavge mod total size (MB) = ',dim2mb*sizarr
      if ( dim2mb*dgnsiz.lt.10000.00 ) then
        write(*,'(a)') '                                  --------'
       else
        write(*,'(a)') '                                 ---------'
      endif
      write(*,203) '    Diagnostics subtotal (MB) = ',dim2mb*dgnsiz
      if ( dim2mb*(dynsiz+dgnsiz).lt.10000.00 ) then
        write(*,'(a)') '                                  --------'
        write(*,203) '  Total common blk  size (MB) = ',
     &               dim2mb*(dynsiz+dgnsiz)
        write(*,'(a)') '                                  ========'
       else
        write(*,'(a)') '                                 ---------'
        write(*,203) '  Total common blk  size (MB) = ',
     &               dim2mb*(dynsiz+dgnsiz)
        write(*,'(a)') '                                 ========='
      endif

      if ( stacks ) then
*     Compute size of local arrays/vectors
*     ====================================
*     Intended to be in order of first execution, to
*     investigate possible stack memory allocation problems
*     (N.B. former "small_local" option now removed)
      print *,' '
      print *,' Stack allocation of local arrays/vectors:'
      print *,' -----------------------------------------'
      print *,' Setup phase:'

*     constr
*     ------
      sizmin = 4
#  ifndef ocean_only
      sizmin = sizmin + nxpa*nypa*nla + 4*nla
#  endif
#  ifndef atmos_only
      sizmin = sizmin + 4*nlo
      sizarr = sizmin
      sizarr = sizarr + nxpo*nypo*nlo
#  endif
      write(*,203) '  constr     local memory   (MB) = ',dim2mb*sizarr
#  ifndef atmos_only
      write(*,203) '  constr     loc. with SAVE (MB) = ',dim2mb*sizmin
#  endif

*     xforc
*     -----
      sizmin = 8
#  ifndef ocean_only
      sizmin = sizmin + 5 + 34 + 1 + ndxr+1
#    ifdef tau_udiff
      sizmin = sizmin + 2
#    endif
      sizarr = sizmin + 2*nxpa*nypa + 4*nxpaor*nypaor + nxtaor*nytaor
     &        + nxto*nyto
#  endif
      write(*,203) '  xforc      local memory   (MB) = ',dim2mb*sizarr
#  ifndef ocean_only
      write(*,203) '  xforc      loc. with SAVE (MB) = ',dim2mb*sizmin
#  endif

*     bilint
*     ------
      sizarr = 7 + 3*nxpo
*     (scalars + 2 x integer + 2 x double)
      write(*,203) '  bilint     local memory   (MB) = ',dim2mb*sizarr

*     auvbli
*     ------
***   sizarr = 5 + 3*nxpaor
***   write(*,203) '  auvbli     local memory   (MB) = ',dim2mb*sizarr

*     auvbcu
*     ------
      sizarr = 5 + 4*16 + 6
      write(*,203) '  auvbcu     local memory   (MB) = ',dim2mb*sizarr

*     local common block bcudata
*     --------------------------
      sizarr = 5*16*16 + 16*(ndxr+1)*(ndxr+1)
      write(*,203) '  bcudata    local memory   (MB) = ',dim2mb*sizarr

*     bcuini
*     ------
      sizarr = 4*4*4*2*2 + 2 + 4
      write(*,203) '  bcuini     local memory   (MB) = ',dim2mb*sizarr

*     wts2bb
*     ------
      sizarr = 2*16*16 + 1 + 5
      write(*,203) '  wts2bb     local memory   (MB) = ',dim2mb*sizarr

*     homsol
*     ------
      sizmin = 7
#  ifndef ocean_only
      sizmin = sizmin + 2*nxpa*nypa + nxta
#  endif
#  ifndef atmos_only
      sizmin = sizmin + nxto
      sizarr = sizmin
#    ifdef cyclic_ocean
      sizarr = sizarr + 2*nxpo*nypo
#    endif
#  endif
      write(*,203) '  homsol     local memory   (MB) = ',dim2mb*sizarr
#  ifndef atmos_only
#    ifdef cyclic_ocean
      write(*,203) '  homsol     loc. with SAVE (MB) = ',dim2mb*sizmin
#    endif
#  endif

*     hsbxoc/hscyoc
*     -------------
#  ifndef atmos_only
#    ifdef cyclic_ocean
      sizarr = lwftoc + 2*(nypo-1) + 1
      write(*,203) '  hscyoc     local memory   (MB) = ',dim2mb*sizarr
#    else
      sizarr = lwftoc + 2*(nypo-1) + 1
      write(*,203) '  hsbxoc     local memory   (MB) = ',dim2mb*sizarr
#    endif
#  endif

*     monnc_comp
*     ----------
      sizmin = 21
#  ifndef ocean_only
      sizmin = sizmin + 3*nxta*nypa + 3*nxpa*nypa + nxpa*nypa
     &        + nxpa*nyta + nxta*nypa + nyta + nla
#  endif
      sizarr = sizmin
#  ifndef atmos_only
      sizmin = sizmin + nyto + nlo + 2
      sizarr = sizmin
      sizarr = sizarr + 4*nxto*nypo + 4*nxpo*nypo
     &        + nxpo*nypo + nxpo*nyto + nxto*nypo
#  endif
      write(*,203) '  monnc_comp local memory   (MB) = ',dim2mb*sizarr
#  ifndef atmos_only
      write(*,203) '  monnc_comp loc. with SAVE (MB) = ',dim2mb*sizmin
#  endif

#  ifdef use_netcdf
#    ifndef atmos_only
*     ocnc_init
      sizarr = nxpo + nypo + nlo
      write(*,203) '  ocnc_init  local arrays   (MB) = ',dim2mb*sizarr
*     ocnc_out
      sizarr = nxpo*nypo*nlo
      write(*,203) '  ocnc_out   local arrays   (MB) = ',dim2mb*sizarr
      sizarr = 0
      write(*,203) '  ocnc_out   loc. with SAVE (MB) = ',dim2mb*sizarr
#  ifdef qoc_diag
*     qocdiag_init
      sizarr = nxpo + nypo + nlo
      write(*,203) '  qocdiag_init loc. arrays (MB) = ',dim2mb*sizarr
*     qocdiag_out
      sizarr = 8*nxpo*nypo + nlo + 8
      write(*,203) '  qocdiag_out local arrays (MB) = ',dim2mb*sizarr
      sizarr = nlo + 8
      write(*,203) '  qocdiag_out   with SAVE  (MB) = ',dim2mb*sizarr
#  endif /* qoc_diag */
#    endif
#    ifndef ocean_only
*     atnc_init
      sizarr = nxpa + nypa + nla
      write(*,203) '  atnc_init local arrays   (MB) = ',dim2mb*sizarr
*     atnc_out
      sizarr = nxpa*nypa*nla
      write(*,203) '  atnc_out  local arrays   (MB) = ',dim2mb*sizarr
#    endif
#  endif
      print *,' Timestepping phase:'
#  ifndef atmos_only
*     oml
      sizarr = 2*nxto*nyto
      write(*,203) '  oml        local arrays   (MB) = ',dim2mb*sizarr
      sizarr = 0
      write(*,203) '  oml        loc. with SAVE (MB) = ',dim2mb*sizarr
*     omladf
      sizarr = (nxto+2)*nyto
      write(*,203) '  omladf     local arrays   (MB) = ',dim2mb*sizarr
      sizarr = 0
      write(*,203) '  omladf     loc. with SAVE (MB) = ',dim2mb*sizarr
*     qgostep
      sizarr = 4 + 2*nlo + nxpo*nypo + nxpo*nypo*nlo
      write(*,203) '  qgostep    local arrays   (MB) = ',dim2mb*sizarr
*     ocadif
      sizarr = 4 + nxpo*nypo
#  ifdef cyclic_ocean
      sizarr = sizarr + 8
#  endif
      write(*,203) '  ocadif     local arrays   (MB) = ',dim2mb*sizarr
*     ocinvq
      sizarr = 3*nlo + nxpo*nypo*nlo + nxto
#    ifdef cyclic_ocean
      sizarr = sizarr + 9*nlo + 2*(nlo-1)
#    else
      sizarr = sizarr + 3*nlo + 3*(nlo-1)
#    endif
      write(*,203) '  ocinvq     local arrays   (MB) = ',dim2mb*sizarr
      sizarr = sizarr - nxpo*nypo*nlo
      write(*,203) '  ocinvq     loc. with SAVE (MB) = ',dim2mb*sizarr
#  endif
#  ifndef ocean_only
*     aml
      sizarr = 2*(nla-1) + 3*nxta*nyta
      write(*,203) '  aml        local arrays   (MB) = ',dim2mb*sizarr
*     amladf
      sizarr = (nxta+2)*nyta
      write(*,203) '  amladf     local arrays   (MB) = ',dim2mb*sizarr
*     qgastep
      sizarr = 2*nla + nxpa*nypa + nxpa*nypa*nla
      write(*,203) '  qgastep    local arrays   (MB) = ',dim2mb*sizarr
*     atadif
      sizarr = 4*nxpa*nypa
      write(*,203) '  atadif     local arrays   (MB) = ',dim2mb*sizarr
*     atinvq
      sizarr = 14*nla + nxpa*nypa*nla + nxta + 2*(nla-1)
      write(*,203) '  atinvq     local arrays   (MB) = ',dim2mb*sizarr
*     tavatm
      sizarr = nxpa
      write(*,203) '  tavatm     local arrays   (MB) = ',dim2mb*sizarr
#  endif
#  ifndef atmos_only
      sizarr = 3*nxpo
      write(*,203) '  tavocn     local arrays   (MB) = ',dim2mb*sizarr
#  endif
      print *,' Finishing phase:'

*     tavout
*     ------
      sizmin = 0.0d0
#  ifndef ocean_only
      sizmin = sizmin + 1 + nxpa*nyta + nxta*nypa
#  endif
#  ifdef use_netcdf
#    ifndef ocean_only
      sizmin = sizmin + nxpa + nypa + nla
#    endif
#    ifndef atmos_only
      sizmin = sizmin + nxpo + nypo + nlo
#    endif
#  endif
      sizarr = sizmin
#  ifndef atmos_only
      sizmin = sizmin + 1
      sizarr = sizmin
      sizarr = sizarr + nxpo*nyto + nxto*nypo
#  endif
      write(*,203) '  tavout     local memory   (MB) = ',dim2mb*sizarr
#  ifndef atmos_only
      write(*,203) '  tavout     loc. with SAVE (MB) = ',dim2mb*sizmin
#  endif
      endif

  201 format(a,i12)
*     Integer format should be large enough to print max. I*4 number
  203 format(a,f10.3)

      END SUBROUTINE memreq
c
c***********************************************************************
c
      SUBROUTINE diffts (nord, nl, coeff, ncoef, dx, rdef, dt)

*     Computes diffusive decay timescales for circular eddies whose radii
*     are the baroclinic Rossby radii, and for two-gridpoint noise.
*     See section 8.6 of the Userguide for derivation of timescales.
*     Also computes the nondimensional timestep stability criterion.

*     Input arguments:
*     nord  : order of the diffusive term
*     nl    : no. of QG layers  (=> nl-1 baroclinic modes)
*     coeff : vector of diffusion coefficients (should be .GE. 0)
*     ncoef : length of coefficient vector
*     dx    : gridlength (m)
*     rdef  : vector of nl modal deformation radii (m)
*             (infinite value for barotropic mode replaced by 0.0)
*     dt    : timestep (sec)
*     (all the above are unchanged on exit)

*     Modules
      USE parameters

      IMPLICIT NONE

*     Subroutine arguments
      integer, INTENT(IN) :: nord,nl,ncoef
      double precision, INTENT(IN) :: coeff(ncoef),dx,rdef(nl),dt
*
*     Local parameters
      DOUBLE PRECISION PIBY2
      PARAMETER ( PIBY2=1.57079632679489662D0 )
      integer nlmax
      double precision secday
      parameter ( nlmax=9, secday=86400.0d0 )
*
*     Local variables
      integer k,m
      double precision tdamp(nlmax),sinfac,dtstab(nlmax)

*     Check internal storage is sufficient
      if ( nl.gt.nlmax ) then
        print *,' diffts has insufficient nlmax = ',nlmax
        print *,' called with nl = ',nl
        print *,' program terminates in diffts'
        stop
      endif

*     Check all diffusion coefficients are non-negative
*     (need positive coeffts for damping)
      do k=1,ncoef
        if ( coeff(k).lt.0.0d0 ) then
          print *,' diffts called with -ve diffusion coefft'
          print *,' coeff vector = ',(coeff(m),m=1,ncoef)
          print *,' program terminates in diffts'
          stop
        endif
      enddo

*     Compute decay timescale(s) for a circular eddy
*     at the deformation radius for each baroclinic mode
*     Deformation radii are usually much greater than
*     the gridlength dx, but for very coarse resolution
*     cases where rdef(m) < dx, use dx instead of
*     rdef(m) to compute the modal decay timescale
      do m=2,nl
        if ( dx.le.rdef(m) ) then
          sinfac = 2.0d0*sin( PIBY2*dx/rdef(m) )/dx
         else
          sinfac = 2.0d0/dx
          write(*,225) '  NOTE: mode',m-1,
     &                 ' defrad < dx; dx used instead'
          write(*,205) '  of defrad for computing damping timescale'
        endif
*       Avoid infinities if coefft = 0
        do k=1,ncoef
          if ( coeff(k).eq.0.0d0 ) then
            tdamp(k) = 0.0d0
           else
            tdamp(k) = 1.0d0/( sinfac**nord*coeff(k)*dble(nord)*secday )
          endif
        enddo
        write(*,225) '  Mode',m-1,' damping time  (days) = ',
     &               (tdamp(k),k=1,ncoef)
      enddo

*     Compute decay timescale for 1-D two-gridpoint noise
*     for each coefft, avoiding infinities if coefft = 0
*     The timescale for 2-D (checkerboard) noise differs only
*     by a factor of nord. This has the shortest decay timescale,
*     and thus determines the overall timestep stability limit.
*     Also compute nondimensional timestep stability factor.
*     This should be < 1 for stable leapfrog timestepping.
      do k=1,ncoef
        if ( coeff(k).eq.0.0d0 ) then
          tdamp(k) = 0.0d0
          dtstab(k) = 0.0d0
         else
          tdamp(k) = (0.5d0*dx)**nord/coeff(k)
          dtstab(k) = nord*dt/tdamp(k)
*         Rescale damping timescale to hours
          tdamp(k) = tdamp(k)/3600.0d0
        endif
      enddo
      write(*,205) '  1-D grid timescale  (hours) = ',
     &             (tdamp(k),k=1,ncoef)
      write(*,205) '  2-D grid timescale  (hours) = ',
     &             (tdamp(k)/nord,k=1,ncoef)
      write(*,205) '  Timestep stability factor   = ',
     &             (dtstab(k),k=1,ncoef)

  205 format(a,9f13.5)
  225 format(a,i2,a,9f13.5)

      END SUBROUTINE diffts
c
c***********************************************************************
c
      SUBROUTINE resave (iounit)

*     Write an unformatted restart dump, either at end
*     of run, or periodically in case of program crash.
*     Dump contains time in years, and all progostic fields

*     Modules
      USE parameters
#ifndef ocean_only
      USE atstate
#endif
#ifndef atmos_only
      USE ocstate
#endif
      USE intrfac, ONLY : sst, sstm, ast, astm, hmixa, hmixam
      USE timinfo, ONLY : tyrs

      IMPLICIT NONE

*     Subroutine arguments
      integer iounit
*
*     Dump restart file on unit iounit
      write (iounit) tyrs
#ifndef atmos_only
      write (iounit) po,pom
#endif
#ifndef ocean_only
      write (iounit) pa,pam
#endif
      write (iounit) sst,sstm
      write (iounit) ast,astm
      write (iounit) hmixa,hmixam

      END SUBROUTINE resave
