c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 25/09/2013
c***********************************************************************
c
c     Copyright 2013 Jeff Blundell, Andy Hogg and Bill Dewar.
c     This file is part of Q-GCM.
c
c     Q-GCM is free software: you can redistribute it and/or modify
c     it under the terms of the GNU General Public License as
c     published by the Free Software Foundation, either version 3
c     of the License, or (at your option) any later version.
c
c     Q-GCM is distributed in the hope that it will be useful,
c     but WITHOUT ANY WARRANTY; without even the implied warranty
c     of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
c     See the GNU General Public License for more details.
c
c     You should have received a copy of the GNU General Public License
c     along with Q-GCM.  If not, see <http://www.gnu.org/licenses/>.
c
c***********************************************************************
c
      MODULE omuse_valsubs

*     Contains subroutines valids, scan2D and scan3D which check the
*     validity of the current model state against plausible ranges
*     for selected variables, in the hope of finding the location
*     of model blow-ups, and preventing degeneration into NaNs.

*     Modules

      IMPLICIT NONE

      PRIVATE

      PUBLIC  :: valids, hcheck
      PRIVATE :: scan2D, scan3D

      logical :: hcheck=.TRUE.

      CONTAINS

c***********************************************************************
c
      SUBROUTINE valids (solnok)

*     Scan prognostic variable arrays and some forcing fields for
*     extreme values, and set logical argument solnok to .false.
*     if any are located. Diagnosis and termination of program can
*     then be controlled by the calling program. Intended to be cheap
*     to call, to enable frequent checking of solution validity.
*     Note that in atmos_only or ocean_only configurations, not
*     all forcing or boundary fields will be checked; it is assumed
*     that this will have been done in setting up these fields.
*
*     This version (v1.5.0 onwards) includes additional
*     checking of full (perturbed) ocean layer thicknesses.

*     Modules
      USE parameters
#ifndef ocean_only
      USE atstate, ONLY : pa, qa, wekta
#endif /* not ocean_only */
#ifndef atmos_only
      USE ocstate, ONLY : po, qo, wekto
*     Extra information for layer thickness checking
      USE occonst, ONLY : gpoc, hoc, dtopoc
#endif /* not atmos_only */
      USE intrfac, ONLY : sst, ast, tauxa, tauya

      IMPLICIT NONE

*     Subroutine arguments
      logical solnok
*
*     Local parameters
      double precision bignum,tauext,wtaext,wtoext,patext,
     &                 qatext,astext,pocext,qocext,sstext
      parameter ( bignum = 1.0d30, tauext = 10.0d0,
     &            wtaext = 1.0d0, wtoext = 1.0d-3,
     &            astext=90.0d0, patext=1.0d7, qatext=0.05d0,
***  &            sstext=50.0d0, pocext=1.0d4, qocext=0.05d0 )
     &            sstext=75.0d0, pocext=1.0d4, qocext=0.05d0 )
*
*     Extra parameters for layer thickness checking:
*     hcheck turns thickness checking on/off
*     thkmin = min. acceptable layer thickness (m)
*     spfail turns single point failure mode on/off
*     if hcheck is .true. and spfail is .false., a failure
*     criterion based on %age of unacceptable points is used
*     critpc = max. acceptable %age of excessively thin points
*     jwid = no. of points (in j-direction) on either side of
*            the point of failure printed out in the diagnostics
*     Equivalent i-width is hard-wired into the code (currently 3)
*     because it has to be consistent with format statements
      logical spfail
      double precision thkmin,critpc
      integer jwid
***   parameter ( hcheck = .true., thkmin = 75.0d0, spfail = .false.,
      parameter ( thkmin = 100.0d0, spfail = .false.,
     &            critpc = 20.0d0, jwid = 4 )
*
*     Local variables
      integer i,j,k
      double precision wekmin,wekmax
#ifndef ocean_only
      double precision patmin,patmax,qatmin,qatmax,astmin,astmax,
     &                 txamin,txamax,tyamin,tyamax
#endif /* not ocean_only */
#ifndef atmos_only
      double precision pocmin,pocmax,qocmin,qocmax,sstmin,sstmax
*
*     Extra variables for layer thickness checking
      logical pcfail,hffail
      integer im,jm,km,ilo,ihi,ioff,jlo,jhi,joff
      double precision rgpoc(nlo-1),etaoc(nlo-1),hfull,hfmint,hfmaxt,
     &                 hfmini,hfmaxi,hfminb,hfmaxb,hfmina,hfmaxa,
     &                 hfarr(7,2*jwid+1),etaarr(7,2*jwid+1,nlo-1),
     &                 hfbad(nlo),wti,wtj
#endif /* not atmos_only */

#ifndef ocean_only
*     Check atmospheric fields for extreme values
*     ===========================================
*     Values in the neighbourhood of invalid values
*     in any of the fields scanned will be printed
      patmin =  bignum
      patmax = -bignum
      qatmin =  bignum
      qatmax = -bignum
      astmin =  bignum
      astmax = -bignum
      wekmin =  bignum
      wekmax = -bignum
      txamin =  bignum
      txamax = -bignum
      tyamin =  bignum
      tyamax = -bignum

!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k)
!$OMP&         SHARED  (patmin,patmax,pa,qatmin,qatmax,qa,
!$OMP&                  astmin,astmax,ast,wekmin,wekmax,wekta,
!$OMP&                  txamin,txamax,tauxa,tyamin,tyamax,tauya)

      do k=1,nla
!$OMP   DO SCHEDULE (STATIC)
!$OMP&     REDUCTION(min:patmin,qatmin) REDUCTION(max:patmax,qatmax)
        do j=1,nypa
          do i=1,nxpa
            patmin = min( patmin, pa(i,j,k) )
            patmax = max( patmax, pa(i,j,k) )
            qatmin = min( qatmin, qa(i,j,k) )
            qatmax = max( qatmax, qa(i,j,k) )
          enddo
        enddo
!$OMP   END DO
      enddo

!$OMP DO SCHEDULE (STATIC)
!$OMP&   REDUCTION(min:astmin,wekmin) REDUCTION(max:astmax,wekmax)
      do j=1,nyta
        do i=1,nxta
          astmin = min( astmin, ast(i,j) )
          astmax = max( astmax, ast(i,j) )
          wekmin = min( wekmin, wekta(i,j) )
          wekmax = max( wekmax, wekta(i,j) )
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP DO SCHEDULE (STATIC)
!$OMP&   REDUCTION(min:txamin,tyamin) REDUCTION(max:txamax,tyamax)
      do j=1,nypa
        do i=1,nxpa
          txamin = min( txamin, tauxa(i,j) )
          txamax = max( txamax, tauxa(i,j) )
          tyamin = min( tyamin, tauya(i,j) )
          tyamax = max( tyamax, tauya(i,j) )
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

**    print *,' pa: min, max = ',patmin,patmax
      if ( abs(patmin).ge.patext .or. abs(patmax).ge.patext ) then
        print *,' '
        print *,' *** Invalid value(s) of pa found by valids ***'
        print *,' pa: allowed range = +/- ',patext
        print *,' pa: min, max = ',patmin,patmax
        solnok = .false.
*       Scan for location of relevant extremum
        if ( abs(patmin).ge.patext ) then
          call scan3D (pa, nxpa, nypa, nla, patmin, 'pa minimum')
        endif
        if ( abs(patmax).ge.patext ) then
          call scan3D (pa, nxpa, nypa, nla, patmax, 'pa maximum')
        endif
      endif
**    print *,' qa: min, max = ',qatmin,qatmax
      if ( abs(qatmin).ge.qatext .or. abs(qatmax).ge.qatext ) then
        print *,' '
        print *,' *** Invalid value(s) of qa found by valids ***'
        print *,' qa: allowed range = +/- ',qatext
        print *,' qa: min, max = ',qatmin,qatmax
        solnok = .false.
*       Scan for location of relevant extremum
        if ( abs(qatmin).ge.qatext ) then
          call scan3D (qa, nxpa, nypa, nla, qatmin, 'qa minimum')
        endif
        if ( abs(qatmax).ge.qatext ) then
          call scan3D (qa, nxpa, nypa, nla, qatmax, 'qa maximum')
        endif
      endif
**    print *,' ast: min, max = ',astmin,astmax
      if ( abs(astmin).ge.astext .or. abs(astmax).ge.astext ) then
        print *,' '
        print *,' *** Invalid value(s) of ast found by valids ***'
        print *,' ast: allowed range = +/- ',astext
        print *,' ast: min, max = ',astmin,astmax
        solnok = .false.
*       Scan for location of relevant extremum
        if ( abs(astmin).ge.astext ) then
          call scan2D (ast, nxta, nyta, astmin, 'AST minimum')
        endif
        if ( abs(astmax).ge.astext ) then
          call scan2D (ast, nxta, nyta, astmax, 'AST maximum')
        endif
      endif
**    print *,' wekta: min, max = ',wekmin,wekmax
      if ( abs(wekmin).ge.wtaext .or. abs(wekmax).ge.wtaext ) then
        print *,' '
        print *,' *** Invalid value(s) of wekta found by valids ***'
        print *,' wekta: allowed range = +/- ',wtaext
        print *,' wekta: min, max = ',wekmin,wekmax
        solnok = .false.
*       Scan for location of relevant extremum
        if ( abs(wekmin).ge.wtaext ) then
          call scan2D (wekta, nxta, nyta, wekmin, 'wekta minimum')
        endif
        if ( abs(wekmax).ge.wtaext ) then
          call scan2D (wekta, nxta, nyta, wekmax, 'wekta maximum')
        endif
      endif
**    print *,' tauxa: min, max = ',txamin,txamax
**    print *,' tauya: min, max = ',tyamin,tyamax
      if ( abs(txamin).ge.tauext .or. abs(txamax).ge.tauext .or.
     &     abs(tyamin).ge.tauext .or. abs(tyamax).ge.tauext ) then
        print *,' '
        print *,' *** Invalid value(s) of tauxa',
     &          ' or tauya found by valids ***'
        print *,' tauxa, tauya: allowed range = +/- ',tauext
        print *,' tauxa: min, max = ',txamin,txamax
        print *,' tauya: min, max = ',tyamin,tyamax
        solnok = .false.
*       Scan for location of relevant extremum
        if ( abs(txamin).ge.tauext ) then
          call scan2D (tauxa, nxpa, nypa, txamin, 'tauxa minimum')
        endif
        if ( abs(txamax).ge.tauext ) then
          call scan2D (tauxa, nxpa, nypa, txamax, 'tauxa maximum')
        endif
        if ( abs(tyamin).ge.tauext ) then
          call scan2D (tauya, nxpa, nypa, tyamin, 'tauya minimum')
        endif
        if ( abs(tyamax).ge.tauext ) then
          call scan2D (tauya, nxpa, nypa, tyamax, 'tauya maximum')
        endif
      endif
#endif /* not ocean_only */

#ifndef atmos_only
*     Check oceanic fields for extreme values
*     =======================================
*     Values in the neighbourhood of invalid values
*     in any of the fields scanned will be printed
      pocmin =  bignum
      pocmax = -bignum
      qocmin =  bignum
      qocmax = -bignum
      sstmin =  bignum
      sstmax = -bignum
      wekmin =  bignum
      wekmax = -bignum

!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k)
!$OMP&         SHARED  (pocmin,pocmax,po,qocmin,qocmax,qo,
!$OMP&                  sstmin,sstmax,sst,wekmin,wekmax,wekto)

      do k=1,nlo
!$OMP   DO SCHEDULE (STATIC)
!$OMP&     REDUCTION(min:pocmin,qocmin) REDUCTION(max:pocmax,qocmax)
        do j=1,nypo
          do i=1,nxpo
            pocmin = min( pocmin, po(i,j,k) )
            pocmax = max( pocmax, po(i,j,k) )
            qocmin = min( qocmin, qo(i,j,k) )
            qocmax = max( qocmax, qo(i,j,k) )
          enddo
        enddo
!$OMP   END DO
      enddo

!$OMP DO SCHEDULE (STATIC)
!$OMP&   REDUCTION(min:sstmin,wekmin) REDUCTION(max:sstmax,wekmax)
      do j=1,nyto
        do i=1,nxto
          sstmin = min( sstmin, sst(i,j) )
          sstmax = max( sstmax, sst(i,j) )
          wekmin = min( wekmin, wekto(i,j) )
          wekmax = max( wekmax, wekto(i,j) )
        enddo
      enddo
!$OMP END DO NOWAIT

!$OMP END PARALLEL

**    print *,' po: min, max = ',pocmin,pocmax
      if ( abs(pocmin).ge.pocext .or. abs(pocmax).ge.pocext ) then
        print *,' '
        print *,' *** Invalid value(s) of po found by valids ***'
        print *,' po: allowed range = +/- ',pocext
        print *,' po: min, max = ',pocmin,pocmax
        solnok = .false.
*       Scan for location of relevant extremum
        if ( abs(pocmin).ge.pocext ) then
          call scan3D (po, nxpo, nypo, nlo, pocmin, 'po minimum')
        endif
        if ( abs(pocmax).ge.pocext ) then
          call scan3D (po, nxpo, nypo, nlo, pocmax, 'po maximum')
        endif
      endif
**    print *,' qo: min, max = ',qocmin,qocmax
      if ( abs(qocmin).ge.qocext .or. abs(qocmax).ge.qocext ) then
        print *,' '
        print *,' *** Invalid value(s) of qo found by valids ***'
        print *,' qo: allowed range = +/- ',qocext
        print *,' qo: min, max = ',qocmin,qocmax
        solnok = .false.
*       Scan for location of relevant extremum
        if ( abs(qocmin).ge.qocext ) then
          call scan3D (qo, nxpo, nypo, nlo, qocmin, 'qo minimum')
        endif
        if ( abs(qocmax).ge.qocext ) then
          call scan3D (qo, nxpo, nypo, nlo, qocmax, 'qo maximum')
        endif
      endif
**    print *,' sst: min, max = ',sstmin,sstmax
      if ( abs(sstmin).ge.sstext .or. abs(sstmax).ge.sstext ) then
        print *,' '
        print *,' *** Invalid value(s) of sst found by valids ***'
        print *,' sst: allowed range = +/- ',sstext
        print *,' sst: min, max = ',sstmin,sstmax
        solnok = .false.
*       Scan for location of relevant extremum
        if ( abs(sstmin).ge.sstext ) then
          call scan2D (sst, nxto, nyto, sstmin, 'SST minimum')
        endif
        if ( abs(sstmax).ge.sstext ) then
          call scan2D (sst, nxto, nyto, sstmax, 'SST maximum')
        endif
      endif
**    print *,' wekto: min, max = ',wekmin,wekmax
      if ( abs(wekmin).ge.wtoext .or. abs(wekmax).ge.wtoext ) then
        print *,' '
        print *,' *** Invalid value(s) of wekto found by valids ***'
        print *,' wekto: allowed range = +/- ',wtoext
        print *,' wekto: min, max = ',wekmin,wekmax
        solnok = .false.
*       Scan for location of relevant extremum
        if ( abs(wekmin).ge.wtoext ) then
          call scan2D (wekto, nxto, nyto, wekmin, 'wekto minimum')
        endif
        if ( abs(wekmax).ge.wtoext ) then
          call scan2D (wekto, nxto, nyto, wekmax, 'wekto maximum')
        endif
      endif

*     Check full oceanic layer thicknesses for extreme values
*     =======================================================
*     Top and bottom layers need to be treated differently from
*     intermediate layer(s), so might as well scan for their
*     extrema separately, and thus get more informative output.
*     Notation: t = top, i = intermediate, b = bottom, a = all.
*     hffail notes whether any invalid thicknesses have been found
*     (in which case the worst case position will be located)
      hffail = .false.
      if ( hcheck ) then
        do k=1,nlo-1
          rgpoc(k) = 1.0d0/gpoc(k)
        enddo

*       Scan for extreme values, treating t, i & b cases separately
*       -----------------------------------------------------------
        hfmint =  bignum
        hfmaxt = -bignum
        hfmini =  bignum
        hfmaxi = -bignum
        hfminb =  bignum
        hfmaxb = -bignum
!$OMP PARALLEL DEFAULT (NONE)
!$OMP&         PRIVATE (i,j,k,etaoc,hfull)
!$OMP&         SHARED  (hoc,hfmint,hfmaxt,hfmini,hfmaxi,
!$OMP&                  hfminb,hfmaxb,rgpoc,po,dtopoc)
!$OMP   DO SCHEDULE (STATIC)
!$OMP&  REDUCTION(min:hfmint,hfmini,hfminb)
!$OMP&  REDUCTION(max:hfmaxt,hfmaxi,hfmaxb)
        do j=1,nypo
          do i=1,nxpo
            do k=1,nlo-1
              etaoc(k) = rgpoc(k)*( po(i,j,k+1) - po(i,j,k) )
            enddo
*           Top layer
            hfull = hoc(1) - etaoc(1)
            hfmint = min( hfmint, hfull )
            hfmaxt = max( hfmaxt, hfull )
*           Intermediate layers
            do k=2,nlo-1
              hfull = hoc(k) - etaoc(k) + etaoc(k-1)
              hfmini = min( hfmini, hfull )
              hfmaxi = max( hfmaxi, hfull )
            enddo
*           Bottom layer
            hfull = hoc(nlo) + etaoc(nlo-1) - dtopoc(i,j)
            hfminb = min( hfminb, hfull )
            hfmaxb = max( hfmaxb, hfull )
          enddo
        enddo
!$OMP END DO NOWAIT
!$OMP END PARALLEL
**      print *,' hfull top: min, max = ',hfmint,hfmaxt
**      print *,' hfull int: min, max = ',hfmini,hfmaxi
**      print *,' hfull bot: min, max = ',hfminb,hfmaxb

*       Check if all layers are valid; if not, investigate in more detail
*       -----------------------------------------------------------------
*       Find extrema across all layers
        hfmina = min( hfmint, hfmini, hfminb )
        hfmaxa = max( hfmaxt, hfmaxi, hfmaxb )
**      print *,' hfull all: min, max = ',hfmina,hfmaxa
        do k=1,nlo
          hfbad(k) = 0.0d0
        enddo
        if ( hfmina.le.thkmin ) then
          hffail = .true.
*         Scan for location of the most extreme extremum,
*         and compute fraction of unacceptable values.
          do j=1,nypo
            wtj = 1.0d0
            if ( j.eq.1 .or. j.eq.nypo ) wtj = 0.5d0
            do i=1,nxpo
              wti = 1.0d0
              if ( i.eq.1 .or. i.eq.nxpo ) wti = 0.5d0
              do k=1,nlo-1
                etaoc(k) = rgpoc(k)*( po(i,j,k+1) - po(i,j,k) )
              enddo
*             Top layer
              hfull = hoc(1) - etaoc(1)
              if ( hfull.eq.hfmina ) then
                im = i
                jm = j
                km = 1
              endif
              if ( hfull.lt.thkmin ) hfbad(1) = hfbad(1) + wti*wtj
*             Intermediate layers
              do k=2,nlo-1
                hfull = hoc(k) - etaoc(k) + etaoc(k-1)
                if ( hfull.eq.hfmina ) then
                  im = i
                  jm = j
                  km = k
                endif
                if ( hfull.lt.thkmin ) hfbad(k) = hfbad(k) + wti*wtj
              enddo
*             Bottom layer
              hfull = hoc(nlo) + etaoc(nlo-1) - dtopoc(i,j)
              if ( hfull.eq.hfmina ) then
                im = i
                jm = j
                km = nlo
              endif
              if ( hfull.lt.thkmin ) hfbad(nlo) = hfbad(nlo) + wti*wtj
            enddo
          enddo
*         Compute bad value fractions, and convert to percentage
          do k=1,nlo
            hfbad(k) = 100.0d0*hfbad(k)*ocnorm
**          print *,' k, hfbad = ',k,hfbad(k)
          enddo
        endif

*       Apply termination criterion
*       ---------------------------
        if ( spfail ) then
*         Terminate if any unacceptable values found
          if ( hffail ) then
            solnok = .false.
            print *,' '
            print *,' *** Invalid value(s) of oceanic layer',
     &              ' thickness found by valids ***'
            print *,' min. acceptable layer thickness (m) = ',thkmin
*           Print percentage of invalid values, just for information
            print *,' Layer  %age invalid'
            do k=1,nlo
              write(*,'(i5,f13.3)') k,hfbad(k)
            enddo
          endif
         else
*         Terminate if an unacceptable fraction is found in any layer
          pcfail = .false.
          do k=1,nlo
            if ( hfbad(k).gt.critpc ) pcfail = .true.
          enddo
          if ( pcfail ) then
            solnok = .false.
            print *,' '
            print *,' *** Excessive %age of invalid values of',
     &              ' oceanic layer thickness found by valids ***'
            print *,' min. acceptable layer thickness (m) = ',thkmin
            write(*,'(a,f8.3)') '  %age threshold = ',critpc
            print *,' Layer  %age invalid'
            do k=1,nlo
              write(*,'(i5,f13.3)') k,hfbad(k)
            enddo
          endif
        endif

*       If terminating the run, print out diagnostics
*       ---------------------------------------------
        if ( .not.solnok ) then
          print *,' '
*         If hfull has any invalid values, even if it isn't the
*         cause of the decision to terminate, then compute and
*         print values in neighbourhood of the worst extremum.
*         The idea is that even if bad values of hfull weren't
*         deemed bad enough to cause the program to terminate,
*         they may cause problems in other fields.
*         Only if hffail is true will im, jm and km have been set
*         Only do for termination case, else may get excessive output
*         (e.g. if spfail is false, and the fraction of invalid
*          values is less than criptc, no output is required,
*          even if hfull contains some invalid values).
          if ( hffail ) then
            print *,' hfull top layer : min, max = ',hfmint,hfmaxt
            if ( nlo.gt.2 ) then
              print *,' hfull int layers: min, max = ',hfmini,hfmaxi
            endif
            print *,' hfull bot layer : min, max = ',hfminb,hfmaxb
            print *,' hfull all layers: min, max = ',hfmina,hfmaxa
            print *,' min. acceptable layer thickness (m) = ',thkmin
            write(*,'(2x,a,3i7)')
     &              'hfull minimum located at i, j, k = ',im,jm,km
            write(*,'(2x,a,f12.4)')
     &              'unperturbed thickness hoc(k) = ',hoc(km)
            ilo = max(  1 ,im-3)
            ihi = min(nxpo,im+3)
            jlo = max(  1 ,jm-jwid)
            jhi = min(nypo,jm+jwid)
            do j=jhi,jlo,-1
              joff = j - jlo + 1
              do i=ilo,ihi
                ioff = i - ilo + 1
                do k=1,nlo-1
                  etaarr(ioff,joff,k) =
     &              rgpoc(k)*( po(i,j,k+1) - po(i,j,k) )
                enddo
                if ( km.eq.1 ) then
*                 Top layer
                  hfarr(ioff,joff) = hoc(1) - etaarr(ioff,joff,1)
                 else if ( km.eq.nlo ) then
*                 Bottom layer
                  hfarr(ioff,joff) = hoc(nlo) +
     &              etaarr(ioff,joff,nlo-1) - dtopoc(i,j)
                 else
*                 Intermediate layer
                  hfarr(ioff,joff) = hoc(km) -
     &              etaarr(ioff,joff,km) + etaarr(ioff,joff,km-1)
                endif
              enddo
            enddo
*           Print neighbourhood of extremum
            print *,' '
            write(*,'(2x,a)') 'hfull(k) in neighbourhood of minimum:'
            write(*,'(3x,7i13)') (i,i=ilo,ihi)
            do j=jhi,jlo,-1
              joff = j - jlo + 1
              write(*,'(i7,1p,7d13.5)')
     &                j,(hfarr(i-ilo+1,joff),i=ilo,ihi)
            enddo
            if ( km.gt.1 ) then
              print *,' '
              write(*,'(2x,a)') 'eta(k-1) in neighbourhood of minimum:'
              write(*,'(3x,7i13)') (i,i=ilo,ihi)
              do j=jhi,jlo,-1
                joff = j - jlo + 1
                write(*,'(i7,1p,7d13.5)')
     &                  j,(etaarr(i-ilo+1,joff,km-1),i=ilo,ihi)
              enddo
            endif
            if ( km.lt.nlo ) then
              print *,' '
              write(*,'(2x,a)') 'eta(k) in neighbourhood of minimum:'
              write(*,'(3x,7i13)') (i,i=ilo,ihi)
              do j=jhi,jlo,-1
                joff = j - jlo + 1
                write(*,'(i7,1p,7d13.5)')
     &                  j,(etaarr(i-ilo+1,joff,km),i=ilo,ihi)
              enddo
            endif
            if ( km.eq.nlo ) then
              print *,' '
              write(*,'(2x,a,a)') 'topography in neighbourhood',
     &                            ' of minimum:'
              write(*,'(3x,7i13)') (i,i=ilo,ihi)
              do j=jhi,jlo,-1
                joff = j - jlo + 1
                write(*,'(i7,1p,7d13.5)') j,(dtopoc(i,j),i=ilo,ihi)
              enddo
            endif
            print *,' '
          endif
        endif

      endif
#endif /* not atmos_only */

      END SUBROUTINE valids
c
c***********************************************************************
c
      SUBROUTINE scan2D (array, nx, ny, extrem, string)

*     Scan a 2D array of dimensions (nx,ny) for the location of the
*     value extrem, and print the location and surrounding values.
*     'string' identifies the extremum being located
*
      IMPLICIT NONE
*
*     Subroutine arguments
      integer, INTENT(IN) :: nx,ny
      double precision, INTENT(IN) :: array(nx,ny),extrem
      character (len=*) :: string
*
*     Local parameters
      integer jwid
      parameter ( jwid = 4 )
*
*     Local variables
      integer i,j,im,jm,ilo,ihi
      logical found

      found = .false.
*     Initialise im and jm (to values outside
*     array bounds) to appease gfortran compiler
      im = -1
      jm = -1
      do j=1,ny
        do i=1,nx
          if ( array(i,j).eq.extrem ) then
            write(*,'(2x,a,a,2i7)') string,' located at i, j = ',i,j
            im = i
            jm = j
            found = .true.
          endif
        enddo
      enddo
*     Print neighbourhood of extremum
      if ( found ) then
        ilo = max( 1,im-3)
        ihi = min(nx,im+3)
        write(*,'(3x,7i13)') (i,i=ilo,ihi)
        do j=min(ny,jm+jwid),max(1,jm-jwid),-1
          write(*,'(i7,1p,7d13.5)') j,(array(i,j),i=ilo,ihi)
        enddo
       else
        print *,' '
        print *,' Problem in scan2D; sought value not found'
      endif

      END SUBROUTINE scan2D
c
c***********************************************************************
c
      SUBROUTINE scan3D (array, nx, ny, nl, extrem, string)

*     Scan a 3D array of dimensions (nx,ny,nl) for the location of the
*     value extrem, and print the location and surrounding values.
*     'string' identifies the extremum being located
*
      IMPLICIT NONE
*
*     Subroutine arguments
      integer, INTENT(IN) :: nx,ny,nl
      double precision, INTENT(IN) :: array(nx,ny,nl),extrem
      character (len=*) :: string
*
*     Local parameters
      integer jwid
      parameter ( jwid = 4 )
*
*     Local variables
      integer i,j,k,im,jm,km,ilo,ihi
      logical found

      found = .false.
*     Initialise im, jm and km (to values outside
*     array bounds) to appease gfortran compiler
      im = -1
      jm = -1
      km = -1
      do k=1,nl
        do j=1,ny
          do i=1,nx
            if ( array(i,j,k).eq.extrem ) then
              write(*,'(2x,a,a,3i7)')
     &              string,' located at i, j, k = ',i,j,k
              im = i
              jm = j
              km = k
              found = .true.
            endif
          enddo
        enddo
      enddo
*     Print neighbourhood of extremum
      if ( found ) then
        ilo = max( 1,im-3)
        ihi = min(nx,im+3)
        write(*,'(3x,7i13)') (i,i=ilo,ihi)
        do j=min(ny,jm+jwid),max(1,jm-jwid),-1
          write(*,'(i7,1p,7d13.5)') j,(array(i,j,km),i=ilo,ihi)
        enddo
       else
        print *,' '
        print *,' Problem in scan3D; sought value not found'
      endif

      END SUBROUTINE scan3D
c
c***********************************************************************
c
      END MODULE omuse_valsubs
c
c***********************************************************************
