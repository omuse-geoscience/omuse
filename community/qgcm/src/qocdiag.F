      SUBROUTINE qt4dif_int(qt2dif_,qt4dif_,qotjac_,qotent_,dqdt_)
      
*     Compute oceanic dq/dt and its component terms. The code to do this
*     is borrowed from subroutines qgostep and ocadif in qgosubs.F.
*     An alternative approach would be to be to save the terms in
*     those routines, but this would (a) complicate them, and (b)
*     entail declaring much more storage for transfer. Here the arrays are
*     invoked only when necessary, at the expense of some recomputation.

*     Modules
      USE parameters
#ifndef atmos_only
      USE occonst
      USE ocstate
#endif
      use mondiag

      IMPLICIT NONE

      double precision, intent(out) :: qt2dif_(nlo),qt4dif_(nlo),qotjac_(nlo),
     &                  qotent_(nlo),dqdt_(nlo)

#ifndef atmos_only
!!    Other variables used locally
      integer i,j,k
      double precision adfaco,bcfaco,fohfac(nlo),bdrfac,ah2fac,
     &                 ah4fac,del2p(nxpo,nypo),del4p(nxpo,nypo),
     &                 qt2dif(nxpo,nypo),qt4dif(nxpo,nypo),d6p,
     &                 qotjac(nxpo,nypo),qotent(nxpo,nypo),
     &                 dqdt(nxpo,nypo),wrk(nxpo,nypo),rdto

*     Compute vorticity derivative and components with code from qgosubs.F
*     --------------------------------------------------------------------
      adfaco = 1.0d0/(12.0d0*dxo*dyo*fnot)
*     Version with nondimensional bccooc
      bcfaco = bccooc*dxom2/( 0.5d0*bccooc + 1.0d0 )
      do k=1,nlo
        fohfac(k) = fnot/hoc(k)
      enddo
      bdrfac = 0.5d0*sign(1.0d0,fnot)*delek/hoc(nlo)
      rdto = 1.0d0/dto

      do k=1,nlo

        ah2fac = ah2oc(k)/fnot
        ah4fac = ah4oc(k)/fnot

!$OMP   PARALLEL DEFAULT (NONE)
!$OMP&           PRIVATE (i,j,d6p)
!$OMP&           SHARED  (k,del2p,bcfaco,pom,dxom2,del4p,qt2dif,qt4dif,
!$OMP&                    qotjac,qotent,dqdt,po,qo,fohfac,wekpo,entoc,
!$OMP&                    bdrfac,ah2fac,ah4fac,adfaco,rdto,qom)

*       Compute Del-sqd(p) at previous time level for dissipation
*       ---------------------------------------------------------
*       N & S boundaries (including corners) - mixed BCs
*       No x-derivative term because p is constant along boundary
!$OMP   DO SCHEDULE (STATIC)
        do i=1,nxpo
          del2p(i,  1 ) = bcfaco*( pom(i,   2  ,k) - pom(i,  1 ,k) )
          del2p(i,nypo) = bcfaco*( pom(i,nypo-1,k) - pom(i,nypo,k) )
        enddo
!$OMP   END DO NOWAIT

!$OMP   DO SCHEDULE (STATIC)
        do j=2,nypo-1
*         Western boundary (i=1)
#  ifdef cyclic_ocean
*         Zonally cyclic ocean
          del2p(1,j) = (  pom(1,j-1,k) + pom(nxpo-1,j,k) + pom(2,j,k)
     &                  + pom(1,j+1,k) - 4.0d0*pom(1,j,k) )*dxom2
#  else
*         Finite box ocean - mixed BCs
          del2p( 1 ,j) = bcfaco*( pom(  2  ,j,k) - pom( 1 ,j,k) )
#  endif /* cyclic_ocean */
*         Inner points (standard case)
          do i=2,nxpo-1
            del2p(i,j) = (  pom(i,j-1,k) + pom(i-1,j,k) + pom(i+1,j,k)
     &                    + pom(i,j+1,k) - 4.0d0*pom(i,j,k) )*dxom2
          enddo
*         Eastern boundary (i=nxpo)
#  ifdef cyclic_ocean
*         Zonally cyclic ocean
          del2p(nxpo,j) = del2p(1,j)
#  else
*         Finite box ocean - mixed BCs
          del2p(nxpo,j) = bcfaco*( pom(nxpo-1,j,k) - pom(nxpo,j,k) )
#  endif /* cyclic_ocean */
        enddo
!$OMP   END DO

*       Compute Del-4th(p) at previous time level for dissipation
*       ---------------------------------------------------------
*       N & S boundaries (including corners) - mixed BCs
*       No x-derivative term because p is constant along boundary
!$OMP   DO SCHEDULE (STATIC)
        do i=1,nxpo
          del4p(i,  1 ) = bcfaco*( del2p(i,   2  ) - del2p(i,  1 ) )
          del4p(i,nypo) = bcfaco*( del2p(i,nypo-1) - del2p(i,nypo) )
        enddo
!$OMP   END DO NOWAIT

!$OMP   DO SCHEDULE (STATIC)
        do j=2,nypo-1
*         Western boundary (i=1)
#  ifdef cyclic_ocean
*         Zonally cyclic ocean
          del4p(1,j) = (  del2p(1,j-1) + del2p(nxpo-1,j) + del2p(2,j)
     &                  + del2p(1,j+1) - 4.0d0*del2p(1,j) )*dxom2
#  else
*         Finite box ocean - mixed BCs
          del4p( 1 ,j) = bcfaco*( del2p(  2  ,j) - del2p( 1 ,j) )
#  endif /* cyclic_ocean */
*         Inner points (standard case)
          do i=2,nxpo-1
            del4p(i,j) = (  del2p(i,j-1) + del2p(i-1,j) + del2p(i+1,j)
     &                    + del2p(i,j+1) - 4.0d0*del2p(i,j) )*dxom2
          enddo
*         Eastern boundary (i=nxpo)
#  ifdef cyclic_ocean
*         Zonally cyclic ocean
          del4p(nxpo,j) = del4p(1,j)
#  else
*         Finite box ocean - mixed BCs
          del4p(nxpo,j) = bcfaco*( del2p(nxpo-1,j) - del2p(nxpo,j) )
#  endif /* cyclic_ocean */
        enddo
!$OMP   END DO

*       Zero all the results arrays as a precaution
*       -------------------------------------------
*       This can be removed when the boundary terms are computed
!$OMP   DO SCHEDULE (STATIC)
        do j=1,nypo
          do i=1,nxpo
            qt2dif(i,j) = 0.0d0
            qt4dif(i,j) = 0.0d0
            qotjac(i,j) = 0.0d0
            qotent(i,j) = 0.0d0
            dqdt(i,j)   = 0.0d0
          enddo
        enddo
!$OMP   END DO

*       Compute advective (Jacobian) and diffusive contributions to dq/dt
*       -----------------------------------------------------------------
*       Compute Del-6th(p) except at N & S boundaries
*       No need to apply N & S boundary conditions
!$OMP   DO SCHEDULE (STATIC)
        do j=2,nypo-1
*         Western boundary (i=1)
#  ifdef cyclic_ocean
*         Zonally cyclic ocean
          d6p = dxom2*(  del4p(1,j-1) + del4p(nxpo-1,j) + del4p(2,j)
     &                 + del4p(1,j+1) - 4.0d0*del4p(1,j) )
          qt2dif(1,j) = ah2fac*del4p(1,j)
          qt4dif(1,j) = -ah4fac*d6p
          qotjac(1,j) = adfaco*
     & ( (qo( 2 , j ,k)-qo(nxpo-1, j ,k))*(po( 1 ,j+1,k)-po(  1 ,j-1,k))
     &  +(qo( 1 ,j-1,k)-qo(  1 ,j+1,k))*(po( 2 , j ,k)-po(nxpo-1, j ,k))
     &  + qo(   2  , j ,k)*( po(   2  ,j+1,k) - po(   2  ,j-1,k) )
     &  - qo(nxpo-1, j ,k)*( po(nxpo-1,j+1,k) - po(nxpo-1,j-1,k) )
     &  - qo(   1  ,j+1,k)*( po(   2  ,j+1,k) - po(nxpo-1,j+1,k) )
     &  + qo(   1  ,j-1,k)*( po(   2  ,j-1,k) - po(nxpo-1,j-1,k) )
     &  + po(   1  ,j+1,k)*( qo(   2  ,j+1,k) - qo(nxpo-1,j+1,k) )
     &  - po(   1  ,j-1,k)*( qo(   2  ,j-1,k) - qo(nxpo-1,j-1,k) )
     &  - po(   2  , j ,k)*( qo(   2  ,j+1,k) - qo(   2  ,j-1,k) )
     &  + po(nxpo-1, j ,k)*( qo(nxpo-1,j+1,k) - qo(nxpo-1,j-1,k) ) )
*         (layer-dependent) entrainment/forcing term
          if ( k.eq.1 ) then
!            qotent(1,j) = fohfac(1)*( wekpo(1,j) - entoc(1,j) )
           else if ( k.eq.2 ) then
            qotent(1,j) = fohfac(2)*entoc(1,j)
           else
            qotent(1,j) = 0.0d0
          endif
*         Add bottom drag contribution
          if ( k.eq.nlo ) then
            qotent(1,j) = qotent(1,j) - bdrfac*del2p(1,j)
          endif
          dqdt(1,j) =  qotjac(1,j) + qt2dif(1,j)
     &               + qt4dif(1,j) + qotent(1,j)
#  else
*         Finite box ocean
          qt2dif(1,j) = 0.0d0
          qt4dif(1,j) = 0.0d0
          qotjac(1,j) = 0.0d0
          qotent(1,j) = 0.0d0
          dqdt(1,j)   = rdto*( qo(1,j,k) - qom(1,j,k) )
#  endif /* cyclic_ocean */
*         Inner points (standard case)
          do i=2,nxpo-1
            d6p = dxom2*(  del4p(i,j-1) + del4p(i-1,j) + del4p(i+1,j)
     &                   + del4p(i,j+1) - 4.0d0*del4p(i,j) )
            qt2dif(i,j) = ah2fac*del4p(i,j)
            qt4dif(i,j) = -ah4fac*d6p
            qotjac(i,j) = adfaco*
     &   ( (qo(i+1, j ,k)-qo(i-1, j ,k))*(po( i ,j+1,k)-po( i ,j-1,k))
     &    +(qo( i ,j-1,k)-qo( i ,j+1,k))*(po(i+1, j ,k)-po(i-1, j ,k))
     &    + qo(i+1, j ,k)*( po(i+1,j+1,k) - po(i+1,j-1,k) )
     &    - qo(i-1, j ,k)*( po(i-1,j+1,k) - po(i-1,j-1,k) )
     &    - qo( i ,j+1,k)*( po(i+1,j+1,k) - po(i-1,j+1,k) )
     &    + qo( i ,j-1,k)*( po(i+1,j-1,k) - po(i-1,j-1,k) )
     &    + po( i ,j+1,k)*( qo(i+1,j+1,k) - qo(i-1,j+1,k) )
     &    - po( i ,j-1,k)*( qo(i+1,j-1,k) - qo(i-1,j-1,k) )
     &    - po(i+1, j ,k)*( qo(i+1,j+1,k) - qo(i+1,j-1,k) )
     &    + po(i-1, j ,k)*( qo(i-1,j+1,k) - qo(i-1,j-1,k) ) )
*           (layer-dependent) entrainment/forcing term
            if ( k.eq.1 ) then
!              qotent(i,j) = fohfac(1)*( wekpo(i,j) - entoc(i,j) )
             else if ( k.eq.2 ) then
              qotent(i,j) = fohfac(2)*entoc(i,j)
             else
              qotent(i,j) = 0.0d0
            endif
*           Add bottom drag contribution
            if ( k.eq.nlo ) then
              qotent(i,j) = qotent(i,j) - bdrfac*del2p(i,j)
            endif
            dqdt(i,j) =  qotjac(i,j) + qt2dif(i,j)
     &                 + qt4dif(i,j) + qotent(i,j)
          enddo
*         Eastern boundary (i=nxpo)
#  ifdef cyclic_ocean
*         Zonally cyclic ocean
          qt2dif(nxpo,j) = qt2dif(1,j)
          qt4dif(nxpo,j) = qt4dif(1,j)
          qotjac(nxpo,j) = qotjac(1,j)
          qotent(nxpo,j) = qotent(1,j)
          dqdt(nxpo,j)   = dqdt(1,j)
#  else
*         Finite box ocean
          qt2dif(nxpo,j) = 0.0d0
          qt4dif(nxpo,j) = 0.0d0
          qotjac(nxpo,j) = 0.0d0
          qotent(nxpo,j) = 0.0d0
          dqdt(nxpo,j)   = rdto*( qo(nxpo,j,k) - qom(nxpo,j,k) )
#  endif /* cyclic_ocean */
        enddo
!$OMP   END DO

*       Add zonal boundary contributions (to dqdt only), by differencing
*       qo in time. This is slightly displaced in time, but only by dto/2
!$OMP   DO SCHEDULE (STATIC)
        do i=1,nxpo
          dqdt(i,  1 ) = rdto*( qo(i,  1 ,k) - qom(i,  1 ,k) )
          dqdt(i,nypo) = rdto*( qo(i,nypo,k) - qom(i,nypo,k) )
        enddo
!$OMP   END DO

!$OMP   END PARALLEL

        wrk=po(:,:,k)*qt2dif
        call genint(wrk, nxpo,nypo, 0.5d0,0.5d0, qt2dif_(k))
        wrk=po(:,:,k)*qt4dif
        call genint(wrk, nxpo,nypo, 0.5d0,0.5d0, qt4dif_(k))
        wrk=po(:,:,k)*qotjac
        call genint(wrk, nxpo,nypo, 0.5d0,0.5d0, qotjac_(k))
        wrk=po(:,:,k)*qotent
        call genint(wrk, nxpo,nypo, 0.5d0,0.5d0, qotent_(k))
        wrk=po(:,:,k)*dqdt
        call genint(wrk, nxpo,nypo, 0.5d0,0.5d0, dqdt_(k))

      enddo

#endif /* not atmos_only */

      END SUBROUTINE qt4dif_int
