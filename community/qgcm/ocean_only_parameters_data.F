c***********************************************************************
c     Q-GCM Version 1.5.0 : last modified 02/07/2013
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
      MODULE parameters

*     Contains all array dimensioning parameters, and the physical
*     parameters f0 and beta which determine the latitude of the domain
*     All other physical parameters are read in from the file input.params

      IMPLICIT NONE

      PUBLIC
      SAVE

*     The parameters below are for the North Atlantic configuration
*     Atmosphere = 30720 km x 7680 km; Ocean = 4800 km x 4800 km.

*     Array dimensioning parameters and position of ocean in atmosphere
*     -----------------------------------------------------------------

*     User settable grid parameters:
      integer :: nxta,nyta,nla,nxaooc,nyaooc,ndxr,nlo
*     Atmosphere
      PARAMETER ( nxta =  384, nyta =  96, nla = 3 )  !! 80 km atmosphere
!!    PARAMETER ( nxta =  512, nyta = 128, nla = 3 )  !! 60 km atmosphere
!!    PARAMETER ( nxta =  768, nyta = 192, nla = 3 )  !! 40 km atmosphere
!!    PARAMETER ( nxta = 1536, nyta = 384, nla = 3 )  !! 20 km atmosphere
*     Ocean
      PARAMETER ( nxaooc =  60, nyaooc =  60, ndxr = 16, nlo = 3 )  !! 80 km atm.
!!    PARAMETER ( nxaooc =  80, nyaooc =  80, ndxr = 12, nlo = 3 )  !! 60 km atm.
!!    PARAMETER ( nxaooc = 120, nyaooc = 120, ndxr =  8, nlo = 3 )  !! 40 km atm.
!!    PARAMETER ( nxaooc = 240, nyaooc = 240, ndxr =  4, nlo = 3 )  !! 20 km atm.

*     nxta, nyta are the numbers of atmospheric gridcells W/E and N/S.
*     nxaooc, nyaooc are the numbers of atmospheric gridcells over
*     the ocean, and thus define the size of the ocean domain.
*     nla, nlo are the numbers of layers in the atmosphere and ocean.
*     ndxr is the ratio of atmospheric/oceanic gridlengths (dxa/dxo).
*     Need to choose nxta to be at least nxaooc;
*     in the cyclic_ocean case nxaooc must equal nxta.
*     Need to choose nyta to be at least nyaooc

*     Need to choose nxta to be an efficient (factorisable) FFT length.
*     Need to choose nxaooc and ndxr so as to make nxto a good FFT length
*     (see below). This means they both need to have small prime factors.

*     For FFTPACK, optimised routines are provided for transforms
*     of lengths 2, 3, 4 and 5. So nxta, ndxr and nxaooc should
*     ideally be chosen to have prime factors of 2, 3 and 5 only.

*     Derived grid parameters (do not alter):
      integer :: nxpa,nypa,nxtaor,nytaor,nxpaor,nypaor,
     &           nxto,nyto,nxpo,nypo,nx1,ny1
      double precision :: atnorm,ocnorm
      PARAMETER ( nxpa = nxta+1, nypa = nyta+1,
     &            nxtaor = nxta*ndxr, nytaor = nyta*ndxr,
     &            nxpaor = nxtaor+1, nypaor = nytaor+1,
     &            nxto = ndxr*nxaooc, nyto = ndxr*nyaooc,
     &            nxpo = nxto+1, nypo = nyto+1,
     &            nx1 = 1 + (nxta-nxaooc)/2, ny1 = 1 + (nyta-nyaooc)/2,
     &            atnorm = 1.0d0/(nxta*nyta),
     &            ocnorm = 1.0d0/(nxto*nyto) )

*     nxpa, nypa are the numbers of atmospheric p points W/E and N/S.
*     nxtaor, nytaor are the numbers of atmos. T points at ocean resolution
*     nxpaor, nypaor are the numbers of atmos. p points at ocean resolution
*     (required for the new windstress formulation introduced at v1.4.0)
*     nxto, nyto are the numbers of oceanic T points W/E and N/S.
*     nxpo, nypo are the numbers of oceanic p points W/E and N/S.
*     nx1, ny1 are starting indices for the ocean in the atmos. grid.
*     We choose them to centre the ocean in the atmospheric domain (if possible).
*     atnorm, ocnorm are normalisation factors for computing mean
*     values, and are both equal to 1/(number of T gridcells).

*     Rotation parameters (for both atmosphere and ocean)
*     ---------------------------------------------------
      double precision :: fnot,beta
      PARAMETER ( fnot = 9.37456D-05, beta = 1.75360D-11 )  !! 40 North

*     fnot is the Coriolis parameter at the
*     central latitude of the domain (rad s^-1)
*     beta is the y-derivative of the Coriolis parameter (rad s^-1 m^-1)

*     Covariance subsampling and array dimensioning parameters
*     --------------------------------------------------------
      integer :: nscvat,nvcvat,nmcvat,nscvoc,nvcvoc,nmcvoc
      PARAMETER ( nscvat =  2, nvcvat = nxta*nyta/(nscvat*nscvat),
     &            nscvoc = 16, nvcvoc = nxto*nyto/(nscvoc*nscvoc),
     &            nmcvat = nvcvat*(nvcvat+1)/2,
     &            nmcvoc = nvcvoc*(nvcvoc+1)/2 )

*     nscvat is the number of atmos. gridpoints (in either
*     direction) averaged over to produce the subsampled
*     results accumulated to compute covariance statistics
*     nscvoc is corresponding no. of ocean gridpoints used
*     nvcvat is the length of the subsampled atmos. state vector
*     nscvat needs to be a factor of both nxta and nyta
*     nscvoc needs to be a factor of both nxto and nyto
*     nvcvoc is the length of the subsampled ocean state vector
*     nmcvat is the number of independent elements
*     in the covariance matrix for the atmosphere
*     (only the lower triangular elements are stored)
*     nmcvoc is the number of independent elements for the ocean
*     These quantities are used for dimensioning storage
*     in the diagnostics module covaria_prog.F
*     N.B. need to ensure that (especially) nmcvoc doesn't
*     exceed the maximum valid size for an integer*4 number

      END MODULE parameters
c
c***********************************************************************
