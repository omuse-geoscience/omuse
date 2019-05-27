!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module domain_size

!BOP
! !MODULE: domain_size
!
! !DESCRIPTION:
!  This module contains parameters for the global model domain size
!  decomposition block size.  It is used by the domain and block
!  modules for decomposing the model domain across processors.
!  Variables are now set in POP\_DomainSizeMod and this routine
!  only provide back compatibility.
!
! !REVISION HISTORY:
!  CVS:$Id: domain_size.F90.test,v 1.2 2003/03/25 13:41:48 pwjones Exp $
!  CVS:$Name: POP_2_0_1 $

! !USES:

   use kinds_mod
   use POP_DomainSizeMod

   implicit none
   private
   save

! !DEFINED PARAMETERS:

   integer (int_kind), parameter, public ::  &  ! model size parameters
      nx_global = POP_nxGlobal, &
      ny_global = POP_nyGlobal, &
      km = POP_km,              &
      nt = POP_nt

   integer (int_kind), parameter, public :: &
      block_size_x = POP_blockSizeX, &
      block_size_y = POP_blockSizeY

   integer (int_kind), parameter, public :: &
      max_blocks_clinic = POP_maxBlocksClinic, &
      max_blocks_tropic = POP_maxBlocksTropic

!EOP
!BOC
!EOC
!***********************************************************************

 end module domain_size

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
