import sys

template= \
"""
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

 module POP_DomainSizeMod

!BOP
! !MODULE: POP_DomainSizeMod
!
! !DESCRIPTION:
!  This module contains parameters for the global model domain size
!  decomposition block size.  It is used by the domain and block
!  modules for decomposing the model domain across processors.
!
! !REVISION HISTORY:
!  SVN:$Id: POP_DomainSizeMod.F90.test 12 2006-08-15 19:57:39Z  $
!  2006-08-14: Phil Jones
!              New domain size module following new naming conventions

! !USES:

   use POP_KindsMod

   implicit none
   private
   save

! !DEFINED PARAMETERS:

   integer (POP_i4), parameter, public ::  &  ! model size parameters
      POP_nxGlobal =  {POP_nxGlobal} ,&       ! extent of horizontal axis in i direction
      POP_nyGlobal =  {POP_nyGlobal} ,&       ! extent of horizontal axis in j direction
      POP_km = {POP_km}          ,&           ! number of vertical levels
      POP_nt =  2                             ! total number of tracers

   integer (POP_i4), parameter, public :: &
      POP_blockSizeX = {POP_blockSizeX}, &! size of block in first  horizontal dimension
      POP_blockSizeY = {POP_blockSizeY}   ! size of block in second horizontal dimension

   !*** The model will inform the user of the correct
   !*** values for the parameters below.  A value higher than
   !*** necessary will not cause the code to fail, but will
   !*** allocate more memory than is necessary.  A value that
   !*** is too low will cause the code to exit.  
   !*** A good initial guess is found using
   !*** max=(nx_global/block_size_x)*(ny_global/block_size_y)/
   !***         num_procs
 
   integer (POP_i4), parameter, public :: &
      POP_maxBlocksClinic = {POP_maxBlocksClinic},  &! max number of blocks per processor
      POP_maxBlocksTropic = {POP_maxBlocksTropic}    !   in each distribution

!EOP
!BOC
!EOC
!***********************************************************************

 end module POP_DomainSizeMod

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
"""

def optimal_blocksize(x):
     if x<60:
        return x//2
     if x<180: 
        return x//3
     candidates=[ 64, 56, 48, 40, 32, 24, 20]
     remainder=[ x%c for c in candidates]
     m=min(remainder)
     index=remainder.index(m)
     return candidates[index]

if __name__=="__main__":
    """ outputs a POP_DomainSizeMod.F90, given input "kxlxm" """
    try:
      sizes=sys.argv[1].split("x")
      assert len(sizes)==3
      sizes=[int(x) for x in sizes]
    except IndexError as ex:
      print("provide size string ixjxk")
      raise ex
    except Exception as ex:
      print("input error")
      raise ex
    
    POP_nxGlobal,POP_nyGlobal, POP_km=sizes
    POP_blockSizeX=optimal_blocksize(POP_nxGlobal)
    POP_blockSizeY=optimal_blocksize(POP_nyGlobal)
  
    nblocks=8
    print("min number of procs:", (POP_nxGlobal//POP_blockSizeX)*(POP_nyGlobal//POP_blockSizeY)//nblocks)
    print("max number of procs:", (POP_nxGlobal//POP_blockSizeX)*(POP_nyGlobal//POP_blockSizeY))

    print("blocksize x,y and # blocks: ",POP_blockSizeX,POP_blockSizeY, nblocks)

    f=open("POP_DomainSizeMod.F90","w")
    f.write(template.format(POP_nxGlobal=POP_nxGlobal,POP_nyGlobal=POP_nyGlobal, 
       POP_km=POP_km,POP_blockSizeX=POP_blockSizeX,POP_blockSizeY=POP_blockSizeY,
       POP_maxBlocksClinic=nblocks,POP_maxBlocksTropic=nblocks))
    f.close()
