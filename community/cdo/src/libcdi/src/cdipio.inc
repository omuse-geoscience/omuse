! This file was automatically generated, don't edit!
!
! Fortran interface for CDI library version 1.7.0
!
! Author:
! -------
! Uwe Schulzweida, MPI-MET, Hamburg,   June 2015
!

!
!  parallel IO IOMode
!
      INTEGER    PIO_NONE              
      PARAMETER (PIO_NONE               =  0)
      INTEGER    PIO_MPI               
      PARAMETER (PIO_MPI                =  1)
      INTEGER    PIO_WRITER            
      PARAMETER (PIO_WRITER             =  2)
      INTEGER    PIO_ASYNCH            
      PARAMETER (PIO_ASYNCH             =  3)
      INTEGER    PIO_FPGUARD           
      PARAMETER (PIO_FPGUARD            =  4)
!
!  parallel IO routines
!
!                     pioEndDef
      EXTERNAL        pioEndDef

!                     pioEndTimestepping
      EXTERNAL        pioEndTimestepping

!                     pioFinalize
      EXTERNAL        pioFinalize

!                     cdiPioNoPostCommSetup
      EXTERNAL        cdiPioNoPostCommSetup

      INTEGER         pioInit
!                                    (INTEGER         commSuper,
!                                     INTEGER         nProcsIO,
!                                     INTEGER         IOMode,
!                                     INTEGER         pioNamespace,
!                                     REAL            partInflate,
!                                     PROCEDURE       postCommSetupActions)
      EXTERNAL        pioInit

!                     pioWriteTimestep
      EXTERNAL        pioWriteTimestep

!                     cdiPioRDMAProgress
      EXTERNAL        cdiPioRDMAProgress

!                     streamWriteVarPart
!                                    (INTEGER         streamID,
!                                     INTEGER         varID,
!                                     CHOICE          data,
!                                     INTEGER         nmiss,
!                                     TYPE(XT_IDXLIST)partDesc)
      EXTERNAL        streamWriteVarPart

!                     streamWriteScatteredVarPart
!                                    (INTEGER         streamID,
!                                     INTEGER         varID,
!                                     CHOICE          data,
!                                     INTEGER         numBlocks,
!                                     INTEGER         blocklengths,
!                                     INTEGER         displacements,
!                                     INTEGER         nmiss,
!                                     TYPE(XT_IDXLIST)partDesc)
      EXTERNAL        streamWriteScatteredVarPart

