PROGRAM collectdata2003
#ifdef USE_MPI
  USE yaxt, ONLY: xt_initialize, xt_finalize, xt_idxlist, xt_idxstripes_new, &
       xt_idxlist_delete, xi => xt_int_kind, xt_stripe
  USE iso_c_binding, ONLY: c_int
#endif

  IMPLICIT NONE

  INCLUDE 'cdi.inc'
  INCLUDE 'cdipio.inc'

#ifdef USE_MPI
  INCLUDE 'mpif.h'

  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(9)
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND(14)
#endif

  ! For parallel IO:
  !  Parameter and variables needed.
  INTEGER, PARAMETER :: nProcsIO     = 3
  ! INTEGER, PARAMETER :: IOMode       = PIO_NONE
  ! INTEGER, PARAMETER :: IOMode       = PIO_MPI
  ! INTEGER, PARAMETER :: IOMode       = PIO_WRITER
  ! INTEGER, PARAMETER :: IOMode       = PIO_ASYNCH
  INTEGER, PARAMETER :: IOMode       = PIO_FPGUARD

  INTEGER :: commModel
#ifdef USE_MPI
  LOGICAL :: run_model
  INTEGER :: commGlob, ierror, pio_namespace
#else
  LOGICAL, PARAMETER :: run_model = .TRUE.
#endif

  ! Start parallel environment
#ifdef USE_MPI
  CALL MPI_INIT ( ierror )
  commGlob = MPI_COMM_WORLD
  CALL xt_initialize(commGlob)

  ! For parallel IO:
  ! Initialize environment.
  commModel = pioInit(commGlob, nProcsIO, IOMode, pio_namespace, 1.1, &
       cdiPioNoPostCommSetup)
  run_model = commModel /= MPI_COMM_NULL
  IF (run_model) CALL namespaceSetActive(pio_namespace)
#else
  commModel = 0
#endif

  IF (run_model) CALL modelrun ( commModel )

#ifdef USE_MPI
  ! For parallel IO:
  ! Cleanup environment.
  IF (run_model) CALL pioFinalize ()
  CALL xt_finalize
  CALL MPI_FINALIZE ( ierror )
#endif

CONTAINS

  SUBROUTINE modelrun ( commModel )

    INTEGER, INTENT(in) :: commModel

    INTEGER, PARAMETER :: nlon     = 12   ! Number of longitudes
    INTEGER, PARAMETER :: nlat     =  6   ! Number of latitudes
    INTEGER, PARAMETER :: nlev     =  5   ! Number of levels
    INTEGER, PARAMETER :: nts      =  3   ! Number of time steps
    INTEGER, PARAMETER :: vdate    = 19850101
    INTEGER, PARAMETER :: vtime    = 120000
    INTEGER, PARAMETER :: filetype = FILETYPE_GRB

    INTEGER :: gridID, zaxisID1, zaxisID2, taxisID
    INTEGER :: vlistID, varID1, varID2, streamID, tsID
    INTEGER :: i, nmiss, status
    DOUBLE PRECISION :: lons ( nlon ), lats ( nlat ), levs ( nlev )
    DOUBLE PRECISION :: var1 ( nlon * nlat ), var2 ( nlon * nlat * nlev )
    CHARACTER(len=256) :: varname
    INTEGER :: last, start
#ifdef USE_MPI
    INTEGER :: rank, comm_size, ierror, chunk
    TYPE var1ddeco
      INTEGER :: start, chunksize
      TYPE(xt_idxlist) :: partdesc
    END TYPE var1ddeco
    TYPE(var1ddeco) :: vardeco1, vardeco2
#endif

    lons = (/0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330/)
    lats = (/-75, -45, -15, 15, 45, 75/)
    levs = (/101300, 92500, 85000, 50000, 20000/)

    nmiss = 0

    !     Create a regular lon/lat grid
    gridID = gridCreate(GRID_LONLAT, nlon*nlat)
    CALL gridDefXsize(gridID, nlon)
    CALL gridDefYsize(gridID, nlat)
    CALL gridDefXvals(gridID, lons)
    CALL gridDefYvals(gridID, lats)

    !     Create a surface level Z-axis
    zaxisID1 = zaxisCreate(ZAXIS_SURFACE, 1)

    !     Create a pressure level Z-axis
    zaxisID2 = zaxisCreate(ZAXIS_PRESSURE, nlev)
    CALL zaxisDefLevels(zaxisID2, levs)

    !     Create a variable list
    vlistID = vlistCreate()

    !     Define the variables
    varID1 = vlistDefVar(vlistID, gridID, zaxisID1, TIME_VARIABLE)
    varID2 = vlistDefVar(vlistID, gridID, zaxisID2, TIME_VARIABLE)

#ifdef USE_MPI
    CALL mpi_comm_rank(commModel, rank, ierror)
    IF (ierror /= mpi_success) STOP 1
    CALL mpi_comm_size(commModel, comm_size, ierror)
    IF (ierror /= mpi_success) STOP 1
    start = uniform_partition_start((/ 1, SIZE(var1) /), comm_size, rank + 1)
    chunk = uniform_partition_start((/ 1, SIZE(var1) /), comm_size, rank + 2) &
         - start
    vardeco1 = var1ddeco(start, chunk, &
         xt_idxstripes_new(xt_stripe(INT(start, xi) - 1_xi, 1_xi, &
         INT(chunk, c_int))))
    start = uniform_partition_start((/ 1, SIZE(var2) /), comm_size, rank + 1)
    chunk = uniform_partition_start((/ 1, SIZE(var2) /), comm_size, rank + 2) &
         - start
    vardeco2 = var1ddeco(start, chunk, &
         xt_idxstripes_new(xt_stripe(INT(start, xi) - 1_xi, 1_xi, &
         INT(chunk, c_int))))
#endif
    !     Define the variable names
    varname = 'varname1'
    CALL vlistDefVarName(vlistID, varID1, varname)
    varname = 'varname2'
    CALL vlistDefVarName(vlistID, varID2, varname)

    !     Create a Time axis
    taxisID = taxisCreate(TAXIS_ABSOLUTE)

    !     Assign the Time axis to the variable list
    CALL vlistDefTaxis(vlistID, taxisID)

    streamID = streamOpenWrite('example.grb', filetype)
    IF ( streamID < 0 ) THEN
      WRITE(0,*) cdiStringError(streamID)
      STOP
    END IF

    !     Assign the variable list to the dataset
    CALL streamDefVlist(streamID, vlistID)

#ifdef USE_MPI
    !     For parallel IO:
    !          End definition stage for CDI resources,
    !          balance load of variables among collecting IO server,
    !          Decompose data on model processes for IO
    !          Transfer resources to the collecting IO server.
    CALL pioEndDef ();
#endif
    !     Loop over the number of time steps
    DO tsID = 0, nts-1
      CALL taxisDefVdate ( taxisID, vdate + tsID )
      CALL taxisDefVtime ( taxisID, vtime )
      status = streamDefTimestep ( streamID, tsID )

      ! For parallel IO:
      ! Inquire start index and chunk for IO transposition of var1
#ifdef USE_MPI
      start = vardeco1%start
      last = start + vardeco1%chunksize - 1
#else
      start = 1
      last = SIZE(var1)
#endif
      ! Init decomposed data for var1
      DO i = start, last
        var1(i) = 1.1
      END DO
      ! Write var1
#ifdef USE_MPI
      CALL streamwritevarpart(streamID, varID1, var1(start:last), nmiss, &
           vardeco1%partdesc)
#else
      CALL streamWriteVar(streamID, varID1, var1, nmiss)
#endif
      ! For parallel IO:
      ! Inquire start index and chunk for IO transposition of var2
#ifdef USE_MPI
      start = vardeco2%start
      last = start + vardeco2%chunksize - 1
#else
      start = 1
      last = SIZE(var2)
#endif
      ! Init decomposed data for var2
      DO i = start, last
        var2(i) = 2.2
      END DO
      ! Write var2
#ifdef USE_MPI
      CALL streamwritevarpart(streamID, varID2, var2(start:last), nmiss, &
           vardeco2%partdesc)
#else
      CALL streamWriteVar(streamID, varID2, var2, nmiss)
#endif

#ifdef USE_MPI
      ! For parallel IO:
      ! Start transmission of all data for output in this timestep to IO server.
      CALL pioWriteTimestep
#endif
    END DO

#ifdef USE_MPI
    ! For parallel IO:
    ! Preparation for local cleanup
    CALL pioEndTimestepping ()
#endif

    !     Close the output stream
    CALL streamClose(streamID)

    !     Destroy the objects
    CALL vlistDestroy(vlistID)
    CALL taxisDestroy(taxisID)
    CALL zaxisDestroy(zaxisID1)
    CALL zaxisDestroy(zaxisID2)
    CALL gridDestroy(gridID)
#ifdef USE_MPI
    CALL xt_idxlist_delete(vardeco1%partdesc)
    CALL xt_idxlist_delete(vardeco2%partdesc)
    CALL mpi_barrier(commModel, ierror)
    IF (ierror /= mpi_success) STOP 1
#endif
  END SUBROUTINE modelrun

#ifdef USE_MPI
  FUNCTION uniform_partition_start(set_interval, nparts, part_idx) &
       RESULT(start)
    INTEGER(i4), INTENT(in) :: nparts
    INTEGER(i4), INTENT(in) :: set_interval(2)
    INTEGER(i4), INTENT(in) :: part_idx

    INTEGER(i4) :: start, part_offset

    part_offset = INT((INT(set_interval(2) - set_interval(1) + 1, i8) &
         &             * INT(part_idx - 1, i8)) / INT(nparts, i8))
    start = set_interval(1) + part_offset
  END FUNCTION uniform_partition_start
#endif
END PROGRAM collectdata2003
