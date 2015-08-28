      PROGRAM CDIWRITEF2003

      USE iso_c_binding
      USE mo_cdi

      IMPLICIT NONE

      INTEGER nlon, nlat, nlev, nts
      PARAMETER (nlon = 12)   ! Number of longitudes
      PARAMETER (nlat =  6)   ! Number of latitudes
      PARAMETER (nlev =  5)   ! Number of levels
      PARAMETER (nts  =  3)   ! Number of time steps

      INTEGER gridID, zaxisID1, zaxisID2, taxisID
      INTEGER vlistID, varID1, varID2, streamID, tsID
      INTEGER i, nmiss, status
      DOUBLE PRECISION lons(nlon), lats(nlat), levs(nlev)
      DOUBLE PRECISION var1(nlon*nlat), var2(nlon*nlat*nlev)
      CHARACTER(len=256, kind=c_char) :: varname
      CHARACTER(kind=c_char,len=1), POINTER :: msg(:)

      DATA lons /0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330/
      DATA lats /-75, -45, -15, 15, 45, 75/
      DATA levs /101300, 92500, 85000, 50000, 20000/

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

!     Define the variable names
      varname = "varname1" // c_null_char
      CALL vlistDefVarName(vlistID, varID1, varname)
      CALL vlistDefVarName(vlistID, varID2, C_CHAR_"varname2"//C_NULL_CHAR)

!     Create a Time axis
      taxisID = taxisCreate(TAXIS_ABSOLUTE)

!     Assign the Time axis to the variable list
      CALL vlistDefTaxis(vlistID, taxisID)

!     Create a dataset in netCDF format
      streamID = streamOpenWrite(C_CHAR_"example.nc"//C_NULL_CHAR, FILETYPE_NC)
      IF ( streamID < 0 ) THEN
         msg => cdiStringError(streamID)
         WRITE(0,'(132a)') msg
         STOP 1
      END IF

!     Assign the variable list to the dataset
      CALL streamDefVlist(streamID, vlistID)

!     Loop over the number of time steps
      DO tsID = 0, nts-1
!        Set the verification date to 1985-01-01 + tsID
         CALL taxisDefVdate(taxisID, 19850101+tsID)
!        Set the verification time to 12:00:00
         CALL taxisDefVtime(taxisID, 120000)
!        Define the time step
         status = streamDefTimestep(streamID, tsID)

!        Init var1 and var2
         DO i = 1, nlon*nlat
            var1(i) = 1.1
         END DO
         DO i = 1, nlon*nlat*nlev
            var2(i) = 2.2
         END DO

!        Write var1 and var2
         CALL streamWriteVar(streamID, varID1, var1, nmiss)
         CALL streamWriteVar(streamID, varID2, var2, nmiss)
      END DO

!     Close the output stream
      CALL streamClose(streamID)

!     Destroy the objects
      CALL vlistDestroy(vlistID)
      CALL taxisDestroy(taxisID)
      CALL zaxisDestroy(zaxisID1)
      CALL zaxisDestroy(zaxisID2)
      CALL gridDestroy(gridID)

      END PROGRAM CDIWRITEF2003
