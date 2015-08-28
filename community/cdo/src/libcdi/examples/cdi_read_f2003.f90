PROGRAM CDIREADF2003
  use iso_c_binding
  use mo_cdi

  IMPLICIT NONE

  INTEGER :: gsize, nlevel, nvars, code
  INTEGER :: vdate, vtime, nmiss, status, ilev
  INTEGER :: streamID, varID, gridID, zaxisID
  INTEGER :: tsID, vlistID, taxisID
  DOUBLE PRECISION, ALLOCATABLE :: field(:,:)
  CHARACTER(kind=c_char), POINTER, DIMENSION(:) :: &
       msg, cdi_version
  CHARACTER(kind=c_char, LEN = cdi_max_name + 1) :: &
       name, longname, units
  INTEGER :: name_c_len, longname_c_len, units_c_len

  cdi_version => cdiLibraryVersion()

  WRITE (0, '(a,132a)') 'cdi version: ', cdi_version

  ! Open the dataset
  streamID = streamOpenRead(C_CHAR_"example.nc"//C_NULL_CHAR)
  IF ( streamID < 0 ) THEN
    PRINT *,  'Could not Read the file.'
    msg => cdiStringError(streamID)
    WRITE(0,'(132a)') msg
    STOP 1
  END IF

  ! Get the variable list of the dataset
  vlistID = streamInqVlist(streamID)

  nvars = vlistNvars(vlistID)

  DO varID = 0, nvars-1
    code = vlistInqVarCode(vlistID, varID)
    CALL vlistInqVarName(vlistID, varID, name)
    CALL vlistInqVarLongname(vlistID, varID, longname)
    CALL vlistInqVarUnits(vlistID, varID, units)

    ! CALL ctrim(name)
    ! CALL ctrim(longname)
    ! CALL ctrim(units)

    longname_c_len = c_len(longname)
    name_c_len = c_len(name)
    units_c_len = c_len(units)
    PRINT '(a,2(i0,a),132a)', 'Parameter: ', varID+1, ' ', code, ' ', &
         name(1:name_c_len), ' ', longname(1:longname_c_len), ' ', &
         units(1:units_c_len), ' |'

  END DO

  ! Get the Time axis form the variable list
  taxisID = vlistInqTaxis(vlistID)

  ! Loop over the time steps
  DO tsID = 0, 999999
    ! Read the time step
    status = streamInqTimestep(streamID, tsID)
    IF ( status == 0 ) exit

    ! Get the verification date and time
    vdate = taxisInqVdate(taxisID)
    vtime = taxisInqVtime(taxisID)

    PRINT '(a,i3,i10,i10)', 'Timestep: ', tsID+1, vdate, vtime

    ! Read the variables at the current timestep
    DO varID = 0, nvars-1
      gridID = vlistInqVarGrid(vlistID, varID)
      gsize = gridInqSize(gridID)
      zaxisID = vlistInqVarZaxis(vlistID, varID)
      nlevel = zaxisInqSize(zaxisID)
      ALLOCATE(field(gsize, nlevel))
      CALL streamReadVar(streamID, varID, field, nmiss)
      DO ilev = 1, nlevel
        PRINT '(a,i3,a,i3,a,f10.5,1x,f10.5)', '   var=', varID+1, &
             ' level=', ilev, ':', &
             MINVAL(field(:,ilev)), MAXVAL(field(:,ilev))
      END DO
      DEALLOCATE(field)
    END DO
  END DO

  ! Close the input stream
  CALL streamClose(streamID)

END PROGRAM CDIREADF2003
