#include <stdio.h>
#include <string.h>

#include <cdi.h>
#include <cdo.h>
#include "dmemory.h"
#include "field.h"


void field_init(field_t *field)
{
  memset(field, 0, sizeof(field_t));
}


field_t **field_allocate(int vlistID, int ptype, int init)
{
  int nvars, nlevel;
  int varID, zaxisID, levelID;
  int gridID, gridsize;
  int nwpv; // number of words per value; real:1  complex:2
  double missval;
  field_t **field;

  nvars = vlistNvars(vlistID);

  field = (field_t **) malloc(nvars*sizeof(field_t *));

  for ( varID = 0; varID < nvars; ++varID )
    {
      nwpv     = vlistInqNWPV(vlistID, varID);
      gridID   = vlistInqVarGrid(vlistID, varID);
      gridsize = gridInqSize(gridID);
      zaxisID  = vlistInqVarZaxis(vlistID, varID);
      nlevel   = zaxisInqSize(zaxisID);
      missval  = vlistInqVarMissval(vlistID, varID);

      field[varID] = (field_t*) malloc(nlevel*sizeof(field_t));
      for ( levelID = 0; levelID < nlevel; ++levelID )
	{
	  field_init(&field[varID][levelID]);

	  field[varID][levelID].nwpv    = nwpv;
	  field[varID][levelID].grid    = gridID;
	  field[varID][levelID].nsamp   = 0;
	  field[varID][levelID].nmiss   = 0;
	  field[varID][levelID].missval = missval;
	  field[varID][levelID].ptr     = NULL;
	  field[varID][levelID].weight  = NULL;

	  if ( ptype == FIELD_ALL || ptype == FIELD_PTR )
	    {
	      field[varID][levelID].ptr = (double*) malloc(nwpv*gridsize*sizeof(double));
	      if ( init ) memset(field[varID][levelID].ptr, 0, nwpv*gridsize*sizeof(double));
	    }

	  if ( ptype == FIELD_ALL || ptype == FIELD_PTR )
	    {
	      field[varID][levelID].weight = (double*) malloc(nwpv*gridsize*sizeof(double));
	      if ( init ) memset(field[varID][levelID].weight, 0, nwpv*gridsize*sizeof(double));
	    }    
	}
    }

  return (field);
}


field_t **field_malloc(int vlistID, int ptype)
{
  return (field_allocate(vlistID, ptype, 0));
}


field_t **field_calloc(int vlistID, int ptype)
{
  return (field_allocate(vlistID, ptype, 1));
}


void field_free(field_t **field, int vlistID)
{
  int nvars, nlevel;
  int varID, levelID;

  nvars = vlistNvars(vlistID);

  for ( varID = 0; varID < nvars; ++varID )
    {
      nlevel = zaxisInqSize(vlistInqVarZaxis(vlistID, varID));
      for ( levelID = 0; levelID < nlevel; ++levelID )
	{
	  if ( field[varID][levelID].ptr )    free(field[varID][levelID].ptr);
       	  if ( field[varID][levelID].weight ) free(field[varID][levelID].weight);
	}

      free(field[varID]);
    }

  free(field);
}
