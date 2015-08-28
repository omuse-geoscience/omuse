#include <cdi.h>
#include "cdo.h"

void nospec(int vlistID)
{
  int gridID, gridtype;
  int varID, nvars;

  nvars = vlistNvars(vlistID);
  for ( varID = 0; varID < nvars; varID++ )
    {
      gridID = vlistInqVarGrid(vlistID, varID);
      gridtype = gridInqType(gridID);
      if ( gridtype == GRID_SPECTRAL )
	cdoAbort("Operator not defined for spectral fields");
    }
}
