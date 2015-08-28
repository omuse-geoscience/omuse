#include <stdio.h>
#include "cdi.h"

int main(void)
{
  const int nlon = 12; // Number of longitudes
  const int nlat =  6; // Number of latitudes 
  const int nlev =  5; // Number of levels    
  const int nts  =  3; // Number of time steps
  int nmiss;
  double var1[nlon*nlat];
  double var2[nlon*nlat*nlev];


  // Open the input dataset
  int streamID1  = streamOpenRead("example.nc");
  if ( streamID1 < 0 )
    {
      fprintf(stderr, "%s\n", cdiStringError(streamID1));
      return 1;
    }

  // Get the variable list of the dataset
  int vlistID1 = streamInqVlist(streamID1);

  // Set the variable IDs
  int varID1 = 0;
  int varID2 = 1;

  // Open the output dataset (GRIB format)
  int streamID2  = streamOpenWrite("example.grb", FILETYPE_GRB);
  if ( streamID2 < 0 )
    {
      fprintf(stderr, "%s\n", cdiStringError(streamID2));
      return 1;
    }

  int vlistID2 = vlistDuplicate(vlistID1);

  streamDefVlist(streamID2, vlistID2);

  // Loop over the number of time steps
  for ( int tsID = 0; tsID < nts; tsID++ )
    {
      // Inquire the input time step
      streamInqTimestep(streamID1, tsID);

      // Define the output time step
      streamDefTimestep(streamID2, tsID);

      // Read var1 and var2 
      streamReadVar(streamID1, varID1, var1, &nmiss);
      streamReadVar(streamID1, varID2, var2, &nmiss);

      // Write var1 and var2
      streamWriteVar(streamID2, varID1, var1, nmiss);
      streamWriteVar(streamID2, varID2, var2, nmiss);
    }

  // Close the streams
  streamClose(streamID1);
  streamClose(streamID2);

  return 0;
}
