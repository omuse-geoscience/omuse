#include <stdio.h>
#include "cdi.h"

int main(void)
{
  const int nlon = 12; // Number of longitudes
  const int nlat =  6; // Number of latitudes
  const int nlev =  5; // Number of levels
  const int nts  =  3; // Number of time steps
  int nmiss, vdate, vtime;
  double var1[nlon*nlat];
  double var2[nlon*nlat*nlev];


  // Open the dataset 
  int streamID = streamOpenRead("example.nc");
  if ( streamID < 0 )
    {
      fprintf(stderr, "%s\n", cdiStringError(streamID));
      return 1;
    }

  // Get the variable list of the dataset 
  int vlistID = streamInqVlist(streamID);

  // Set the variable IDs 
  int varID1 = 0;
  int varID2 = 1;

  // Get the Time axis from the variable list 
  int taxisID = vlistInqTaxis(vlistID);

  // Loop over the number of time steps 
  for ( int tsID = 0; tsID < nts; tsID++ )
    {
      // Inquire the time step 
      streamInqTimestep(streamID, tsID);

      // Get the verification date and time 
      vdate = taxisInqVdate(taxisID);
      vtime = taxisInqVtime(taxisID);
      printf("read timestep %d:  date=%d  time=%d\n", tsID+1, vdate, vtime);

      // Read var1 and var2 
      streamReadVar(streamID, varID1, var1, &nmiss);
      streamReadVar(streamID, varID2, var2, &nmiss);
    }

  // Close the input stream 
  streamClose(streamID);

  return 0;
}
