#include <stdio.h>
#include "cdi.h"

#define  nlon   12 // Number of longitudes
#define  nlat    6 // Number of latitudes 
#define  nhlev   6 // Number of half hybrid levels    
#define  nflev   5 // Number of full hybrid levels    
#define  nts     3 // Number of time steps

int main(void)
{
  int gridID, zaxisID1, zaxisID2, zaxisID3, zaxisID4, taxisID;
  int vlistID, varID1, varID2, varID3, varID4, streamID, tsID;
  int i, k, nmiss = 0;
  double lons[nlon] = {0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330};
  double lats[nlat] = {-75, -45, -15, 15, 45, 75};
  double levs[nhlev] = {1, 2, 3, 4, 5, 6};
  double levs2[nhlev] = {2, 3, 4, 5, 6, -1};
  double vct[2*nhlev] = {0, 6000, 10000, 16000, 8000, 0, 0, 0.0004, 0.03, 0.2, 0.7, 1.0};
  double var1[nlon*nlat];
  double var2[nlon*nlat*nflev];
  double var3[nlon*nlat*nhlev];
  double var4[nlon*nlat*nflev];


  // Create a regular lon/lat grid
  gridID = gridCreate(GRID_LONLAT, nlon*nlat);
  gridDefXsize(gridID, nlon);
  gridDefYsize(gridID, nlat);
  gridDefXvals(gridID, lons);
  gridDefYvals(gridID, lats);

  // Create a surface level Z-axis
  zaxisID1 = zaxisCreate(ZAXIS_SURFACE, 1);

  // Create a hybrid level Z-axis
  zaxisID2 = zaxisCreate(ZAXIS_HYBRID, nflev);
  zaxisDefLevels(zaxisID2, levs);
  zaxisDefVct(zaxisID2, nhlev*2, vct);

  // Create a hybrid half level Z-axis
  zaxisID3 = zaxisCreate(ZAXIS_HYBRID, nhlev);
  zaxisDefLevels(zaxisID3, levs);
  zaxisDefVct(zaxisID3, nhlev*2, vct);

  // Create a hybrid level Z-axis
  zaxisID4 = zaxisCreate(ZAXIS_HYBRID, nflev);
  zaxisDefLevels(zaxisID4, levs);
  zaxisDefLbounds(zaxisID4, levs);
  zaxisDefUbounds(zaxisID4, levs2);
  zaxisDefVct(zaxisID4, nhlev*2, vct);
 
  // Create a variable list
  vlistID = vlistCreate();

  // Define the variables
  varID1 = vlistDefVar(vlistID, gridID, zaxisID1, TSTEP_INSTANT);
  varID2 = vlistDefVar(vlistID, gridID, zaxisID2, TSTEP_INSTANT);
  varID3 = vlistDefVar(vlistID, gridID, zaxisID3, TSTEP_INSTANT);
  varID4 = vlistDefVar(vlistID, gridID, zaxisID4, TSTEP_INSTANT);

  // Define the variable names
  vlistDefVarName(vlistID, varID1, "sp");
  vlistDefVarName(vlistID, varID2, "t");
  vlistDefVarName(vlistID, varID3, "w");
  vlistDefVarName(vlistID, varID4, "u");

  // Create a Time axis
  taxisID = taxisCreate(TAXIS_ABSOLUTE);

  // Assign the Time axis to the variable list
  vlistDefTaxis(vlistID, taxisID);

  // Create a dataset in netCDF format
  streamID = streamOpenWrite("example.nc", FILETYPE_NC);
  if ( streamID < 0 )
    {
      fprintf(stderr, "%s\n", cdiStringError(streamID));
      return(1);
    }

  // Assign the variable list to the dataset
  streamDefVlist(streamID, vlistID);

  // Loop over the number of time steps
  for ( tsID = 0; tsID < nts; tsID++ )
    {
      // Set the verification date to 1985-01-01 + tsID
      taxisDefVdate(taxisID, 19850101+tsID);
      // Set the verification time to 12:00:00
      taxisDefVtime(taxisID, 120000);
      // Define the time step
      streamDefTimestep(streamID, tsID);

      // Init var1 and var2
      for ( i = 0; i < nlon*nlat; ++i ) var1[i] = 1.1;
      for ( k = 0; k < nflev; ++k )
	for ( i = 0; i < nlon*nlat; ++i ) var2[i+k*nlon*nlat] = 2.2+k;
      for ( k = 0; k < nhlev; ++k )
	for ( i = 0; i < nlon*nlat; ++i ) var3[i+k*nlon*nlat] = -2.2-k;
      for ( k = 0; k < nflev; ++k )
	for ( i = 0; i < nlon*nlat; ++i ) var4[i+k*nlon*nlat] = 100+k;
 
      // Write var1 and var2
      streamWriteVar(streamID, varID1, var1, nmiss);
      streamWriteVar(streamID, varID2, var2, nmiss);
      streamWriteVar(streamID, varID3, var3, nmiss);
      streamWriteVar(streamID, varID4, var4, nmiss);
    }

  // Close the output stream
  streamClose(streamID);

  // Destroy the objects
  vlistDestroy(vlistID);
  taxisDestroy(taxisID);
  zaxisDestroy(zaxisID1);
  zaxisDestroy(zaxisID2);
  gridDestroy(gridID);

  return 0;
}
