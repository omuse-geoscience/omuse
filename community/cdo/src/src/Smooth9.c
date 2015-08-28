/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2015 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

/*
   This module contains the following operators:

      Smoothstat       smooth9             running 9-point-average
*/


#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"

void *Smooth9(void *argument)
{
  int operatorID;
  int streamID1, streamID2;
  int gridsize;
  int gridID;
  int nrecs, recID;
  int tsID;
  int varID, levelID;
  int vlistID1, vlistID2;
  int nmiss, nmiss2;
  int i, j , i2;
  double avg,divavg;
  double missval1, missval2;
  double *array1, *array2;
  int taxisID1, taxisID2;
  int nlon, nlat;
  int gridtype;
  int nvars;
  int grid_is_cyclic;
  int *varIDs = NULL, *mask =  NULL; 

  cdoInitialize(argument);

  cdoOperatorAdd("smooth9",   0,   0, NULL);
 
  operatorID = cdoOperatorID();
  UNUSED(operatorID);

  streamID1 = streamOpenRead(cdoStreamName(0));

  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = vlistDuplicate(vlistID1);

  taxisID1 = vlistInqTaxis(vlistID1);
  taxisID2 = taxisDuplicate(taxisID1);
  vlistDefTaxis(vlistID2, taxisID2);

  nvars = vlistNvars(vlistID1);
  varIDs  = (int*) malloc(nvars*sizeof(int)); 

  for ( varID = 0; varID < nvars; ++varID )
    {
      gridID = vlistInqVarGrid(vlistID1, varID);
      gridtype = gridInqType(gridID);
      if ( gridtype == GRID_GAUSSIAN ||
           gridtype == GRID_LONLAT   ||
           gridtype == GRID_CURVILINEAR )
	{
	  varIDs[varID] = 1;
	}
      else
	{
	  varIDs[varID] = 0;
	  cdoWarning("Unsupported grid for varID %d", varID);
	}
    }

  gridsize = vlistGridsizeMax(vlistID1);
  array1 = (double*) malloc(gridsize*sizeof(double));
  array2 = (double*) malloc(gridsize *sizeof(double));
  mask   = (int*) malloc(gridsize *sizeof(int));
 
  streamID2 = streamOpenWrite(cdoStreamName(1), cdoFiletype());

  streamDefVlist(streamID2, vlistID2);

  tsID = 0;
  while ( (nrecs = streamInqTimestep(streamID1, tsID)) )
    {
      taxisCopyTimestep(taxisID2, taxisID1);

      streamDefTimestep(streamID2, tsID);

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID1, &varID, &levelID);
	  streamReadRecord(streamID1, array1, &nmiss);
	
	  if ( varIDs[varID] )
	    {	    
	      missval1 = vlistInqVarMissval(vlistID1, varID);
	      missval2 = missval1;

	      gridID = vlistInqVarGrid(vlistID1, varID);
	      gridsize = gridInqSize(gridID);
	      nlon = gridInqXsize(gridID);	 
	      nlat = gridInqYsize(gridID);
	      grid_is_cyclic = gridIsCircular(gridID);
	   
	      for ( i = 0; i < gridsize; i++) 
		{		
		  if ( DBL_IS_EQUAL(missval1, array1[i]) ) mask[i] = 0;
		  else mask[i] = 1;
		}
 
	      nmiss2=0;
	      for ( i = 0; i < nlat; i++ )
		{
		  for ( i2 = 0; i2 < nlon; i2++ )
		    {		      
		      avg = 0; divavg = 0; 	  		     			

		      if ( ( i == 0 ) || ( i2 == 0 )  ||
                           ( i == ( nlat - 1 ) ) || ( i2 == ( nlon -1 ) ) )
			{
			  j = i2+nlon*i;
			  if ( mask[j] )
			    {
			      avg += array1[j];  divavg+= 1;					     		       
			      /* upper left corner */
			      if ( (  i!=0 ) && ( i2!=0 ) ) 
				{ 
				  j = ((i-1)*nlon)+i2-1;
				  if ( mask[j] ) 
				    { avg +=   0.3*array1[j]; divavg+=0.3;}
				}
			      else if ( i != 0 && grid_is_cyclic ) 
				{ 
				  j = (i-1)*nlon+i2-1+nlon;
				  if ( mask[j] ) 
				    { avg+=0.3*array1[j]; divavg+=0.3;}
				}
			      
			      /* upper cell */
			      if ( i!=0 ) 
				{
				  j = ((i-1)*nlon)+i2;
				  if ( mask[j] )
				    { avg +=  0.5*array1[j];  divavg+= 0.5; }
				}
			      
			      /* upper right corner */
			      if ( ( i!=0) && ( i2!=(nlon-1) ) ) 
				{
				  j = ((i-1)*nlon)+i2+1;
				  if ( mask[j] )
				    { avg +=   0.3*array1[j]; divavg+= 0.3; }
				}
			      else if ( i!= 0 && grid_is_cyclic )
				{ 
				  j = (i-1)*nlon+i2+1-nlon;
				  if ( mask[j] ) 
				    {avg+=0.3*array1[j]; divavg+=0.3;}
				}
			      
			      /* left cell */
			      if  ( i2!=0 ) 
				{
				  j = ((i)*nlon)+i2-1;
				  if ( mask[j] )
				    {  avg +=   0.5*array1[j];  divavg+= 0.5;}
				}
			      else if ( grid_is_cyclic )
				{
				  j = i*nlon-1+nlon;
				  if ( mask[j] ) 
				    { avg+=0.5*array1[j]; divavg+=0.5;}
				}
			      
			      /* right cell */
			      if ( i2!=(nlon-1) ) 
				{ 
				  j = (i*nlon)+i2+1;
				  if ( mask[j] )
				    {  avg +=   0.5*array1[j]; divavg+= 0.5; } 
				}
			      else if ( grid_is_cyclic )
				{
				  j = i*nlon+i2+1-nlon;
				  if (mask[j])
				    { avg+=0.5*array1[j]; divavg+=0.5;}
				}
			      
			      /* lower left corner */
			      if ( mask[j] &&  ( (i!=(nlat-1))&& (i2!=0) ) )
				{	       
				  j = ((i+1)*nlon+i2-1);
				  if ( mask[j] )
				    { avg +=   0.3*array1[j];  divavg+= 0.3; }
				}
			      else if ( i!= nlat-1 && grid_is_cyclic ) 
				{
				  j= (i+1)*nlon-1+nlon; 
				  if ( mask[j] ) 
				    { avg+= 0.3*array1[j]; divavg+=0.3; }
				}
			      
			      /* lower cell */
			      if  ( i!=(nlat-1) ) 
				{
				  j = ((i+1)*nlon)+i2;
				  if ( mask[j] ) 
				    { avg += 0.5*array1[j];  divavg+= 0.5;  }
				}
			      
			      /* lower right corner */
			      if ( i!=(nlat-1) && (i2!=(nlon-1) ) )
				{
				  j = ((i+1)*nlon)+i2+1;
				  if ( mask[j] )
				    {  avg += 0.3*array1[j]; divavg+= 0.3; }	
				}
			      else if ( i != (nlat-1) && grid_is_cyclic )
				{
				  j= ((i+1)*nlon)+i2+1-nlon;
				  if ( mask[j] )
				    {avg+=0.3*array1[j]; divavg+=0.3;}
				}
			    }
			}
		      else if ( mask[i2+nlon*i] )
			{			 
			  avg += array1[i2+nlon*i]; divavg+= 1;
			    
			  j = ((i-1)*nlon)+i2-1;
			  if ( mask[j] )
			    { avg += 0.3*array1[j]; divavg+= 0.3; }

			  j = ((i-1)*nlon)+i2;
			  if ( mask[j] )
			    { avg += 0.5*array1[j]; divavg+= 0.5; }

			  j = ((i-1)*nlon)+i2+1;
			  if ( mask[j] )
			    { avg += 0.3*array1[j]; divavg+= 0.3; }

			  j = ((i)*nlon)+i2-1;
			  if ( mask[j] )
			    { avg += 0.5*array1[j]; divavg+= 0.5; }

			  j = (i*nlon)+i2+1;
			  if ( mask[j] )
			    { avg += 0.5*array1[j]; divavg+= 0.5; } 

			  j = ((i+1)*nlon+i2-1);	        
			  if ( mask[j] )
			    { avg += 0.3*array1[j]; divavg+= 0.3; }

			  j = ((i+1)*nlon)+i2;
			  if ( mask[j] )
			    { avg += 0.5*array1[j]; divavg+= 0.5; }

			  j = ((i+1)*nlon)+i2+1;
			  if ( mask[j] )
			    { avg += 0.3*array1[j]; divavg+= 0.3; }
			}

		      if ( fabs(divavg) > 0 )
			{
			  array2[i*nlon+i2] = avg/divavg;			
			}
		      else 
			{
			  array2[i*nlon+i2] = missval2;					
			  nmiss2++;
			}
		    }			    	     
		}    
	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, array2, nmiss2);		
	    }     	   
	  else 
	    {
	      streamDefRecord(streamID2, varID, levelID);
	      streamWriteRecord(streamID2, array1, nmiss);
	    }
	}
      tsID++;
    }

  streamClose(streamID2);
  streamClose(streamID1);

  free(varIDs);
  if ( array2 ) free(array2);
  if ( array1 ) free(array1);

  cdoFinish();

  return (0);
}
