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
 
 Eofcoeff             eofcoeff             process eof coefficients
*/
#define WEIGHTS 1

#include "cdi.h"
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"
#include "grid.h"


// NO MISSING VALUE SUPPORT ADDED SO FAR

void *Eofcoeff3d(void * argument)
{
  char eof_name[6], oname[1024], filesuffix[32];
  const char *refname;
  //double *w;
  double missval1 = -999, missval2 = -999;
  double *xvals, *yvals, *zvals;  
  field_t ***eof;  
  field_t in;  
  field_t **out;
  //int operatorID, operfunc;  
  int gridsize;
  int i, varID, recID, levelID, tsID, eofID, *varID3;    
  int gridID1,gridID3;
  int nrecs, nvars, nlevs, neof, nchars, nmiss; 
  int streamID1, streamID2, *streamIDs;
  int taxisID2, taxisID3;
  int vlistID1, vlistID2, vlistID3;
  int zaxisID3;
   
  cdoInitialize(argument);
  cdoOperatorAdd("eofcoeff3d",  0,       0, NULL);
  //  operatorID = cdoOperatorID();
  //  operfunc = cdoOperatorFunc(operatorID);
     
  streamID1 = streamOpenRead(cdoStreamName(0));
  streamID2 = streamOpenRead(cdoStreamName(1));
  
  vlistID1 = streamInqVlist(streamID1);
  vlistID2 = streamInqVlist(streamID2);
  
  //taxisID1 = vlistInqTaxis(vlistID1);  
  taxisID2 = vlistInqTaxis(vlistID2); 
  taxisID3 = taxisDuplicate(taxisID2);
  
  gridID1 = vlistInqVarGrid(vlistID1, 0);
  //gridID2 = vlistInqVarGrid(vlistID2, 0);
  
  if ( vlistGridsizeMax(vlistID1)==vlistGridsizeMax(vlistID2) )
    gridsize = vlistGridsizeMax(vlistID1);  
  else 
    {
      gridsize = -1;
      cdoAbort("Gridsize of input files does not match");
    }
      
  
  if ( vlistNgrids(vlistID2) > 1 || vlistNgrids(vlistID1) > 1 )
    cdoAbort("Too many grids in input");
  
  nvars = vlistNvars(vlistID1)==vlistNvars(vlistID2) ? vlistNvars(vlistID1) : -1;
  nrecs = vlistNrecs(vlistID1); 
  nlevs = zaxisInqSize(vlistInqVarZaxis(vlistID1, 0));
  //w = (double*) malloc(gridsize*sizeof(double));
  //gridWeights(gridID2, w);
  
  
  if (vlistGridsizeMax(vlistID2)   != gridsize ||
      vlistInqVarGrid(vlistID2, 0) != gridID1 )
    cdoAbort("EOFs (%s) and data (%s) defined on different grids", cdoStreamName(0)->args, cdoStreamName(1)->args);    
 
  strcpy(oname, cdoStreamName(2)->args);
  nchars = strlen(oname);
  
  refname = cdoStreamName(0)->argv[cdoStreamName(0)->argc-1];
  filesuffix[0] = 0;
  cdoGenFileSuffix(filesuffix, sizeof(filesuffix), streamInqFiletype(streamID1), vlistID1, refname);
 
  eof = (field_t***) malloc(nvars * sizeof(field_t**));
  for ( varID=0; varID<nvars; varID++)
    eof[varID] = (field_t**) malloc(nlevs*sizeof(field_t*));

  eofID = 0;
  while ( 1 )       
   {     
     nrecs = streamInqTimestep(streamID1, eofID);
     if ( nrecs == 0) break;

     for ( recID = 0; recID < nrecs; recID++ )
       {         
         streamInqRecord(streamID1, &varID, &levelID);
         missval1 = vlistInqVarMissval(vlistID1, varID);
         if ( eofID == 0 )
           eof[varID][levelID] = (field_t*) malloc(1*sizeof(field_t));
         else
           eof[varID][levelID] = (field_t*) realloc(eof[varID][levelID], (eofID+1)*sizeof(field_t));
         eof[varID][levelID][eofID].grid   = gridID1;
         eof[varID][levelID][eofID].nmiss  = 0;
         eof[varID][levelID][eofID].missval= missval1;
         eof[varID][levelID][eofID].ptr    = (double*) malloc(gridsize*sizeof(double));
         memset(&eof[varID][levelID][eofID].ptr[0], missval1, gridsize*sizeof(double));

         if ( varID >= nvars )
           cdoAbort("Internal error - too high varID");
         if ( levelID >= nlevs )
           cdoAbort("Internal error - too high levelID");
         
         streamReadRecord(streamID1, eof[varID][levelID][eofID].ptr, 
                          &eof[varID][levelID][eofID].nmiss);
       }
     eofID++;
   }
  neof = eofID;  
  
  if ( cdoVerbose ) cdoPrint("%s contains %i eof's", cdoStreamName(0)->args, neof);
  // Create 1x1 Grid for output
  gridID3 = gridCreate(GRID_LONLAT, 1);
  gridDefXsize(gridID3, 1);
  gridDefYsize(gridID3, 1);
  xvals= (double*) malloc(1*sizeof(double));
  yvals= (double*) malloc(1*sizeof(double));
  xvals[0]=0;
  yvals[0]=0;
  gridDefXvals(gridID3, xvals);
  gridDefYvals(gridID3, yvals);
  
  zvals = (double *) malloc( 1* sizeof(double ));
  zvals[0] = 0.;
  zaxisID3 = zaxisCreate(ZAXIS_GENERIC,1);
  zaxisDefLevels(zaxisID3,zvals);
  zaxisDefName(zaxisID3,"zaxis_Reduced");
  zaxisDefLongname(zaxisID3,"Reduced zaxis from EOF3D - only one coefficient per 3D eigenvector and time step");
  
  vlistID3 = vlistCreate();
  vlistDefTaxis(vlistID3,taxisID3);
  varID3 = (int*) malloc( nvars * sizeof(int));
  for ( varID=0; varID<nvars; varID++ )
    varID3[varID] = vlistDefVar(vlistID3, gridID3, zaxisID3, TSTEP_INSTANT);
  
  // open streams for eofcoeff output
  streamIDs = (int*) malloc(neof*sizeof(int)); 
  eofID = 0;
  for ( eofID = 0; eofID < neof; eofID++)
    {
      oname[nchars] = '\0';                       

      sprintf(eof_name, "%5.5i", eofID);
      strcat(oname, eof_name);
      if ( filesuffix[0] )
        strcat(oname, filesuffix);
      
      argument_t *fileargument = file_argument_new(oname);
      streamIDs[eofID] = streamOpenWrite(fileargument, cdoFiletype());
      file_argument_free(fileargument);

      if (cdoVerbose) 
        cdoPrint("opened %s ('w')  as stream%i for %i. eof", oname, streamIDs[eofID], eofID+1);
      
      streamDefVlist(streamIDs[eofID], vlistID3);
    }
  
  // ALLOCATE temporary fields for data read and write
  in.ptr = (double*) malloc(gridsize*sizeof(double));
  in.grid = gridID1;  
  out = (field_t**) malloc(nvars*sizeof(field_t*));
  for ( varID = 0; varID < nvars; varID++ ) {
    out[varID] = (field_t*) malloc( neof * sizeof(field_t));
    for ( eofID=0; eofID<neof; eofID++ ) {
      out[varID][eofID].missval = missval1;
      out[varID][eofID].nmiss = 0;
      out[varID][eofID].ptr = (double*) malloc(1*sizeof(double));
    }
  }

  tsID=0;
  while ( 1 )
    {      
      nrecs = streamInqTimestep(streamID2, tsID);
      if ( nrecs == 0 ) break;

      for ( varID = 0; varID < nvars; varID++ )
	for ( eofID = 0; eofID < neof; eofID++ ) {
	  out[varID][eofID].ptr[0]  = 0;
	  out[varID][eofID].grid    = gridID3;
	  out[varID][eofID].missval = missval2;            
	}

      taxisCopyTimestep(taxisID3, taxisID2);

      for ( recID =0; recID< nrecs; recID++ )
        {
          streamInqRecord(streamID2, &varID, &levelID);
          missval2 = vlistInqVarMissval(vlistID2, varID);
          streamReadRecord(streamID2, in.ptr, &in.nmiss);  

          
          for (eofID = 0; eofID < neof; eofID++ )
            {
              if ( recID == 0 ) streamDefTimestep(streamIDs[eofID],tsID);

	      nmiss = 0;
              for(i=0;i<gridsize;i++)
                {                  
                  if (! DBL_IS_EQUAL(in.ptr[i],missval2) && 
                      ! DBL_IS_EQUAL(eof[varID][levelID][eofID].ptr[i],missval1 ) )
                    {
                      // double tmp = w[i]*in.ptr[i]*eof[varID][levelID][eofID].ptr[i];
                      double tmp = in.ptr[i]*eof[varID][levelID][eofID].ptr[i];
                      out[varID][eofID].ptr[0] += tmp;                   
                    }
		  else
		    nmiss += 1;
                }            

              if ( nmiss ) {
		out[varID][eofID].nmiss=1;
		out[varID][eofID].ptr[0]=missval2;
	      }
            }

          if ( varID >= nvars )
            cdoAbort("Internal error - too high varID");
          if ( levelID >= nlevs )
            cdoAbort("Internal error - too high levelID");                              
        }

      for ( eofID = 0; eofID < neof; eofID++ ) {
	for ( varID = 0; varID < nvars; varID++ ) {
	  streamDefRecord(streamIDs[eofID], varID, 0);
	  streamWriteRecord(streamIDs[eofID],out[varID][eofID].ptr,out[varID][eofID].nmiss);
	}
      }

      tsID++;
    }
  
  for ( eofID = 0; eofID < neof; eofID++) streamClose(streamIDs[eofID]);

  streamClose(streamID2);
  streamClose(streamID1);
  
  cdoFinish();

  return (0);
}

