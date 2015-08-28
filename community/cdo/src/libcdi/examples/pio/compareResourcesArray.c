#if defined (HAVE_CONFIG_H)
#include "config.h"
#endif

#include <stdio.h>

#include <mpi.h>
#include <yaxt.h>
#include "cdi.h"
#include "cdipio.h"
#include "dmemory.h"
#include "pio_util.h"
#include "resource_handle.h"
#include "resource_unpack.h"

extern int    reshListCompare      ( int, int );

enum {
  IOMode           = PIO_NONE,
  nProcsIO         = 1,
  DOUBLE_PRECISION = 8,
  nlon             = 12,
  nlat             = 6,
  nlev             = 5,
  ntsteps          = 3 };

double lons[nlon] = {0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330};
double lats[nlat] = {-75, -45, -15, 15, 45, 75};
double levs[nlev] = {101300, 92500, 85000, 50000, 20000};

static int defineGrid()
{
  int gridID = CDI_UNDEFID;
  int mask_vec[nlon*nlat];
  const int * mp = &mask_vec[0];
  double area_vec[nlon*nlat];
  const double * ap = &area_vec[0];
  int i;
  
  gridID = gridCreate(GRID_LONLAT, nlon*nlat);
  gridDefXsize(gridID, nlon);
  gridDefYsize(gridID, nlat);
  gridDefXvals(gridID, lons);
  gridDefYvals(gridID, lats);
  gridDefNvertex ( gridID, 1 );
  gridDefXbounds ( gridID, lons );
  gridDefYbounds ( gridID, lats );
  for ( i = 0; i < nlon*nlat; i++ )    
    mask_vec[i] = i % 2 ;
  gridDefMaskGME ( gridID, mp );
  for ( i = 0; i < nlon*nlat; i++ )    
    mask_vec[i] = 1;
  gridDefMask ( gridID, mp );
  gridDefXname ( gridID, "myXname" );
  gridDefXlongname ( gridID, "myXlongname" );
  gridDefXunits ( gridID, "myXunits" );
  gridDefYname ( gridID, "myYname" );
  gridDefYlongname ( gridID, "myYlongname" );
  gridDefYunits ( gridID, "myYunits" );
  gridDefPrec ( gridID, DOUBLE_PRECISION );
  gridDefXpole ( gridID, 90.0 );
  gridDefYpole ( gridID, 180.0 );
  gridDefAngle ( gridID, 360.0 );
  gridDefTrunc ( gridID, 1 );
  gridDefGMEnd ( gridID, 2 );
  gridDefGMEni ( gridID, 3 );
  gridDefGMEni2 ( gridID, 4 );
  gridDefGMEni3 ( gridID, 5 );
  gridDefNumber ( gridID, 6 );
  gridDefPosition ( gridID, 7 );
  gridDefReference ( gridID, "myReference" );
/* gridDefLCC ( gridID, double originLon, double originLat,  */
/* 	   double lonParY, double lat1, double lat2, double xinc, double yinc, int projflag, int scanflag); */
/* gridDefLcc2 ( gridID, double earth_radius, double lon_0,  */
/* 	    double lat_0, double lat_1,double lat_2);*/
/* gridDefLaea ( gridID, double earth_radius, double lon_0, double lat_0); */
  for ( i = 0; i < nlon*nlat; i++ )    
    area_vec[i] = 0.1 * i;
  gridDefArea ( gridID, ap );
  for ( i = 0; i < nlon*nlat; i++ )    
    mask_vec[i] = i;
  gridDefRowlon ( gridID, nlon*nlat, mp );
  gridDefComplexPacking ( gridID, 1 );	       

  return gridID;
}

static int defineZaxis ()
{
  int zaxisID = CDI_UNDEFID;
  double vct[3] = { 3.0, 3.3, 3.6 };

  zaxisID = zaxisCreate(ZAXIS_PRESSURE, nlev);
  zaxisDefLevels(zaxisID, levs);
  zaxisDefLevel ( zaxisID, 2, 8507.3 );
  zaxisDefName ( zaxisID, "myName" );
  zaxisDefLongname ( zaxisID, "myLongname" );
  zaxisDefUnits ( zaxisID, "myUnits" );
  zaxisDefPrec ( zaxisID, DOUBLE_PRECISION );
  zaxisDefLtype ( zaxisID, 1 );
  zaxisDefVct ( zaxisID, 3, vct );
  zaxisDefLbounds ( zaxisID, &levs[0] );
  zaxisDefUbounds ( zaxisID, &levs[0] );
  zaxisDefWeights ( zaxisID, &levs[0] );

  return zaxisID;
}

static int defineTaxis ()
{
  int taxisID = CDI_UNDEFID;

  taxisID = taxisCreate(TAXIS_ABSOLUTE);
  
  taxisDefType  ( taxisID, 0 );
  taxisDefVdate ( taxisID, 1 );
  taxisDefVtime ( taxisID, 2 );
  taxisDefRdate ( taxisID, 3 );
  taxisDefRtime ( taxisID, 4 );
  taxisDefVdateBounds ( taxisID, 5, 6 );
  taxisDefVtimeBounds ( taxisID, 7, 8 );
  taxisDefCalendar ( taxisID, 1 );
  taxisDefTunit ( taxisID, 1 );
  taxisDefNumavg ( taxisID, 1 );

  return taxisID;
}

void defineStream ( int streamID, int vlistID )
{
  streamDefByteorder ( streamID, 1 );
  streamDefCompType  ( streamID, 2 );
  streamDefCompLevel ( streamID, 3 );
  streamDefVlist(streamID, vlistID);
}

int defineVlist ( int gridID, int zaxisID, int taxisID )
{
  int vlistID = CDI_UNDEFID;
  int zaxisID2 = zaxisCreate(ZAXIS_SURFACE, 1);
  int varID1, varID2;

  vlistID = vlistCreate();
  varID1 = vlistDefVar(vlistID, gridID, zaxisID, TIME_VARIABLE);
  varID2 = vlistDefVar(vlistID, gridID, zaxisID2, TIME_VARIABLE);
  vlistDefVarName(vlistID, varID1, "varname1");
  {
    int globfac[] = { 23, 42 };
    vlistDefAttInt(vlistID, varID1, "seer's globule factors", DATATYPE_INT16,
                   2, globfac);
  }
  vlistDefVarName(vlistID, varID2, "varname2");
  vlistDefAttTxt(vlistID, varID2, "txt demo", 6, "banana");
  vlistDefTaxis(vlistID, taxisID);
  return vlistID;
}

static int
defineInstitute ()
{
  int instID = CDI_UNDEFID;

  instID = institutDef( 0, 0,"MYINSTITUTE", "myInstitute");

  return instID;
}

static void
defineModel(int instID)
{
  modelDef(instID, 0, "myModel");
}

static void
printResources()
{
  FILE *fp = fopen("reshArrayModel.txt", "w");
  if (!fp)
    xabort("%s", "could not open file");
  reshListPrint(fp);
  fclose(fp);
}


static void modelRun ( MPI_Comm comm )
{
  int gridID, zaxisID, taxisID, instID, vlistID, streamID;

  char * recvBuffer, * sendBuffer;
  int bufferSize, differ;
  MPI_Status status;

  namespaceSetActive ( 0 );

  gridID  = defineGrid      ();
  zaxisID = defineZaxis     ();
  taxisID = defineTaxis     ();
  instID  = defineInstitute ();
  defineModel(instID);
  vlistID = defineVlist     ( gridID, zaxisID, taxisID);
  streamID = streamOpenWrite("example.grb", FILETYPE_GRB);
  if ( streamID < 0 ) xabort ( "Could not open file" );
  defineStream ( streamID, vlistID );

  reshPackBufferCreate(&sendBuffer, &bufferSize, &comm);
  xmpi ( MPI_Send ( sendBuffer, bufferSize, MPI_PACKED, 0, 0, comm ));
  recvBuffer = xmalloc((size_t)bufferSize);
  xmpi ( MPI_Recv ( recvBuffer, bufferSize, MPI_PACKED, 0,
		    0, comm, &status ));

  namespaceSetActive ( 1 );
  reshUnpackResources(recvBuffer, bufferSize, &comm);
  free ( recvBuffer );
  reshPackBufferDestroy ( &sendBuffer );

  differ = reshListCompare ( 0, 1 );
  printf ( "The resource arrays %s.\n", differ ? "differ" : "are equal" );
  printResources();

  namespaceSetActive ( 0 );
  streamClose(streamID);
  return;
}

int main (int argc, char *argv[])
{
  int sizeGlob, pioNamespace;
  MPI_Comm commGlob, commModel;

  MPI_Init(&argc, &argv);
  commGlob = MPI_COMM_WORLD;
  xt_initialize(commGlob);
  xmpi ( MPI_Comm_set_errhandler ( commGlob, MPI_ERRORS_RETURN ));
  xmpi ( MPI_Comm_size ( commGlob, &sizeGlob ));

  if ( sizeGlob != 1 )
      xabort ( "test transition of resource array only with 1 PE." );

  if ( nProcsIO != 1 )
    xabort ( "bad distribution of tasks on PEs" );

  commModel = pioInit(commGlob, nProcsIO, IOMode, &pioNamespace, 1.0f,
                      cdiPioNoPostCommSetup);
  namespaceSetActive(pioNamespace);

  modelRun ( commModel );

  xt_finalize();
  MPI_Finalize ();

  return 0;
}

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
