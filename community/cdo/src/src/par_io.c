#if defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#if defined(HAVE_LIBPTHREAD)
#  include <pthread.h>
#endif

#include <string.h> /* memcpy */

#include "cdo.h"
#include "cdo_int.h"
#include "par_io.h"
#include "pstream.h"


void *readRecord(void *arg)
{
  int streamID;
  int *varID, *levelID, *nmiss;
  double *array;
  read_arg_t *read_arg = (read_arg_t *) arg;

  streamID = read_arg->streamID;
  varID    = read_arg->varID;
  levelID  = read_arg->levelID;
  nmiss    = read_arg->nmiss;
  array    = read_arg->array;

  /* fprintf(stderr, "streamInqRecord: streamID = %d\n", streamID); */
  streamInqRecord(streamID, varID, levelID);
  streamReadRecord(streamID, array, nmiss);
  /* fprintf(stderr, "readRecord: varID %d levelID %d\n", *varID, *levelID); */

  return (NULL);
}


void parReadRecord(int streamID, int *varID, int *levelID, double *array, int *nmiss, par_io_t *parIO)
{
  int lpario = FALSE;
  int recID = 0, nrecs = 0;
#if defined(HAVE_LIBPTHREAD)
  pthread_t thrID;
  /* pthread_attr_t attr; */
  int rval;
#endif

#if defined(HAVE_LIBPTHREAD)
  if ( parIO )
    {
      lpario = TRUE;
      recID = parIO->recID;
      nrecs = parIO->nrecs;
      thrID = parIO->thrID;
    }
#endif

  if ( recID == 0 || lpario == FALSE )
    {
      read_arg_t read_arg;
      read_arg.streamID = streamID;
      read_arg.varID    = varID;
      read_arg.levelID  = levelID;
      read_arg.nmiss    = nmiss;
      read_arg.array    = array;

      readRecord(&read_arg);
    }
#if defined(HAVE_LIBPTHREAD)
  else
    {
      /* fprintf(stderr, "parIO1: %ld streamID %d %d %d\n", (long)thrID, streamID, recID, nrecs); */
      rval = pthread_join(thrID, NULL);
      if ( rval != 0 ) cdoAbort("pthread_join failed!");

      *varID    = parIO->varID;
      *levelID  = parIO->levelID;
      *nmiss    = parIO->nmiss;
      /* fprintf(stderr, "parIO2: %ld streamID %d %d %d\n", (long)thrID, streamID, *varID, *levelID); */
      memcpy(array, parIO->array, parIO->array_size*sizeof(double));
    }

  if ( lpario && nrecs > 1 )
    {
      read_arg_t *read_arg = &(parIO->read_arg);
      if ( (recID+1) < nrecs )
	{
	  if ( recID == 0 )
	    {
	      pthread_attr_init(&parIO->attr);
	      pthread_attr_setdetachstate(&parIO->attr, PTHREAD_CREATE_JOINABLE);
	    }

	  read_arg->streamID = streamID;
	  read_arg->varID    = &parIO->varID;
	  read_arg->levelID  = &parIO->levelID;
	  read_arg->nmiss    = &parIO->nmiss;
	  read_arg->array    = parIO->array;

	  /* fprintf(stderr, "pthread_create: streamID %d %d\n", read_arg->streamID,streamID); */
	  rval = pthread_create(&thrID, &parIO->attr, readRecord, read_arg);
	  if ( rval != 0 ) cdoAbort("pthread_create failed!");

	  /* fprintf(stderr, "thrID = %ld\n", (long) thrID); */
	  parIO->thrID = thrID;
	}
      else
	pthread_attr_destroy(&parIO->attr);
    }
#endif
}
