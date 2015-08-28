#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#ifndef _XOPEN_SOURCE
#define _XOPEN_SOURCE 600
#endif

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <errno.h>
#include <unistd.h>
#include <string.h>

#include "pio_impl.h"
#include "pio_util.h"

bool localDebug    =false;


int dbuffer_init ( struct dBuffer **dbuffer, size_t size )
{
  struct dBuffer *db;
  int status;
  size_t pagesize;


#ifndef _SX
  pagesize = ( size_t ) sysconf ( _SC_PAGESIZE );

  if ( localDebug ) 
    fprintf ( stdout, "dbuffer_init(): pagesize = %zu bytes, size = %zu \n", pagesize, size );

  if ( dbuffer == NULL || size < pagesize )
    {

      fprintf ( stdout, "dbuffer_init: dbuffer=NULL\n" );
      return 1;
    }
  db = ( struct dBuffer * ) malloc ( sizeof ( struct dBuffer ));

  if ( db == NULL )
    {
      perror ( "Not enough memory" );
      return 1;
    }
  
  db->size = pagesize;
  while ( db->size < size )
    {
      db->size <<= 1;
      if ( localDebug ) 
	fprintf ( stdout,"size correction: %zu\n", db->size );
    }
  
  db->wr_pointer = 0;

  if ( ( status = posix_memalign ( ( void ** ) &db->buffer, pagesize, 
				   sizeof ( char ) * ( db->size ))) != 0 ) 
    {
      switch ( status )
	{
	case EINVAL:
	  fprintf ( stderr, 
		    "The alignment argument was not a power of two, or was not a multiple of sizeof(void *).\n" );
	  break;
	case ENOMEM:
	  fprintf ( stderr, 
		    "There was insufficient memory to fulfill the allocation request.\n" );
	  break;
	}
    }
#else

  if ( dbuffer == NULL )
    {
      fprintf ( stdout, "dbuffer_init: dbuffer=NULL\n" );
      return 1;
    }

  db = ( struct dBuffer * ) malloc ( sizeof ( struct dBuffer ));
  
  if ( db == NULL )
    {
      perror ( "Not enough memory" );
      return 1;
    }
  
  db->size = size;
  
  db->wr_pointer = 0;

  db->buffer = ( unsigned char * ) malloc ( sizeof ( unsigned char ) * ( db->size ));
  if ( db->buffer == NULL )
  {
      perror ( "Not enough memory" );
      return 1 ;
  }
#endif

  *dbuffer = db;
  
  return 0;
}

void dbuffer_cleanup ( struct dBuffer **dbuffer )
{
  struct dBuffer *db;

  db = *dbuffer;

  free ( db->buffer );
  free ( db );

  return;
}

size_t dbuffer_data_size ( struct dBuffer *dbuffer )
{
  size_t data_size;

  data_size = ( size_t )( dbuffer->wr_pointer & ( dbuffer->size-1 ));

  return data_size;
}

static size_t
dbuffer_freesize(struct dBuffer *dbuffer)
{
  size_t free_size;

  free_size = ( size_t )( dbuffer->size - 1 - dbuffer_data_size ( dbuffer ));

  return free_size;
}

int dbuffer_reset ( struct dBuffer *dbuffer )
{
  dbuffer->wr_pointer = 0;
  
  return 0;
}

int
dbuffer_push(struct dBuffer *dbuffer, const void *buffer, size_t len)
{
  size_t space_left;
  size_t wr_ptr;

  space_left = dbuffer_freesize(dbuffer);
  if ( len > space_left )
    {
      return 1; /* not enough space left */
    }
  
  wr_ptr = dbuffer->wr_pointer;
  memcpy ( dbuffer->buffer + wr_ptr, buffer, len );
  dbuffer->wr_pointer = wr_ptr + len;

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
