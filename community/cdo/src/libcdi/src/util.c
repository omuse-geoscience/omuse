#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#define _XOPEN_SOURCE 600

#include <errno.h>
#include <stdio.h>
#include <string.h>
#include <float.h>
#include <sys/types.h>

#include "cdi.h"
#include "cdi_int.h"
#include "dmemory.h"
#include "binary.h"


#undef  IsBigendian
#define IsBigendian()  ( u_byteorder.c[sizeof(long) - 1] )

void cdiPrintDatatypes(void)
{
  /* IsBigendian returns 1 for big endian byte order */
  static union {unsigned long l; unsigned char c[sizeof(long)];} u_byteorder = {1};

  fprintf (stderr, "+-------------+-------+\n");
  fprintf (stderr, "| types       | bytes |\n");
  fprintf (stderr, "+-------------+-------+\n");
  fprintf (stderr, "| void *      |   %3d |\n", (int) sizeof(void *));
  fprintf (stderr, "+-------------+-------+\n");
  fprintf (stderr, "| char        |   %3d |\n", (int) sizeof(char));
  fprintf (stderr, "+-------------+-------+\n");
  fprintf (stderr, "| short       |   %3d |\n", (int) sizeof(short));
  fprintf (stderr, "| int         |   %3d |\n", (int) sizeof(int));
  fprintf (stderr, "| long        |   %3d |\n", (int) sizeof(long));
  fprintf (stderr, "| long long   |   %3d |\n", (int) sizeof(long long));
  fprintf (stderr, "| size_t      |   %3d |\n", (int) sizeof(size_t));
  fprintf (stderr, "| off_t       |   %3d |\n", (int) sizeof(off_t));
  fprintf (stderr, "+-------------+-------+\n");
  fprintf (stderr, "| float       |   %3d |\n", (int) sizeof(float));
  fprintf (stderr, "| double      |   %3d |\n", (int) sizeof(double));
  fprintf (stderr, "| long double |   %3d |\n", (int) sizeof(long double));
  fprintf (stderr, "+-------------+-------+\n\n");
#define XSTRING(x)	#x
#define STRING(x)	XSTRING(x)
  fprintf (stderr, "+-------------+-----------+\n");
  fprintf (stderr, "| INT32       | %-9s |\n", STRING(INT32));
  fprintf (stderr, "| INT64       | %-9s |\n", STRING(INT64));
  fprintf (stderr, "| FLT32       | %-9s |\n", STRING(FLT32));
  fprintf (stderr, "| FLT64       | %-9s |\n", STRING(FLT64));
  fprintf (stderr, "+-------------+-----------+\n");

  if ( IsBigendian() )
    fprintf (stderr, "\n  byte ordering is BIGENDIAN\n\n");
  else
    fprintf (stderr, "\n  byte ordering is LITTLEENDIAN\n\n");
}

static char uuidFmt[] = "%02hhx%02hhx%02hhx%02hhx-"
  "%02hhx%02hhx-%02hhx%02hhx-%02hhx%02hhx-"
  "%02hhx%02hhx%02hhx%02hhx%02hhx%02hhx";

enum {
  uuidNumHexChars = 36,
};

void uuid2str(const unsigned char *uuid, char *uuidstr)
{

  if ( uuid == NULL || uuidstr == NULL ) return;

  int iret = sprintf(uuidstr, uuidFmt,
                     uuid[0], uuid[1], uuid[2], uuid[3],
                     uuid[4], uuid[5], uuid[6], uuid[7],
                     uuid[8], uuid[9], uuid[10], uuid[11],
                     uuid[12], uuid[13], uuid[14], uuid[15]);

  if ( iret != uuidNumHexChars ) uuidstr[0] = 0;
}


int str2uuid(const char *uuidstr, unsigned char *uuid)
{
  if ( uuid == NULL || uuidstr == NULL || strlen(uuidstr) != uuidNumHexChars)
    return -1;

  int iret = sscanf(uuidstr, uuidFmt,
                    &uuid[0], &uuid[1], &uuid[2], &uuid[3],
                    &uuid[4], &uuid[5], &uuid[6], &uuid[7],
                    &uuid[8], &uuid[9], &uuid[10], &uuid[11],
                    &uuid[12], &uuid[13], &uuid[14], &uuid[15]);
  if ( iret != CDI_UUID_SIZE ) return -1;
  return iret;
}

//Returns a malloc'ed string that escapes all spaces and backslashes with backslashes.
char* cdiEscapeSpaces(const char* string)
{
  //How much memory do we need?
  size_t escapeCount = 0, length = 0;
  for(; string[length]; ++length)
    escapeCount += string[length] == ' ' || string[length] == '\\';

  char* result = malloc(length + escapeCount + 1);
  if(!result) return NULL;

  //Do the escaping.
  for(size_t in = 0, out = 0; in < length; ++out, ++in)
    {
      if(string[in] == ' ' || string[in] == '\\') result[out++] = '\\';
      result[out] = string[in];
    }
  result[length + escapeCount] = 0;     //termination!
  return result;
}

//input: a space terminated string that may contain escaped characters
//output: a new zero terminated string with the escape characters removed
//*outStringEnd points to the terminating character upon return.
char* cdiUnescapeSpaces(const char* string, const char** outStringEnd)
{
  //How much memory do we need?
  size_t escapeCount = 0, length = 0;
  for(const char* current = string; *current && *current != ' '; current++)
    {
      if(*current == '\\')
        {
          current++, escapeCount++;
          if(!current) return NULL;
        }
      length++;
    }

  char* result = malloc(length + 1);
  if(!result) return NULL;

  //Do the unescaping.
  for(size_t in = 0, out = 0; out < length;)
    {
      if(string[in] == '\\') in++;
      result[out++] = string[in++];
    }
  result[length] = 0;   //termination!
  if(outStringEnd) *outStringEnd = &string[length + escapeCount];
  return result;
}

#ifdef HAVE_DECL_UUID_GENERATE
#include <sys/time.h>
#include <uuid/uuid.h>
void
create_uuid(unsigned char *uuid)
{
  static int uuid_seeded = 0;
  static char uuid_rand_state[31 * sizeof (long)];
  char *caller_rand_state;
  if (uuid_seeded)
    caller_rand_state = setstate(uuid_rand_state);
  else
    {
      struct timeval tv;
      int status = gettimeofday(&tv, NULL);
      if (status != 0)
        {
          perror("uuid random seed generation failed!");
          exit(1);
        }
      unsigned seed = (unsigned)(tv.tv_sec ^ tv.tv_usec);
      caller_rand_state = initstate(seed, uuid_rand_state,
                                    sizeof (uuid_rand_state));
      uuid_seeded = 1;
    }
  uuid_generate(uuid);
  setstate(caller_rand_state);
}
#elif defined (HAVE_DECL_UUID_CREATE)
typedef uint8_t u_int8_t;
typedef uint16_t u_int16_t;
typedef uint32_t u_int32_t;
#include <uuid.h>
void
create_uuid(unsigned char *uuid)
{
  unsigned32 status;
  uuid_create((uuid_t *)uuid, &status);
  if (status == -1)
    {
      perror("uuid generation failed!");
      exit(1);
    }
}
#else
#include <sys/time.h>
void
create_uuid(unsigned char *uuid)
{
  static int uuid_seeded = 0;
  static char uuid_rand_state[31 * sizeof (long)];
  char *caller_rand_state;
  if (uuid_seeded)
    caller_rand_state = setstate(uuid_rand_state);
  else
    {
      struct timeval tv;
      int status = gettimeofday(&tv, NULL);
      if (status != 0)
        {
          perror("failed seed generation!");
          exit(1);
        }
      unsigned seed = tv.tv_sec ^ tv.tv_usec;
      caller_rand_state = initstate(seed, uuid_rand_state,
                                    sizeof (uuid_rand_state));
      uuid_seeded = 1;
    }
  for (size_t i = 0; i < CDI_UUID_SIZE; ++i)
    uuid[i] = (unsigned char)random();
  /* encode variant into msb of octet 8 */
  uuid[8] = (unsigned char)((uuid[8] & 0x3f) | (1 << 7));
  /* encode version 4 ((pseudo-)random uuid) into msb of octet 7 */
  uuid[7] = (unsigned char)((uuid[7] & 0x0f) | (4 << 4));
  setstate(caller_rand_state);
}
#endif

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
