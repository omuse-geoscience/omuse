#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <ctype.h>
#include <stddef.h>
#include <string.h>

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"

#undef  UNDEFID
#define UNDEFID -1

/*int TableDefine = 0; */ /* Define new table also if the entry already exist */
                          /* This is needed for createtable */

#include "tablepar.h"
#include "table.h"

#define MAX_TABLE  256
#define MAX_PARS   1024

typedef struct
{
  int    used;
  PAR   *pars;
  int    npars;
  int    modelID;
  int    number;
  char  *name;
}
PARTAB;

static PARTAB parTable[MAX_TABLE];
static int  parTableSize = MAX_TABLE;
static int  parTableNum  = 0;
static int  ParTableInit = 0;

static char *tablePath = NULL;

static void tableDefModelID(int tableID, int modelID);
static void tableDefNum(int tableID, int tablenum);


void tableDefEntry(int tableID, int id, const char *name,
		   const char *longname, const char *units)
{
  int item;

  if ( tableID >= 0 && tableID < MAX_TABLE && parTable[tableID].used) { } else
    Error("Invalid table ID %d", tableID);
  item = parTable[tableID].npars++;
  parTable[tableID].pars[item].id       = id;
  parTable[tableID].pars[item].dupflags = 0;
  parTable[tableID].pars[item].name     = NULL;
  parTable[tableID].pars[item].longname = NULL;
  parTable[tableID].pars[item].units    = NULL;

  if ( name && strlen(name) > 0 )
    {
      parTable[tableID].pars[item].name     = strdupx(name);
      parTable[tableID].pars[item].dupflags |= TABLE_DUP_NAME;
    }
  if ( longname && strlen(longname) > 0 )
    {
      parTable[tableID].pars[item].longname = strdupx(longname);
      parTable[tableID].pars[item].dupflags |= TABLE_DUP_LONGNAME;
    }
  if ( units && strlen(units) > 0 )
    {
      parTable[tableID].pars[item].units    = strdupx(units);
      parTable[tableID].pars[item].dupflags |= TABLE_DUP_UNITS;
    }
}

static void tableLink(int tableID, const PAR *pars, int npars)
{
  int item;

  for ( item = 0; item < npars; item++ )
    {
      parTable[tableID].pars[item].id       = pars[item].id;
      parTable[tableID].pars[item].dupflags = 0;
      parTable[tableID].pars[item].name     = pars[item].name;
      parTable[tableID].pars[item].longname = pars[item].longname;
      parTable[tableID].pars[item].units    = pars[item].units;
    }

  parTable[tableID].npars = npars;
}

static void parTableInitEntry(int tableID)
{
  parTable[tableID].used    = 0;
  parTable[tableID].pars    = NULL;
  parTable[tableID].npars   = 0;
  parTable[tableID].modelID = UNDEFID;
  parTable[tableID].number  = UNDEFID;
  parTable[tableID].name    = NULL;
}

static void tableGetPath(void)
{
  char *path;

  path = getenv("TABLEPATH");

  if ( path ) tablePath = strdupx(path);
  /*
  printf("tablePath = %s\n", tablePath);
  */
}

static void parTableFinalize(void)
{
  for (int tableID = 0; tableID < MAX_TABLE; ++tableID)
    if (parTable[tableID].used)
      {
        int npars = parTable[tableID].npars;
        for (int item = 0; item < npars; ++item)
          {
            if (parTable[tableID].pars[item].dupflags & TABLE_DUP_NAME)
              free((void *)parTable[tableID].pars[item].name);
            if (parTable[tableID].pars[item].dupflags & TABLE_DUP_LONGNAME)
              free((void *)parTable[tableID].pars[item].longname);
            if (parTable[tableID].pars[item].dupflags & TABLE_DUP_UNITS)
              free((void *)parTable[tableID].pars[item].units);
          }
        free(parTable[tableID].pars);
        free(parTable[tableID].name);
      }
}

static void parTableInit(void)
{
  ParTableInit = 1;

  atexit(parTableFinalize);
  if ( cdiPartabIntern )
    tableDefault();

  tableGetPath();
}

static int tableNewEntry()
{
  int tableID = 0;
  static int init = 0;

  if ( ! init )
    {
      for ( tableID = 0; tableID < parTableSize; tableID++ )
	parTableInitEntry(tableID);
      init = 1;
    }

  /*
    Look for a free slot in parTable.
  */
  for ( tableID = 0; tableID < parTableSize; tableID++ )
    {
      if ( ! parTable[tableID].used ) break;
    }

  if ( tableID == parTableSize )
    Error("no more entries!");

  parTable[tableID].used = 1;
  parTableNum++;

  return (tableID);
}

static int
decodeForm1(char *pline, char *name, char *longname, char *units)
{
  char *pstart, *pend;

  /* FIXME: parse success isn't verified */
  /* long level =  */strtol(pline, &pline, 10);
  while ( isspace((int) *pline) ) pline++;

  pstart = pline;
  while ( ! (isspace((int) *pline) || *pline == 0) ) pline++;
  size_t len = (size_t)(pline - pstart);
  if ( len > 0 )
    {
      memcpy(name, pstart, len);
      name[len] = 0;
    }
  else
    return (0);

  len = strlen(pline);
  if ( len == 0 ) return (0);

  /* Format 1 : code name add mult longname [units] */
  /* FIXME: successful parse isn't verified */
  /* double add  =  */strtod(pline, &pline);
  /* FIXME: successful parse isn't verified */
  /* double mult =  */strtod(pline, &pline);

  while ( isspace((int) *pline) ) pline++;

  len = strlen(pline);
  if ( len > 0 )
    {
      pstart = pline;
      pend = strrchr(pline, '[');
      if ( pend == pstart )
        len = 0;
      else
        {
          if ( pend )
            pend--;
          else
            pend = pstart + len;
          while ( isspace((int) *pend) ) pend--;
          len = (size_t)(pend - pstart + 1);
        }
      if ( len > 0 )
	{
	  memcpy(longname, pstart, len);
	  longname[len] = 0;
	}
      pstart = strrchr(pline, '[');
      if ( pstart )
	{
	  pstart++;
	  while ( isspace((int) *pstart) ) pstart++;
	  pend = strchr(pstart, ']');
	  if ( ! pend ) return (0);
	  pend--;
	  while ( isspace((int) *pend) ) pend--;
	  len = (size_t)(pend - pstart + 1);
	  if ( len > 0 )
	    {
	      memcpy(units, pstart, len);
	      units[len] = 0;
	    }
	}
    }

  return (0);
}

static int
decodeForm2(char *pline, char *name, char *longname, char *units)
{
  /* Format 2 : code | name | longname | units */
  char *pend;
  size_t len;

  pline = strchr(pline, '|');
  pline++;

  while ( isspace((int) *pline) ) pline++;
  if (*pline != '|')
    {
      pend = strchr(pline, '|');
      if ( ! pend )
        {
          pend = pline;
          while ( ! isspace((int) *pend) ) pend++;
          len = (size_t)(pend - pline);
          if ( len > 0 )
            {
              memcpy(name, pline, len);
              name[len] = 0;
            }
          return (0);
        }
      else
        {
          pend--;
          while ( isspace((int) *pend) ) pend--;
          len = (size_t)(pend - pline + 1);
          if ( len > 0 )
            {
              memcpy(name, pline, len);
              name[len] = 0;
            }
        }
    }
  else
    name[0] = '\0';

  pline = strchr(pline, '|');
  pline++;
  while ( isspace((int) *pline) ) pline++;
  pend = strchr(pline, '|');
  if ( !pend ) pend = strchr(pline, 0);
  pend--;
  while ( isspace((int) *pend) ) pend--;
  len = (size_t)(pend - pline + 1);
  if ( len > 0 )
    {
      memcpy(longname, pline, len);
      longname[len] = 0;
    }

  pline = strchr(pline, '|');
  if ( pline )
    {
      pline++;
      while ( isspace((int) *pline) ) pline++;
      pend = strchr(pline, '|');
      if ( !pend ) pend = strchr(pline, 0);
      pend--;
      while ( isspace((int) *pend) ) pend--;
      ptrdiff_t len = pend - pline + 1;
      if ( len < 0 ) len = 0;
      memcpy(units, pline, (size_t)len);
      units[len] = 0;
    }

  return (0);
}

int tableRead(const char *tablefile)
{
  char line[1024], *pline;
  int lnr = 0;
  int id;
  char name[256], longname[256], units[256];
  int tableID = UNDEFID;
  int err;
  char *tablename;
  FILE *tablefp;

  tablefp = fopen(tablefile, "r");
  if ( tablefp == NULL ) return (tableID);

  tablename = strrchr(tablefile, '/');
  if ( tablename == 0 ) tablename = (char *) tablefile;
  else                  tablename++;

  tableID = tableDef(-1, 0, tablename);

  while ( fgets(line, 1023, tablefp) )
    {
      size_t len = strlen(line);
      if ( line[len-1] == '\n' ) line[len-1] = '\0';
      lnr++;
      id       = CDI_UNDEFID;
      name[0]     = 0;
      longname[0] = 0;
      units[0]    = 0;
      if ( line[0] == '#' ) continue;
      pline = line;

      len = strlen(pline);
      if ( len < 4 ) continue;
      while ( isspace((int) *pline) ) pline++;
      id = atoi(pline);
      /*
      if ( id > 255 ) id -= 256;
      */
      if ( id == 0 ) continue;

      while ( isdigit((int) *pline) ) pline++; 

      if ( strchr(pline, '|') )
	err = decodeForm2(pline, name, longname, units);
      else
	err = decodeForm1(pline, name, longname, units);

      if ( err ) continue;

      if ( strlen(name) == 0 ) sprintf(name, "var%d", id);

      tableDefEntry(tableID, id, name, longname, units);
    }

  return (tableID);
}

static int tableFromEnv(int modelID, int tablenum)
{
  int tableID = UNDEFID;
  char tablename[256] = {'\0'};
  int tablenamefound = 0;

  const char *modelName;
  if ( (modelName = modelInqNamePtr(modelID)) )
    {
      strcpy(tablename, modelName);
      if ( tablenum )
	{
	  size_t len = strlen(tablename);
	  sprintf(tablename+len, "_%03d", tablenum);
	}
      tablenamefound = 1;
    }
  else
    {
      int instID = modelInqInstitut(modelID);
      if ( instID != UNDEFID )
	{
          const char *instName;
	  if ( (instName = institutInqNamePtr(instID)) )
	    {
	      strcpy(tablename, instName);
	      if ( tablenum )
		{
		  size_t len = strlen(tablename);
		  sprintf(tablename+len, "_%03d", tablenum);
		}
	      tablenamefound = 1;
	    }
	}
    }

  if ( tablenamefound )
    {
      size_t lenp = 0, lenf;
      char *tablefile = NULL;
      if ( tablePath )
	lenp = strlen(tablePath);
      lenf = strlen(tablename);
      /* if (tablePath) printf("tablePath = %s\n", tablePath); */
      /* if (tablename) printf("tableName = %s\n", tablename); */
      tablefile = (char *) malloc(lenp+lenf+3);
      if ( tablePath )
	{
	  strcpy(tablefile, tablePath);
	  strcat(tablefile, "/");
	}
      else
	tablefile[0] = '\0';
      strcat(tablefile, tablename);
      /* if (tablefile) printf("tableFile = %s\n", tablefile); */

      tableID = tableRead(tablefile);
      if ( tableID != UNDEFID )
	{
	  tableDefModelID(tableID, modelID);
	  tableDefNum(tableID, tablenum);
	}
      /* printf("tableID = %d %s\n", tableID, tablefile); */

      free(tablefile);
    }

  return (tableID);
}

int tableInq(int modelID, int tablenum, const char *tablename)
{
  int tableID = UNDEFID;
  int modelID2 = UNDEFID;
  char tablefile[256] = {'\0'};

  if ( ! ParTableInit ) parTableInit();

  if ( tablename )
    {
      size_t len;
      strcpy(tablefile, tablename);
      /*
      printf("tableInq: tablefile = >%s<\n", tablefile);
      */
      /* search for internal table */
      for ( tableID = 0; tableID < MAX_TABLE; tableID++ )
	{
	  if ( parTable[tableID].used && parTable[tableID].name )
	    {
	      /* len = strlen(parTable[tableID].name); */
	      len = strlen(tablename);
	      if ( memcmp(parTable[tableID].name, tablename, len) == 0 ) break;
	    }
	}
      if ( tableID == MAX_TABLE ) tableID = UNDEFID;
      if ( CDI_Debug )
	Message("tableID = %d tablename = %s", tableID, tablename);
    }
  else
    {
      for ( tableID = 0; tableID < MAX_TABLE; tableID++ )
	{
	  if ( parTable[tableID].used )
	    {
	      if ( parTable[tableID].modelID == modelID &&
		   parTable[tableID].number  == tablenum ) break;
	    }
	}

      if ( tableID == MAX_TABLE ) tableID = UNDEFID;

      if ( tableID == UNDEFID )
	{
	  if ( modelID != UNDEFID )
	    {
              const char *modelName;
	      if ( (modelName = modelInqNamePtr(modelID)) )
		{
		  strcpy(tablefile, modelName);
		  size_t len = strlen(tablefile);
		  for ( size_t i = 0; i < len; i++)
		    if ( tablefile[i] == '.' ) tablefile[i] = '\0';
		  modelID2 = modelInq(-1, 0, tablefile);
		}
	    }
	  if ( modelID2 != UNDEFID )
	    for ( tableID = 0; tableID < MAX_TABLE; tableID++ )
	      {
		if ( parTable[tableID].used )
		  {
		    if ( parTable[tableID].modelID == modelID2 &&
			 parTable[tableID].number  == tablenum ) break;
		  }
	      }
	}

      if ( tableID == MAX_TABLE ) tableID = UNDEFID;

      if ( tableID == UNDEFID && modelID != UNDEFID )
	tableID = tableFromEnv(modelID, tablenum);

      if ( CDI_Debug )
	if ( tablename )
	  Message("tableID = %d tablename = %s", tableID, tablename);
    }

  return (tableID);
}

int tableDef(int modelID, int tablenum, const char *tablename)
{
  int tableID = UNDEFID;

  if ( ! ParTableInit ) parTableInit();
  /*
  if ( ! (modelID == UNDEFID && tablenum == 0) )
    tableID = tableInq(modelID, tablenum, tablename);
    */
  if ( tableID == UNDEFID )
    {
      tableID = tableNewEntry();

      parTable[tableID].modelID = modelID;
      parTable[tableID].number  = tablenum;
      if ( tablename )
	parTable[tableID].name = strdupx(tablename);

      parTable[tableID].pars = (PAR *) malloc(MAX_PARS * sizeof(PAR));
    }

  return (tableID);
}

static void tableDefModelID(int tableID, int modelID)
{
  parTable[tableID].modelID = modelID;
}

static void tableDefNum(int tableID, int tablenum)
{
  parTable[tableID].number  = tablenum;
}

int tableInqNum(int tableID)
{
  int number = 0;

  if ( tableID >= 0 && tableID < MAX_TABLE )
    number = parTable[tableID].number;

  return (number);
}

int tableInqModel(int tableID)
{
  int modelID = -1;

  if ( tableID >= 0 && tableID < MAX_TABLE )
    modelID = parTable[tableID].modelID;

  return (modelID);
}

static void partabCheckID(int item)
{
  if ( item < 0 || item >= parTableSize )
    Error("item %d undefined!", item);

  if ( ! parTable[item].name )
    Error("item %d name undefined!", item);
}

const char *tableInqNamePtr(int tableID)
{
  const char *tablename = NULL;

  if ( CDI_Debug )
    Message("tableID = %d", tableID);

  if ( ! ParTableInit ) parTableInit();

  if ( tableID >= 0 && tableID < parTableSize )
    if ( parTable[tableID].name )
      tablename = parTable[tableID].name;

  return (tablename);
}

void tableWrite(const char *ptfile, int tableID)
{
  int item, npars;
  size_t maxname = 4, maxlname = 10, maxunits = 2;
  FILE *ptfp;
  int tablenum, modelID, instID = CDI_UNDEFID;
  int center = 0, subcenter = 0;
  const char *instnameptr = NULL, *modelnameptr = NULL;

  if ( CDI_Debug )
    Message("write parameter table %d to %s", tableID, ptfile);

  if ( tableID == UNDEFID )
    {
      Warning("parameter table ID undefined");
      return;
    }

  partabCheckID(tableID);

  ptfp = fopen(ptfile, "w");

  npars = parTable[tableID].npars;

  for ( item = 0; item < npars; item++)
    {
      if ( parTable[tableID].pars[item].name )
	{
	  size_t lenname = strlen(parTable[tableID].pars[item].name);
	  if ( lenname  > maxname )  maxname  = lenname;
	}

      if ( parTable[tableID].pars[item].longname )
	{
	  size_t lenlname = strlen(parTable[tableID].pars[item].longname);
	  if ( lenlname > maxlname ) maxlname = lenlname;
	}

      if ( parTable[tableID].pars[item].units )
	{
	  size_t lenunits = strlen(parTable[tableID].pars[item].units);
	  if ( lenunits > maxunits ) maxunits = lenunits;
	}
    }

  tablenum = tableInqNum(tableID);
  modelID = parTable[tableID].modelID;
  if ( modelID != CDI_UNDEFID )
    {
      modelnameptr = modelInqNamePtr(modelID);
      instID = modelInqInstitut(modelID);
    }
  if ( instID != CDI_UNDEFID )
    {
      center = institutInqCenter(instID);
      subcenter = institutInqSubcenter(instID);
      instnameptr = institutInqNamePtr(instID);
    }

  fprintf(ptfp, "# Parameter table\n");
  fprintf(ptfp, "#\n");
  if ( tablenum )
    fprintf(ptfp, "# TABLE_ID=%d\n", tablenum);
  fprintf(ptfp, "# TABLE_NAME=%s\n", parTable[tableID].name);
  if ( modelnameptr )
    fprintf(ptfp, "# TABLE_MODEL=%s\n", modelnameptr);
  if ( instnameptr )
    fprintf(ptfp, "# TABLE_INSTITUT=%s\n", instnameptr);
  if ( center )
    fprintf(ptfp, "# TABLE_CENTER=%d\n", center);
  if ( subcenter )
    fprintf(ptfp, "# TABLE_SUBCENTER=%d\n", subcenter);
  fprintf(ptfp, "#\n");
  fprintf(ptfp, "#\n");
  fprintf(ptfp, "# id       = parameter ID\n");
  fprintf(ptfp, "# name     = variable name\n");
  fprintf(ptfp, "# title    = long name (description)\n");
  fprintf(ptfp, "# units    = variable units\n");
  fprintf(ptfp, "#\n");
  fprintf(ptfp, "# The format of each record is:\n");
  fprintf(ptfp, "#\n");
  fprintf(ptfp, "# id | %-*s | %-*s | %-*s\n",
	  (int)maxname,  "name",
	  (int)maxlname, "title",
	  (int)maxunits, "units");
	  
  for ( item = 0; item < npars; item++)
    {
      const char *name = parTable[tableID].pars[item].name,
        *longname = parTable[tableID].pars[item].longname,
        *units = parTable[tableID].pars[item].units;
      if ( name == NULL ) name = " ";
      if ( longname == NULL ) longname = " ";
      if ( units == NULL ) units = " ";
      fprintf(ptfp, "%4d | %-*s | %-*s | %-*s\n",
	      parTable[tableID].pars[item].id,
	      (int)maxname, name,
	      (int)maxlname, longname,
	      (int)maxunits, units);
    }

  fclose(ptfp);
}


void tableWriteC(const char *filename, int tableID)
{
  FILE *ptfp = fopen(filename, "w");
  if (!ptfp)
    Error("failed to open file \"%s\"!", filename);
  if ( CDI_Debug )
    Message("write parameter table %d to %s", tableID, filename);
  tableFWriteC(ptfp, tableID);
  fclose(ptfp);
}

void tableFWriteC(FILE *ptfp, int tableID)
{
  const char chelp[] = "";
  int item, npars;
  size_t maxname = 0, maxlname = 0, maxunits = 0;
  char tablename[256];


  if ( tableID == UNDEFID )
    {
      Warning("parameter table ID undefined");
      return;
    }

  partabCheckID(tableID);

  npars = parTable[tableID].npars;

  for ( item = 0; item < npars; item++)
    {
      if ( parTable[tableID].pars[item].name )
	{
	  size_t lenname = strlen(parTable[tableID].pars[item].name);
	  if ( lenname  > maxname )  maxname  = lenname;
	}

      if ( parTable[tableID].pars[item].longname )
	{
	  size_t lenlname = strlen(parTable[tableID].pars[item].longname);
	  if ( lenlname > maxlname ) maxlname = lenlname;
	}

      if ( parTable[tableID].pars[item].units )
	{
	  size_t lenunits = strlen(parTable[tableID].pars[item].units);
	  if ( lenunits > maxunits ) maxunits = lenunits;
	}
    }

  strncpy(tablename, parTable[tableID].name, sizeof (tablename));
  tablename[sizeof (tablename) - 1] = '\0';
  {
    size_t len = strlen(tablename);
    for (size_t i = 0; i < len; i++ )
      if ( tablename[i] == '.' ) tablename[i] = '_';
  }
  fprintf(ptfp, "static const PAR %s[] = {\n", tablename);

  for ( item = 0; item < npars; item++ )
    {
      size_t len = strlen(parTable[tableID].pars[item].name),
        llen = parTable[tableID].pars[item].longname
        ? strlen(parTable[tableID].pars[item].longname) : 0,
        ulen = parTable[tableID].pars[item].units
        ? strlen(parTable[tableID].pars[item].units) : 0;
      fprintf(ptfp, "  {%4d, 0, \"%s\", %-*s%c%s%s, %-*s%c%s%s %-*s},\n",
	      parTable[tableID].pars[item].id,
	      parTable[tableID].pars[item].name, (int)(maxname-len), chelp,
              llen?'"':' ',
              llen?parTable[tableID].pars[item].longname:"NULL",
              llen?"\"":"",
              (int)(maxlname-(llen?llen:3)), chelp,
              ulen?'"':' ',
              ulen?parTable[tableID].pars[item].units:"NULL",
              ulen?"\"":"",
              (int)(maxunits-(ulen?ulen:3)), chelp);
    }

  fprintf(ptfp, "};\n\n");
}


int tableInqParCode(int tableID, char *varname, int *code)
{
  int err = 1;

  if ( tableID != UNDEFID && varname != NULL )
    {
      int npars = parTable[tableID].npars;
      for ( int item = 0; item < npars; item++ )
	{
	  if ( parTable[tableID].pars[item].name
               && strcmp(parTable[tableID].pars[item].name, varname) == 0 )
            {
              *code = parTable[tableID].pars[item].id;
              err = 0;
              break;
            }
	}
    }

  return (err);
}


int tableInqParName(int tableID, int code, char *varname)
{
  int err = 1;

  if ( tableID != UNDEFID )
    {
      int npars = parTable[tableID].npars;
      for ( int item = 0; item < npars; item++ )
	{
	  if ( parTable[tableID].pars[item].id == code )
	    {
	      if ( parTable[tableID].pars[item].name )
		strcpy(varname, parTable[tableID].pars[item].name);     //FIXME: This may overrun the supplied buffer!
              err = 0;
	      break;
	    }
	}
    }

  return (err);
}


const char *tableInqParNamePtr(int tableID, int code)
{
  const char *name = NULL;

  if ( tableID != UNDEFID )
    {
      int npars = parTable[tableID].npars;
      for ( int item = 0; item < npars; item++ )
	{
	  if ( parTable[tableID].pars[item].id == code )
	    {
	      name = parTable[tableID].pars[item].name;
	      break;
	    }
	}
    }

  return (name);
}


const char *tableInqParLongnamePtr(int tableID, int code)
{
  const char *longname = NULL;

  if ( tableID != UNDEFID )
    {
      int npars = parTable[tableID].npars;
      for ( int item = 0; item < npars; item++ )
	{
	  if ( parTable[tableID].pars[item].id == code )
	    {
	      longname = parTable[tableID].pars[item].longname;
	      break;
	    }
	}
    }

  return (longname);
}


const char *tableInqParUnitsPtr(int tableID, int code)
{
  const char *units = NULL;

  if ( tableID != UNDEFID )
    {
      int npars = parTable[tableID].npars;
      for ( int item = 0; item < npars; item++ )
	{
	  if ( parTable[tableID].pars[item].id == code )
	    {
	      units = parTable[tableID].pars[item].units;
	      break;
	    }
	}
    }

  return (units);
}


int tableInqParLongname(int tableID, int code, char *longname)
{
  if ( ((tableID >= 0) & (tableID < MAX_TABLE)) | (tableID == UNDEFID) ) { } else
    Error("Invalid table ID %d", tableID);

  int err = 1;

  if ( tableID != UNDEFID )
    {
      int npars = parTable[tableID].npars;
      for ( int item = 0; item < npars; item++ )
	{
	  if ( parTable[tableID].pars[item].id == code )
	    {
	      if ( parTable[tableID].pars[item].longname )
		strcpy(longname, parTable[tableID].pars[item].longname);
              err = 0;
	      break;
	    }
	}
    }

  return (err);
}


int tableInqParUnits(int tableID, int code, char *units)
{

  if ( ((tableID >= 0) & (tableID < MAX_TABLE)) | (tableID == UNDEFID) ) { } else
    Error("Invalid table ID %d", tableID);

  int err = 1;

  if ( tableID != UNDEFID )
    {
      int npars = parTable[tableID].npars;
      for ( int item = 0; item < npars; item++ )
	{
	  if ( parTable[tableID].pars[item].id == code )
	    {
	      if ( parTable[tableID].pars[item].units )
		strcpy(units, parTable[tableID].pars[item].units);
              err = 0;
	      break;
	    }
	}
    }

  return (err);
}


void tableInqPar(int tableID, int code, char *name, char *longname, char *units)
{

  if ( ((tableID >= 0) & (tableID < MAX_TABLE)) | (tableID == UNDEFID) ) { } else
    Error("Invalid table ID %d", tableID);

  int npars = parTable[tableID].npars;

  for ( int item = 0; item < npars; item++ )
    {
      if ( parTable[tableID].pars[item].id == code )
	{
	  if ( parTable[tableID].pars[item].name )
	    strcpy(name, parTable[tableID].pars[item].name);
	  if ( parTable[tableID].pars[item].longname )
	    strcpy(longname, parTable[tableID].pars[item].longname);
	  if ( parTable[tableID].pars[item].units )
	    strcpy(units, parTable[tableID].pars[item].units);
	  break;
	}
    }
}

int tableInqNumber(void)
{
  if ( ! ParTableInit ) parTableInit();

  return (parTableNum);
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
