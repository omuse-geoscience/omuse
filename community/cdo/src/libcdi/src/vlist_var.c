#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include <limits.h>

#include "dmemory.h"
#include "cdi.h"
#include "cdi_int.h"
#include "vlist.h"
#include "vlist_var.h"
#include "resource_handle.h"
#include "vlist_att.h"
#include "namespace.h"
#include "serialize.h"
#include "error.h"

#if  defined  (HAVE_LIBGRIB_API)
#  include "file.h"

#  include <grib_api.h>
#endif


static
void vlistvarInitEntry(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistptr->vars[varID].fvarID        = varID;
  vlistptr->vars[varID].mvarID        = varID;
  vlistptr->vars[varID].flag          = 0;
  vlistptr->vars[varID].param         = 0;
  vlistptr->vars[varID].datatype      = CDI_UNDEFID;
  vlistptr->vars[varID].tsteptype     = TSTEP_INSTANT;
  vlistptr->vars[varID].timave        = 0;
  vlistptr->vars[varID].timaccu       = 0;
  vlistptr->vars[varID].typeOfGeneratingProcess   = 0;
  vlistptr->vars[varID].productDefinitionTemplate = -1;
  vlistptr->vars[varID].chunktype     = cdiChunkType;
  vlistptr->vars[varID].xyz           = 321;
  vlistptr->vars[varID].gridID        = CDI_UNDEFID;
  vlistptr->vars[varID].zaxisID       = CDI_UNDEFID;
  vlistptr->vars[varID].subtypeID     = CDI_UNDEFID;
  vlistptr->vars[varID].instID        = CDI_UNDEFID;
  vlistptr->vars[varID].modelID       = CDI_UNDEFID;
  vlistptr->vars[varID].tableID       = CDI_UNDEFID;
  vlistptr->vars[varID].missvalused   = FALSE;
  vlistptr->vars[varID].missval       = cdiDefaultMissval;
  vlistptr->vars[varID].addoffset     = 0.0;
  vlistptr->vars[varID].scalefactor   = 1.0;
  vlistptr->vars[varID].name          = NULL;
  vlistptr->vars[varID].longname      = NULL;
  vlistptr->vars[varID].stdname       = NULL;
  vlistptr->vars[varID].units         = NULL;
  vlistptr->vars[varID].extra         = NULL;
  vlistptr->vars[varID].levinfo       = NULL;
  vlistptr->vars[varID].comptype      = COMPRESS_NONE;
  vlistptr->vars[varID].complevel     = 1;
  vlistptr->vars[varID].atts.nalloc   = MAX_ATTRIBUTES;
  vlistptr->vars[varID].atts.nelems   = 0;
  vlistptr->vars[varID].lvalidrange   = 0;
  vlistptr->vars[varID].validrange[0] = VALIDMISS;
  vlistptr->vars[varID].validrange[1] = VALIDMISS;
  vlistptr->vars[varID].ensdata       = NULL;
  vlistptr->vars[varID].iorank        = CDI_UNDEFID;
  vlistptr->vars[varID].opt_grib_kvpair_size = 0;
  vlistptr->vars[varID].opt_grib_kvpair      = NULL;
  vlistptr->vars[varID].opt_grib_nentries    = 0;
}



static
int vlistvarNewEntry(int vlistID)
{
  int varID = 0;
  int vlistvarSize;
  var_t *vlistvar;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistvarSize = vlistptr->varsAllocated;
  vlistvar     = vlistptr->vars;
  /*
    Look for a free slot in vlistvar.
    (Create the table the first time through).
  */
  if ( ! vlistvarSize )
    {
      vlistvarSize = 2;
      vlistvar = (var_t *)xmalloc((size_t)vlistvarSize * sizeof (var_t));
      for (int i = 0; i < vlistvarSize; i++ )
	vlistvar[i].isUsed = FALSE;
    }
  else
    {
      while (varID < vlistvarSize && vlistvar[varID].isUsed)
        ++varID;
    }
  /*
    If the table overflows, double its size.
  */
  if ( varID == vlistvarSize )
    {
      int i;

      vlistvarSize = 2*vlistvarSize;
      vlistvar = (var_t *)xrealloc(vlistvar,
                                   (size_t)vlistvarSize * sizeof(var_t));
      varID = vlistvarSize/2;

      for ( i = varID; i < vlistvarSize; i++ )
	vlistvar[i].isUsed = FALSE;
    }

  vlistptr->varsAllocated = vlistvarSize;
  vlistptr->vars          = vlistvar;

  vlistvarInitEntry(vlistID, varID);

  vlistptr->vars[varID].isUsed = TRUE;

  return (varID);
}

void vlistCheckVarID(const char *caller, int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if ( vlistptr == NULL )
    Errorc("vlist undefined!");

  if ( varID < 0 || varID >= vlistptr->nvars )
    Errorc("varID %d undefined!", varID);

  if ( ! vlistptr->vars[varID].isUsed )
    Errorc("varID %d undefined!", varID);
}


int vlistDefVarTiles(int vlistID, int gridID, int zaxisID, int tsteptype, int tilesetID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  if ( CDI_Debug )
    Message("gridID = %d  zaxisID = %d  tsteptype = %d", gridID, zaxisID, tsteptype);

  int varID = vlistvarNewEntry(vlistID);

  vlistptr->nvars++;
  vlistptr->vars[varID].gridID    = gridID;
  vlistptr->vars[varID].zaxisID   = zaxisID;
  vlistptr->vars[varID].tsteptype = tsteptype;
  vlistptr->vars[varID].subtypeID = tilesetID;

  if ( tsteptype < 0 )
    {
      Message("Unexpected tstep type %d, set to TSTEP_INSTANT!", tsteptype);
      vlistptr->vars[varID].tsteptype = TSTEP_INSTANT;
    }

  vlistAdd2GridIDs(vlistptr, gridID);
  vlistAdd2ZaxisIDs(vlistptr, zaxisID);
  vlistAdd2SubtypeIDs(vlistptr, tilesetID);

  vlistptr->vars[varID].param = cdiEncodeParam(-(varID + 1), 255, 255);
  reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
  return (varID);
}

/*
@Function  vlistDefVar
@Title     Define a Variable

@Prototype int vlistDefVar(int vlistID, int gridID, int zaxisID, int tsteptype)
@Parameter
    @Item  vlistID   Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  gridID    Grid ID, from a previous call to @fref{gridCreate}.
    @Item  zaxisID   Z-axis ID, from a previous call to @fref{zaxisCreate}.
    @Item  tsteptype One of the set of predefined CDI timestep types.
                     The valid CDI timestep types are @func{TSTEP_CONSTANT} and @func{TSTEP_INSTANT}.

@Description
The function @func{vlistDefVar} adds a new variable to vlistID.

@Result
@func{vlistDefVar} returns an identifier to the new variable.

@Example
Here is an example using @func{vlistCreate} to create a variable list
and add a variable with @func{vlistDefVar}.

@Source
#include "cdi.h"
   ...
int vlistID, varID;
   ...
vlistID = vlistCreate();
varID = vlistDefVar(vlistID, gridID, zaxisID, TIME_INSTANT);
   ...
streamDefVlist(streamID, vlistID);
   ...
vlistDestroy(vlistID);
   ...
@EndSource
@EndFunction
*/
int vlistDefVar(int vlistID, int gridID, int zaxisID, int tsteptype)
{
  /* call "vlistDefVarTiles" with a trivial tile index: */
  return vlistDefVarTiles(vlistID, gridID, zaxisID, tsteptype, CDI_UNDEFID);
}

void
cdiVlistCreateVarLevInfo(vlist_t *vlistptr, int varID)
{
  xassert(varID >= 0 && varID < vlistptr->nvars
          && vlistptr->vars[varID].levinfo == NULL);
  int zaxisID = vlistptr->vars[varID].zaxisID;
  size_t nlevs = (size_t)zaxisInqSize(zaxisID);

  vlistptr->vars[varID].levinfo
    = (levinfo_t*)xmalloc((size_t)nlevs * sizeof(levinfo_t));

  for (size_t levID = 0; levID < nlevs; levID++ )
    vlistptr->vars[varID].levinfo[levID] = DEFAULT_LEVINFO((int)levID);
}

/*
@Function  vlistDefVarParam
@Title     Define the parameter number of a Variable

@Prototype void vlistDefVarParam(int vlistID, int varID, int param)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  param    Parameter number.

@Description
The function @func{vlistDefVarParam} defines the parameter number of a variable.

@EndFunction
*/
void vlistDefVarParam(int vlistID, int varID, int param)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if (vlistptr->vars[varID].param != param)
    {
      vlistptr->vars[varID].param = param;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

/*
@Function  vlistDefVarCode
@Title     Define the code number of a Variable

@Prototype void vlistDefVarCode(int vlistID, int varID, int code)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  code     Code number.

@Description
The function @func{vlistDefVarCode} defines the code number of a variable.

@EndFunction
*/
void vlistDefVarCode(int vlistID, int varID, int code)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  int param = vlistptr->vars[varID].param;
  int pnum, pcat, pdis;
  cdiDecodeParam(param, &pnum, &pcat, &pdis);
  int newParam = cdiEncodeParam(code, pcat, pdis);
  if (vlistptr->vars[varID].param != newParam)
    {
      vlistptr->vars[varID].param = newParam;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


void vlistInqVar(int vlistID, int varID, int *gridID, int *zaxisID, int *tsteptype)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  *gridID    = vlistptr->vars[varID].gridID;
  *zaxisID   = vlistptr->vars[varID].zaxisID;
  *tsteptype = vlistptr->vars[varID].tsteptype;

  return;
}

/*
@Function  vlistInqVarGrid
@Title     Get the Grid ID of a Variable

@Prototype int vlistInqVarGrid(int vlistID, int varID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier.

@Description
The function @func{vlistInqVarGrid} returns the grid ID of a variable.

@Result
@func{vlistInqVarGrid} returns the grid ID of the variable.

@EndFunction
*/
int vlistInqVarGrid(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].gridID);
}

/*
@Function  vlistInqVarZaxis
@Title     Get the Zaxis ID of a Variable

@Prototype int vlistInqVarZaxis(int vlistID, int varID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier.

@Description
The function @func{vlistInqVarZaxis} returns the zaxis ID of a variable.

@Result
@func{vlistInqVarZaxis} returns the zaxis ID of the variable.

@EndFunction
*/
int vlistInqVarZaxis(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].zaxisID);
}


/*
@Function  vlistInqVarSubtype
@Title     Get the Subtype ID of a Variable

@Description
The function @func{vlistInqVarSubtype} returns the subtype ID of a variable.

@Result
@func{vlistInqVarSubtype} returns the subtype ID of the variable.

@EndFunction
*/
int vlistInqVarSubtype(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);
  return (vlistptr->vars[varID].subtypeID);
}


/*
@Function  vlistInqVarParam
@Title     Get the parameter number of a Variable

@Prototype int vlistInqVarParam(int vlistID, int varID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier.

@Description
The function @func{vlistInqVarParam} returns the parameter number of a variable.

@Result
@func{vlistInqVarParam} returns the parameter number of the variable.

@EndFunction
*/
int vlistInqVarParam(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].param);
}

/*
@Function  vlistInqVarCode
@Title     Get the Code number of a Variable

@Prototype int vlistInqVarCode(int vlistID, int varID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier.

@Description
The function @func{vlistInqVarCode} returns the code number of a variable.

@Result
@func{vlistInqVarCode} returns the code number of the variable.

@EndFunction
*/
int vlistInqVarCode(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  int param = vlistptr->vars[varID].param;
  int pdis, pcat, pnum;
  cdiDecodeParam(param, &pnum, &pcat, &pdis);
  int code = pnum;
  if ( pdis != 255 ) code = -varID-1; // GRIB2 Parameter

  if ( code < 0 && vlistptr->vars[varID].tableID != -1 && vlistptr->vars[varID].name != NULL )
    {
      tableInqParCode(vlistptr->vars[varID].tableID, vlistptr->vars[varID].name, &code);
    }

  return (code);
}


const char *vlistInqVarNamePtr(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].name);
}


const char *vlistInqVarLongnamePtr(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].longname);
}


const char *vlistInqVarStdnamePtr(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].stdname);
}


const char *vlistInqVarUnitsPtr(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].units);
}

/*
@Function  vlistInqVarName
@Title     Get the name of a Variable

@Prototype void vlistInqVarName(int vlistID, int varID, char *name)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier.
    @Item  name     Returned variable name. The caller must allocate space for the
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{vlistInqVarName} returns the name of a variable.

@Result
@func{vlistInqVarName} returns the name of the variable to the parameter name if available,
otherwise the result is an empty string.

@EndFunction
*/
void vlistInqVarName(int vlistID, int varID, char *name)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if ( vlistptr->vars[varID].name == NULL )
    {
      int param = vlistptr->vars[varID].param;
      int pdis, pcat, pnum;
      cdiDecodeParam(param, &pnum, &pcat, &pdis);
      if ( pdis == 255 )
	{
	  int code = pnum;
	  int tableID = vlistptr->vars[varID].tableID;
	  if ( tableInqParName(tableID, code, name) != 0 )
	    sprintf(name, "var%d", code);
	}
      else
	{
	  sprintf(name, "param%d.%d.%d", pnum, pcat, pdis);
	}
    }
  else
    strcpy(name, vlistptr->vars[varID].name);   //FIXME: This may overrun the provided buffer.

  return;
}

/*
@Function vlistCopyVarName
@Tatle    Get the name of a Variable in a safe way

@Prototype char* vlistCopyVarName(int vlistId, int varId)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier.

@Return A pointer to a malloc'ed string. Must be cleaned up with free().

@Description
This is the buffer overflow immune version of vlistInqVarName().
The memory for the returned string is allocated to fit the string via malloc().

@EndFunction
*/
char* vlistCopyVarName(int vlistId, int varId)
{
  vlist_t* vlistptr = vlist_to_pointer(vlistId);
  vlistCheckVarID(__func__, vlistId, varId);

  //If a name is set in the variable description, use that.
  const char* name = vlistptr->vars[varId].name;
  if(name) return strdup(name);

  //Otherwise we check if we should use the table of parameter descriptions.
  int param = vlistptr->vars[varId].param;
  int discipline, category, number;
  cdiDecodeParam(param, &number, &category, &discipline);
  char *result = NULL;
  if(discipline == 255)
    {
      int tableId = vlistptr->vars[varId].tableID;
      if(( name = tableInqParNamePtr(tableId, number) ))
        result = strdup(name);
      {
        //No luck, fall back to outputting a name of the format "var<num>".
        result = xmalloc(3 + 3 * sizeof (int) * CHAR_BIT / 8 + 2);
        sprintf(result, "var%d", number);
      }
    }
  else
    {
      result = xmalloc(5 + 2 + 3 * (3 * sizeof (int) * CHAR_BIT + 1) + 1);
      sprintf(result, "param%d.%d.%d", number, category, discipline);
    }
  //Finally, we fall back to outputting a name of the format "param<num>.<cat>.<dis>".
  return result;
}

/*
@Function  vlistInqVarLongname
@Title     Get the longname of a Variable

@Prototype void vlistInqVarLongname(int vlistID, int varID, char *longname)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier.
    @Item  longname Long name of the variable. The caller must allocate space for the
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{vlistInqVarLongname} returns the longname of a variable if available,
otherwise the result is an empty string.

@Result
@func{vlistInqVaeLongname} returns the longname of the variable to the parameter longname.

@EndFunction
*/
void vlistInqVarLongname(int vlistID, int varID, char *longname)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  longname[0] = '\0';

  if ( vlistptr->vars[varID].longname == NULL )
    {
      int param = vlistptr->vars[varID].param;
      int pdis, pcat, pnum;
      cdiDecodeParam(param, &pnum, &pcat, &pdis);
      if ( pdis == 255 )
	{
	  int code = pnum;
          int tableID = vlistptr->vars[varID].tableID;
	  if ( tableInqParLongname(tableID, code, longname) != 0 )
	    longname[0] = '\0';
	}
    }
  else
    strcpy(longname, vlistptr->vars[varID].longname);

  return;
}

/*
@Function  vlistInqVarStdname
@Title     Get the standard name of a Variable

@Prototype void vlistInqVarStdname(int vlistID, int varID, char *stdname)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier.
    @Item  stdname  Standard name of the variable. The caller must allocate space for the
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{vlistInqVarStdname} returns the standard name of a variable if available,
otherwise the result is an empty string.

@Result
@func{vlistInqVarName} returns the standard name of the variable to the parameter stdname.

@EndFunction
*/
void vlistInqVarStdname(int vlistID, int varID, char *stdname)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if ( vlistptr->vars[varID].stdname == NULL )
    {
      stdname[0] = '\0';
    }
  else
    strcpy(stdname, vlistptr->vars[varID].stdname);

  return;
}

/*
@Function  vlistInqVarUnits
@Title     Get the units of a Variable

@Prototype void vlistInqVarUnits(int vlistID, int varID, char *units)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier.
    @Item  units    Units of the variable. The caller must allocate space for the
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{vlistInqVarUnits} returns the units of a variable if available,
otherwise the result is an empty string.

@Result
@func{vlistInqVarUnits} returns the units of the variable to the parameter units.

@EndFunction
*/
void vlistInqVarUnits(int vlistID, int varID, char *units)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  units[0] = '\0';

  if ( vlistptr->vars[varID].units == NULL )
    {
      int param = vlistptr->vars[varID].param;
      int pdis, pcat, pnum;
      cdiDecodeParam(param, &pnum, &pcat, &pdis);
      if ( pdis == 255 )
	{
	  int code = pnum;
	  int tableID = vlistptr->vars[varID].tableID;
	  if ( tableInqParUnits(tableID, code, units) != 0 )
	    units[0] = '\0';
	}
    }
  else
    strcpy(units, vlistptr->vars[varID].units);

  return;
}

/* used in MPIOM ! */
int vlistInqVarID(int vlistID, int code)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  for ( int varID = 0; varID < vlistptr->nvars; varID++ )
    {
      int param = vlistptr->vars[varID].param;
      int pdis, pcat, pnum;
      cdiDecodeParam(param, &pnum, &pcat, &pdis);
      if ( pnum == code ) return (varID);
    }

  return (CDI_UNDEFID);
}


int vlistInqVarSize(int vlistID, int varID)
{
  vlistCheckVarID(__func__, vlistID, varID);

  int zaxisID, gridID;
  int tsteptype;
  vlistInqVar(vlistID, varID, &gridID, &zaxisID, &tsteptype);

  int nlevs = zaxisInqSize(zaxisID);

  int gridsize = gridInqSize(gridID);

  int size = gridsize*nlevs;

  return (size);
}

/*
@Function  vlistInqVarDatatype
@Title     Get the data type of a Variable

@Prototype int vlistInqVarDatatype(int vlistID, int varID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier.

@Description
The function @func{vlistInqVarDatatype} returns the data type of a variable.

@Result
@func{vlistInqVarDatatype} returns an identifier to the data type of the variable.
The valid CDI data types are @func{DATATYPE_PACK8}, @func{DATATYPE_PACK16}, @func{DATATYPE_PACK24},
@func{DATATYPE_FLT32}, @func{DATATYPE_FLT64}, @func{DATATYPE_INT8}, @func{DATATYPE_INT16} and 
@func{DATATYPE_INT32}.

@EndFunction
*/
int vlistInqVarDatatype(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].datatype);
}


int vlistInqVarNumber(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  int number = CDI_REAL;
  if ( vlistptr->vars[varID].datatype == DATATYPE_CPX32 ||
       vlistptr->vars[varID].datatype == DATATYPE_CPX64 )
    number = CDI_COMP;

  return (number);
}

/*
@Function  vlistDefVarDatatype
@Title     Define the data type of a Variable

@Prototype void vlistDefVarDatatype(int vlistID, int varID, int datatype)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  datatype The data type identifier.
                    The valid CDI data types are @func{DATATYPE_PACK8}, @func{DATATYPE_PACK16},
                    @func{DATATYPE_PACK24}, @func{DATATYPE_FLT32}, @func{DATATYPE_FLT64},
                    @func{DATATYPE_INT8}, @func{DATATYPE_INT16} and @func{DATATYPE_INT32}.

@Description
The function @func{vlistDefVarDatatype} defines the data type of a variable.

@EndFunction
*/
void vlistDefVarDatatype(int vlistID, int varID, int datatype)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if (vlistptr->vars[varID].datatype != datatype)
    {
      vlistptr->vars[varID].datatype = datatype;

      if ( vlistptr->vars[varID].missvalused == FALSE )
        switch (datatype)
          {
          case DATATYPE_INT8:   vlistptr->vars[varID].missval = -SCHAR_MAX; break;
          case DATATYPE_UINT8:  vlistptr->vars[varID].missval =  UCHAR_MAX; break;
          case DATATYPE_INT16:  vlistptr->vars[varID].missval = -SHRT_MAX;  break;
          case DATATYPE_UINT16: vlistptr->vars[varID].missval =  USHRT_MAX; break;
          case DATATYPE_INT32:  vlistptr->vars[varID].missval = -INT_MAX;   break;
          case DATATYPE_UINT32: vlistptr->vars[varID].missval =  UINT_MAX;  break;
          }
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


void vlistDefVarInstitut(int vlistID, int varID, int instID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  if (vlistptr->vars[varID].instID != instID)
    {
      vlistptr->vars[varID].instID = instID;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


int vlistInqVarInstitut(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  return (vlistptr->vars[varID].instID);
}


void vlistDefVarModel(int vlistID, int varID, int modelID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  if (vlistptr->vars[varID].modelID != modelID)
    {
      vlistptr->vars[varID].modelID = modelID;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


int vlistInqVarModel(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  return (vlistptr->vars[varID].modelID);
}


void vlistDefVarTable(int vlistID, int varID, int tableID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if (vlistptr->vars[varID].tableID != tableID)
    {
      vlistptr->vars[varID].tableID = tableID;
      int tablenum = tableInqNum(tableID);

      int param = vlistptr->vars[varID].param;

      int pnum, pcat, pdis;
      cdiDecodeParam(param, &pnum, &pcat, &pdis);
      vlistptr->vars[varID].param = cdiEncodeParam(pnum, tablenum, pdis);
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


int vlistInqVarTable(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  return (vlistptr->vars[varID].tableID);
}

/*
@Function  vlistDefVarName
@Title     Define the name of a Variable

@Prototype void vlistDefVarName(int vlistID, int varID, const char *name)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  name     Name of the variable.

@Description
The function @func{vlistDefVarName} defines the name of a variable.

@EndFunction
*/
void vlistDefVarName(int vlistID, int varID, const char *name)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if ( name )
    {
      if ( vlistptr->vars[varID].name )
	{
	  free(vlistptr->vars[varID].name);
	  vlistptr->vars[varID].name = NULL;
	}

      vlistptr->vars[varID].name = strdupx(name);
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

/*
@Function  vlistDefVarLongname
@Title     Define the long name of a Variable

@Prototype void vlistDefVarLongname(int vlistID, int varID, const char *longname)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  longname Long name of the variable.

@Description
The function @func{vlistDefVarLongname} defines the long name of a variable.

@EndFunction
*/
void vlistDefVarLongname(int vlistID, int varID, const char *longname)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if ( longname )
    {
      if ( vlistptr->vars[varID].longname )
	{
	  free(vlistptr->vars[varID].longname);
	  vlistptr->vars[varID].longname = 0;
	}

      vlistptr->vars[varID].longname = strdupx(longname);
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

/*
@Function  vlistDefVarStdname
@Title     Define the standard name of a Variable

@Prototype void vlistDefVarStdname(int vlistID, int varID, const char *stdname)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  stdname  Standard name of the variable.

@Description
The function @func{vlistDefVarStdname} defines the standard name of a variable.

@EndFunction
*/
void vlistDefVarStdname(int vlistID, int varID, const char *stdname)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if ( stdname )
    {
      if ( vlistptr->vars[varID].stdname )
	{
	  free(vlistptr->vars[varID].stdname);
	  vlistptr->vars[varID].stdname = 0;
	}

      vlistptr->vars[varID].stdname = strdupx(stdname);
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

/*
@Function  vlistDefVarUnits
@Title     Define the units of a Variable

@Prototype void vlistDefVarUnits(int vlistID, int varID, const char *units)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  units    Units of the variable.

@Description
The function @func{vlistDefVarUnits} defines the units of a variable.

@EndFunction
*/
void vlistDefVarUnits(int vlistID, int varID, const char *units)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if ( units )
    {
      if ( vlistptr->vars[varID].units )
	{
	  free(vlistptr->vars[varID].units);
	  vlistptr->vars[varID].units = 0;
	}

      vlistptr->vars[varID].units = strdupx(units);
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

/*
@Function  vlistInqVarMissval
@Title     Get the missing value of a Variable

@Prototype double vlistInqVarMissval(int vlistID, int varID)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier.

@Description
The function @func{vlistInqVarMissval} returns the missing value of a variable.

@Result
@func{vlistInqVarMissval} returns the missing value of the variable.

@EndFunction
*/
double vlistInqVarMissval(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].missval);
}

/*
@Function  vlistDefVarMissval
@Title     Define the missing value of a Variable

@Prototype void vlistDefVarMissval(int vlistID, int varID, double missval)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  missval  Missing value.

@Description
The function @func{vlistDefVarMissval} defines the missing value of a variable.

@EndFunction
*/
void vlistDefVarMissval(int vlistID, int varID, double missval)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  vlistptr->vars[varID].missval = missval;
  vlistptr->vars[varID].missvalused = TRUE;
}

/*
@Function  vlistDefVarExtra
@Title     Define extra information of a Variable

@Prototype void vlistDefVarExtra(int vlistID, int varID, const char *extra)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate}.
    @Item  varID    Variable identifier.
    @Item  extra    Extra information.

@Description
The function @func{vlistDefVarExtra} defines the extra information of a variable.

@EndFunction
*/
void vlistDefVarExtra(int vlistID, int varID, const char *extra)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if ( extra )
    {
      if ( vlistptr->vars[varID].extra )
	{
	  free(vlistptr->vars[varID].extra);
	  vlistptr->vars[varID].extra = NULL;
	}

      vlistptr->vars[varID].extra = strdupx(extra);
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}

/*
@Function  vlistInqVarExtra
@Title     Get extra information of a Variable

@Prototype void vlistInqVarExtra(int vlistID, int varID, char *extra)
@Parameter
    @Item  vlistID  Variable list ID, from a previous call to @fref{vlistCreate} or @fref{streamInqVlist}.
    @Item  varID    Variable identifier.
    @Item  extra    Returned variable extra information. The caller must allocate space for the
                    returned string. The maximum possible length, in characters, of
                    the string is given by the predefined constant @func{CDI_MAX_NAME}.

@Description
The function @func{vlistInqVarExtra} returns the extra information of a variable.

@Result
@func{vlistInqVarExtra} returns the extra information of the variable to the parameter extra if available,
otherwise the result is an empty string.

@EndFunction
*/
void vlistInqVarExtra(int vlistID, int varID, char *extra)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if ( vlistptr->vars[varID].extra == NULL )
      sprintf(extra, "-");
  else
    strcpy(extra, vlistptr->vars[varID].extra);

  return;
}


int vlistInqVarValidrange(int vlistID, int varID, double *validrange)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if ( validrange != NULL && vlistptr->vars[varID].lvalidrange )
    {
      validrange[0] = vlistptr->vars[varID].validrange[0];
      validrange[1] = vlistptr->vars[varID].validrange[1];
    }

  return (vlistptr->vars[varID].lvalidrange);
}


void vlistDefVarValidrange(int vlistID, int varID, const double *validrange)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  vlistptr->vars[varID].validrange[0] = validrange[0];
  vlistptr->vars[varID].validrange[1] = validrange[1];
  vlistptr->vars[varID].lvalidrange = TRUE;
  reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
}


double vlistInqVarScalefactor(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].scalefactor);
}


double vlistInqVarAddoffset(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].addoffset);
}


void vlistDefVarScalefactor(int vlistID, int varID, double scalefactor)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if ( IS_NOT_EQUAL(vlistptr->vars[varID].scalefactor, scalefactor) )
    {
      vlistptr->vars[varID].scalefactor = scalefactor;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


void vlistDefVarAddoffset(int vlistID, int varID, double addoffset)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if ( IS_NOT_EQUAL(vlistptr->vars[varID].addoffset, addoffset))
    {
      vlistptr->vars[varID].addoffset = addoffset;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


void vlistDefVarTsteptype(int vlistID, int varID, int tsteptype)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  if (vlistptr->vars[varID].tsteptype != tsteptype)
    {
      vlistptr->vars[varID].tsteptype = tsteptype;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


int vlistInqVarTsteptype(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  return (vlistptr->vars[varID].tsteptype);
}


void vlistDefVarTimave(int vlistID, int varID, int timave)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  if (vlistptr->vars[varID].timave != timave)
    {
      vlistptr->vars[varID].timave = timave;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


int vlistInqVarTimave(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  return (vlistptr->vars[varID].timave);
}


void vlistDefVarTimaccu(int vlistID, int varID, int timaccu)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  if (vlistptr->vars[varID].timaccu != timaccu)
    {
      vlistptr->vars[varID].timaccu = timaccu;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


int vlistInqVarTimaccu(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  return (vlistptr->vars[varID].timaccu);
}


void vlistDefVarTypeOfGeneratingProcess(int vlistID, int varID, int typeOfGeneratingProcess)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  if (vlistptr->vars[varID].typeOfGeneratingProcess != typeOfGeneratingProcess)
    {
      vlistptr->vars[varID].typeOfGeneratingProcess = typeOfGeneratingProcess;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


int vlistInqVarTypeOfGeneratingProcess(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  return (vlistptr->vars[varID].typeOfGeneratingProcess);
}


void vlistDefVarProductDefinitionTemplate(int vlistID, int varID, int productDefinitionTemplate)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if (vlistptr->vars[varID].productDefinitionTemplate != productDefinitionTemplate)
    {
      vlistptr->vars[varID].productDefinitionTemplate = productDefinitionTemplate;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


int vlistInqVarProductDefinitionTemplate(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  return (vlistptr->vars[varID].productDefinitionTemplate);
}


void vlistDestroyVarName(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  if ( vlistptr->vars[varID].name )
    {
      free(vlistptr->vars[varID].name);
      vlistptr->vars[varID].name = NULL;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


void vlistDestroyVarLongname(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if ( vlistptr->vars[varID].longname )
    {
      free(vlistptr->vars[varID].longname);
      vlistptr->vars[varID].longname = NULL;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


void vlistDestroyVarStdname(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if ( vlistptr->vars[varID].stdname )
    {
      free(vlistptr->vars[varID].stdname);
      vlistptr->vars[varID].stdname = NULL;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


void vlistDestroyVarUnits(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if ( vlistptr->vars[varID].units )
    {
      free(vlistptr->vars[varID].units);
      vlistptr->vars[varID].units = NULL;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


int vlistInqVarMissvalUsed(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  return (vlistptr->vars[varID].missvalused);
}


void vlistDefFlag(int vlistID, int varID, int levID, int flag)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  levinfo_t li = DEFAULT_LEVINFO(levID);
  if (vlistptr->vars[varID].levinfo)
    ;
  else if (flag != li.flag)
    cdiVlistCreateVarLevInfo(vlistptr, varID);
  else
    return;

  vlistptr->vars[varID].levinfo[levID].flag = flag;

  vlistptr->vars[varID].flag = 0;

  int nlevs = zaxisInqSize(vlistptr->vars[varID].zaxisID);
  for ( int levelID = 0; levelID < nlevs; levelID++ )
    {
      if ( vlistptr->vars[varID].levinfo[levelID].flag )
        {
          vlistptr->vars[varID].flag = 1;
          break;
        }
    }

  reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
}


int vlistInqFlag(int vlistID, int varID, int levID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if (vlistptr->vars[varID].levinfo)
    return (vlistptr->vars[varID].levinfo[levID].flag);
  else
    {
      levinfo_t li = DEFAULT_LEVINFO(levID);
      return li.flag;
    }
}


int vlistFindVar(int vlistID, int fvarID)
{
  int varID;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  for ( varID = 0; varID < vlistptr->nvars; varID++ )
    {
      if ( vlistptr->vars[varID].fvarID == fvarID ) break;
    }

  if ( varID == vlistptr->nvars )
    {
      varID = -1;
      Message("varID not found for fvarID %d in vlistID %d!", fvarID, vlistID);
    }

  return (varID);
}


int vlistFindLevel(int vlistID, int fvarID, int flevelID)
{
  int levelID = -1;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  int varID = vlistFindVar(vlistID, fvarID);

  if ( varID != -1 )
    {
      int nlevs = zaxisInqSize(vlistptr->vars[varID].zaxisID);
      for ( levelID = 0; levelID < nlevs; levelID++ )
	{
	  if ( vlistptr->vars[varID].levinfo[levelID].flevelID == flevelID ) break;
	}

      if ( levelID == nlevs )
	{
	  levelID = -1;
	  Message("levelID not found for fvarID %d and levelID %d in vlistID %d!",
		  fvarID, flevelID, vlistID);
	}
    }

  return (levelID);
}


int vlistMergedVar(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  return (vlistptr->vars[varID].mvarID);
}


int vlistMergedLevel(int vlistID, int varID, int levelID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if (vlistptr->vars[varID].levinfo)
    return vlistptr->vars[varID].levinfo[levelID].mlevelID;
  else
    {
      levinfo_t li = DEFAULT_LEVINFO(levelID);
      return li.mlevelID;
    }
}


void vlistDefIndex(int vlistID, int varID, int levelID, int index)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  levinfo_t li = DEFAULT_LEVINFO(levelID);
  if (vlistptr->vars[varID].levinfo)
    ;
  else if (index != li.index)
    cdiVlistCreateVarLevInfo(vlistptr, varID);
  else
    return;
  vlistptr->vars[varID].levinfo[levelID].index = index;
  reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
}


int vlistInqIndex(int vlistID, int varID, int levelID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  if (vlistptr->vars[varID].levinfo)
    return (vlistptr->vars[varID].levinfo[levelID].index);
  else
    {
      levinfo_t li = DEFAULT_LEVINFO(levelID);
      return li.index;
    }
}


void vlistChangeVarZaxis(int vlistID, int varID, int zaxisID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  int nlevs1 = zaxisInqSize(vlistptr->vars[varID].zaxisID);
  int nlevs2 = zaxisInqSize(zaxisID);

  if ( nlevs1 != nlevs2 ) Error("Number of levels must not change!");

  int nvars = vlistptr->nvars;
  int found = 0;
  int oldZaxisID = vlistptr->vars[varID].zaxisID;
  for ( int i = 0; i < varID; ++i)
    found |= (vlistptr->vars[i].zaxisID == oldZaxisID);
  for ( int i = varID + 1; i < nvars; ++i)
    found |= (vlistptr->vars[i].zaxisID == oldZaxisID);

  if (found)
    {
      int nzaxis = vlistptr->nzaxis;
      for (int i = 0; i < nzaxis; ++i)
	if (vlistptr->zaxisIDs[i] == oldZaxisID )
	  vlistptr->zaxisIDs[i] = zaxisID;
    }
  else
    vlistAdd2ZaxisIDs(vlistptr, zaxisID);

  vlistptr->vars[varID].zaxisID = zaxisID;
  reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
}


void vlistChangeVarGrid(int vlistID, int varID, int gridID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  int nvars = vlistptr->nvars;
  int index;
  for ( index = 0; index < nvars; index++ )
    if ( index != varID )
      if ( vlistptr->vars[index].gridID == vlistptr->vars[varID].gridID ) break;

  if ( index == nvars )
    {
      for ( index = 0; index < vlistptr->ngrids; index++ )
	if ( vlistptr->gridIDs[index] == vlistptr->vars[varID].gridID )
	  vlistptr->gridIDs[index] = gridID;
    }
  else
    vlistAdd2GridIDs(vlistptr, gridID);

  vlistptr->vars[varID].gridID = gridID;
  reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
}


void vlistDefVarCompType(int vlistID, int varID, int comptype)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if (vlistptr->vars[varID].comptype != comptype)
    {
      vlistptr->vars[varID].comptype = comptype;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


int vlistInqVarCompType(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].comptype);
}


void vlistDefVarCompLevel(int vlistID, int varID, int complevel)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if (vlistptr->vars[varID].complevel != complevel)
    {
      vlistptr->vars[varID].complevel = complevel;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


int vlistInqVarCompLevel(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].complevel);
}


void  vlistDefVarChunkType(int vlistID, int varID, int chunktype)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if (vlistptr->vars[varID].chunktype != chunktype)
    {
      vlistptr->vars[varID].chunktype = chunktype;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


int vlistInqVarChunkType(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].chunktype);
}

static
int vlistEncodeXyz(int (*dimorder)[3])
{
  return (*dimorder)[0]*100 + (*dimorder)[1]*10 + (*dimorder)[2];
}


static
void vlistDecodeXyz(int xyz, int (*outDimorder)[3])
{
  (*outDimorder)[0] = xyz/100, xyz -= (*outDimorder)[0]*100;
  (*outDimorder)[1] = xyz/10, xyz -= (*outDimorder)[1]*10;
  (*outDimorder)[2] = xyz;
}


void  vlistDefVarXYZ(int vlistID, int varID, int xyz)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if ( xyz == 3 ) xyz = 321;

  /* check xyz dimension order */
  {
    int dimorder[3];
    vlistDecodeXyz(xyz, &dimorder);
    int dimx = 0, dimy = 0, dimz = 0;
    for ( int id = 0; id < 3; ++id )
      {
        switch ( dimorder[id] )
          {
            case 1: dimx++; break;
            case 2: dimy++; break;
            case 3: dimz++; break;
            default: dimorder[id] = 0; break;    //Ensure that we assign a valid dimension to this position.
          }
     }
    if ( dimz > 1 || dimy > 1 || dimx > 1 ) xyz = 321; // ZYX
    else
      {
        if ( dimz == 0 ) for ( int id = 0; id < 3; ++id ) if ( dimorder[id] == 0 ) {dimorder[id] = 3; break;}
        if ( dimy == 0 ) for ( int id = 0; id < 3; ++id ) if ( dimorder[id] == 0 ) {dimorder[id] = 2; break;}
        if ( dimx == 0 ) for ( int id = 0; id < 3; ++id ) if ( dimorder[id] == 0 ) {dimorder[id] = 1; break;}
        xyz = vlistEncodeXyz(&dimorder);
      }
  }

  assert(xyz == 123 || xyz == 312 || xyz == 231 || xyz == 321 || xyz == 132 || xyz == 213);

  vlistptr->vars[varID].xyz = xyz;
  reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
}


void vlistInqVarDimorder(int vlistID, int varID, int (*outDimorder)[3])
{
  vlist_t *vlistptr;

  vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  vlistDecodeXyz(vlistptr->vars[varID].xyz, outDimorder);
}


int vlistInqVarXYZ(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return (vlistptr->vars[varID].xyz);
}

/* Ensemble Info Routines */
void vlistDefVarEnsemble(int vlistID, int varID, int ensID, int ensCount, int forecast_type )
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if ( vlistptr->vars[varID].ensdata == NULL )
    vlistptr->vars[varID].ensdata
      = (ensinfo_t *)xmalloc( sizeof( ensinfo_t ) );

  vlistptr->vars[varID].ensdata->ens_index          = ensID;
  vlistptr->vars[varID].ensdata->ens_count          = ensCount;
  vlistptr->vars[varID].ensdata->forecast_init_type = forecast_type;
  reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
}


int vlistInqVarEnsemble( int vlistID, int varID, int *ensID, int *ensCount, int *forecast_type )
{
  int status = 0;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  if ( vlistptr->vars[varID].ensdata )
    {
      *ensID = vlistptr->vars[varID].ensdata->ens_index;
      *ensCount = vlistptr->vars[varID].ensdata->ens_count;
      *forecast_type = vlistptr->vars[varID].ensdata->forecast_init_type;

      status = 1;
    }

  return (status);
}



/* vlistDefVarIntKey: Set an arbitrary keyword/integer value pair for GRIB API */
void vlistDefVarIntKey(int vlistID, int varID, const char *name, int value)
{
#if  defined  (HAVE_LIBGRIB_API)
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  if (vlistptr == NULL)  Error("Internal error!");
  int idx;

  if ( vlistptr->locked != 0 )
    Error("User defined vlist object (vlistID=%d) isn't allowed!\n"
          "Need a CDI internal vlist object from streamInqVlist(streamID).", vlistID);

  for ( idx=0; idx<vlistptr->vars[varID].opt_grib_nentries; idx++)
    if ( (strcmp(name, vlistptr->vars[varID].opt_grib_kvpair[idx].keyword) == 0 ) &&
         (vlistptr->vars[varID].opt_grib_kvpair[idx].data_type == t_int) )  break;

  if ( idx < vlistptr->vars[varID].opt_grib_nentries )
    {
      vlistptr->vars[varID].opt_grib_kvpair[idx].int_val = value;
      vlistptr->vars[varID].opt_grib_kvpair[idx].update  = TRUE;
    }
  else
    {
      resize_opt_grib_entries(&vlistptr->vars[varID], vlistptr->vars[varID].opt_grib_nentries+1);
      vlistptr->vars[varID].opt_grib_nentries += 1;
      idx = vlistptr->vars[varID].opt_grib_nentries -1;
      vlistptr->vars[varID].opt_grib_kvpair[idx].data_type   = t_int;
      vlistptr->vars[varID].opt_grib_kvpair[idx].int_val     = value;
      vlistptr->vars[varID].opt_grib_kvpair[idx].update      = TRUE;
      if ( name )
        vlistptr->vars[varID].opt_grib_kvpair[idx].keyword = strdupx(name);
      else
        Error("Internal error, name undefined!");
    }

  if ( CDI_Debug )
    {
      Message("define additional GRIB2 key \"%s\" (integer): %d", name, value);
      Message("total list of registered, additional GRIB2 keys (total: %d):",
              vlistptr->vars[varID].opt_grib_nentries);
      for ( idx=0; idx<vlistptr->vars[varID].opt_grib_nentries; idx++)
        if (vlistptr->vars[varID].opt_grib_kvpair[idx].data_type == t_int)
          Message("%s -> integer %d",
                  vlistptr->vars[varID].opt_grib_kvpair[idx].keyword,
                  vlistptr->vars[varID].opt_grib_kvpair[idx].int_val);
        else if (vlistptr->vars[varID].opt_grib_kvpair[idx].data_type == t_double)
          Message("%s -> double %d",
                  vlistptr->vars[varID].opt_grib_kvpair[idx].keyword,
                  vlistptr->vars[varID].opt_grib_kvpair[idx].dbl_val);
        else
          Message("%s -> unknown", vlistptr->vars[varID].opt_grib_kvpair[idx].keyword);
    }

  reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
#else
  (void)vlistID;
  (void)varID;
  (void)name;
  (void)value;
#endif
}

/* vlistDefVarDblKey: Set an arbitrary keyword/double value pair for GRIB API */
void vlistDefVarDblKey(int vlistID, int varID, const char *name, double value)
{
#if  defined  (HAVE_LIBGRIB_API)
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  if (vlistptr == NULL)  Error("Internal error!");
  int idx;

  if ( vlistptr->locked != 0 )
    Error("User defined vlist object (vlistID=%d) isn't allowed!\n"
          "Need a CDI internal vlist object from streamInqVlist(streamID).", vlistID);

  for ( idx=0; idx<vlistptr->vars[varID].opt_grib_nentries; idx++)
    if ( (strcmp(name, vlistptr->vars[varID].opt_grib_kvpair[idx].keyword) == 0 ) &&
         (vlistptr->vars[varID].opt_grib_kvpair[idx].data_type == t_double) )  break;

  if ( idx < vlistptr->vars[varID].opt_grib_nentries )
    {
      vlistptr->vars[varID].opt_grib_kvpair[idx].dbl_val = value;
      vlistptr->vars[varID].opt_grib_kvpair[idx].update  = TRUE;
    }
  else
    {
      resize_opt_grib_entries(&vlistptr->vars[varID], vlistptr->vars[varID].opt_grib_nentries+1);
      vlistptr->vars[varID].opt_grib_nentries += 1;
      idx = vlistptr->vars[varID].opt_grib_nentries - 1;
      vlistptr->vars[varID].opt_grib_kvpair[idx].data_type = t_double;
      vlistptr->vars[varID].opt_grib_kvpair[idx].dbl_val   = value;
      vlistptr->vars[varID].opt_grib_kvpair[idx].update    = TRUE;
      if ( name )
        vlistptr->vars[varID].opt_grib_kvpair[idx].keyword = strdupx(name);
      else
        Error("Internal error, name undefined!");
    }

  if ( CDI_Debug )
    {
      Message("define additional GRIB2 key \"%s\" (double): %d", name, value);
      Message("total list of registered, additional GRIB2 keys (total: %d):",
              vlistptr->vars[varID].opt_grib_nentries);
      for ( idx=0; idx<vlistptr->vars[varID].opt_grib_nentries; idx++)
        if (vlistptr->vars[varID].opt_grib_kvpair[idx].data_type == t_int)
          Message("%s -> integer %d",
                  vlistptr->vars[varID].opt_grib_kvpair[idx].keyword,
                  vlistptr->vars[varID].opt_grib_kvpair[idx].int_val);
        else if (vlistptr->vars[varID].opt_grib_kvpair[idx].data_type == t_double)
          Message("%s -> double %d",
                  vlistptr->vars[varID].opt_grib_kvpair[idx].keyword,
                  vlistptr->vars[varID].opt_grib_kvpair[idx].dbl_val);
        else
          Message("%s -> unknown", vlistptr->vars[varID].opt_grib_kvpair[idx].keyword);
    }

  reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
#else
  (void)vlistID;
  (void)varID;
  (void)name;
  (void)value;
#endif
}


/* cdiClearAdditionalKeys: Clears the list of additional GRIB keys. */
void cdiClearAdditionalKeys()
{
#if  defined  (HAVE_LIBGRIB_API)
  for (int i=0; i<cdiNAdditionalGRIBKeys; i++)  free(cdiAdditionalGRIBKeys[i]);
  cdiNAdditionalGRIBKeys = 0;
#endif
}

/* cdiDefAdditionalKey: Register an additional GRIB key which is read when file is opened. */
void cdiDefAdditionalKey(const char *name)
{
#if  defined  (HAVE_LIBGRIB_API)
  int idx = cdiNAdditionalGRIBKeys;
  cdiNAdditionalGRIBKeys++;
  if ( idx >= MAX_OPT_GRIB_ENTRIES ) Error("Too many additional keywords!");
  if ( name )
    cdiAdditionalGRIBKeys[idx] = strdupx(name);
  else
    Error("Internal error!");
#else
  (void)name;
#endif
}

/* vlistHasVarKey: returns 1 if meta-data key was read, 0 otherwise. */
int vlistHasVarKey(int vlistID, int varID, const char* name)
{
#if  defined  (HAVE_LIBGRIB_API)
  /* check if the GRIB key was previously read and is stored */
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  for (int i=0; i<vlistptr->vars[varID].opt_grib_nentries; i++)
    {
      if ( strcmp(name, vlistptr->vars[varID].opt_grib_kvpair[i].keyword) == 0 )
	return 1;
    }
#else
  (void)vlistID;
  (void)varID;
  (void)name;
#endif
  return 0;
}

/* vlistInqVarDblKey: raw access to GRIB meta-data */
double vlistInqVarDblKey(int vlistID, int varID, const char* name)
{
  double value = 0;
#if  defined  (HAVE_LIBGRIB_API)
  /* check if the GRIB key was previously read and is stored in "opt_grib_dbl_val" */
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  for (int i=0; i<vlistptr->vars[varID].opt_grib_nentries; i++)
    {
      int isub = subtypeInqActiveIndex(vlistptr->vars[varID].subtypeID);
      if ( (strcmp(name, vlistptr->vars[varID].opt_grib_kvpair[i].keyword) == 0 ) &&
           (vlistptr->vars[varID].opt_grib_kvpair[i].data_type == t_double)       &&
           (vlistptr->vars[varID].opt_grib_kvpair[i].subtype_index == isub) )
        return vlistptr->vars[varID].opt_grib_kvpair[i].dbl_val;
    }
#else
  (void)vlistID;
  (void)varID;
  (void)name;
#endif
  return value;
}


/* vlistInqVarIntKey: raw access to GRIB meta-data */
int vlistInqVarIntKey(int vlistID, int varID, const char* name)
{
  long int value = 0;
#if  defined  (HAVE_LIBGRIB_API)
  /* check if the GRIB key was previously read and is stored in "opt_grib_int_val" */
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  for (int i=0; i<vlistptr->vars[varID].opt_grib_nentries; i++)
    {
      int isub = subtypeInqActiveIndex(vlistptr->vars[varID].subtypeID);
      if ( (strcmp(name, vlistptr->vars[varID].opt_grib_kvpair[i].keyword) == 0 ) &&
           (vlistptr->vars[varID].opt_grib_kvpair[i].data_type == t_int)          &&
           (vlistptr->vars[varID].opt_grib_kvpair[i].subtype_index == isub) )
        return vlistptr->vars[varID].opt_grib_kvpair[i].int_val;
    }

#else
  (void)vlistID;
  (void)varID;
  (void)name;
#endif
  return (int) value;
}


void vlistDefVarIOrank(int vlistID, int varID, int iorank)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID );

  vlistCheckVarID ( __func__, vlistID, varID );

  if (vlistptr->vars[varID].iorank != iorank)
    {
      vlistptr->vars[varID].iorank = iorank;
      reshSetStatus(vlistID, &vlistOps, RESH_DESYNC_IN_USE);
    }
}


int vlistInqVarIOrank(int vlistID, int varID)
{
  vlist_t *vlistptr = vlist_to_pointer(vlistID);

  vlistCheckVarID(__func__, vlistID, varID);

  return vlistptr->vars[varID].iorank;
}


int vlistVarCompare(vlist_t *a, int varIDA, vlist_t *b, int varIDB)
{
  xassert(a && b
          && varIDA >= 0 && varIDA < a->nvars
          && varIDB >= 0 && varIDB < b->nvars);
  var_t *pva = a->vars + varIDA, *pvb = b->vars + varIDB;
#define FCMP(f) ((pva->f) != (pvb->f))
#define FCMPFLT(f) (IS_NOT_EQUAL((pva->f), (pvb->f)))
#define FCMPSTR(fs) ((pva->fs) != (pvb->fs) && strcmp((pva->fs), (pvb->fs)))
#define FCMP2(f) (namespaceResHDecode(pva->f).idx       \
                  != namespaceResHDecode(pvb->f).idx)
  int diff = FCMP(fvarID) | FCMP(mvarID) | FCMP(flag) | FCMP(param)
    | FCMP(datatype) | FCMP(tsteptype) | FCMP(timave) | FCMP(timaccu)
    | FCMP(chunktype) | FCMP(xyz) | FCMP2(gridID) | FCMP2(zaxisID)
    | FCMP2(instID) | FCMP2(modelID) | FCMP2(tableID) | FCMP(missvalused)
    | FCMPFLT(missval) | FCMPFLT(addoffset) | FCMPFLT(scalefactor) | FCMPSTR(name)
    | FCMPSTR(longname) | FCMPSTR(stdname) | FCMPSTR(units) | FCMPSTR(extra)
    | FCMP(comptype) | FCMP(complevel) | FCMP(lvalidrange)
    | FCMPFLT(validrange[0]) | FCMPFLT(validrange[1]);
#undef FCMP
#undef FCMPFLT
#undef FCMP2
  if ((diff |= ((pva->levinfo == NULL) ^ (pvb->levinfo == NULL))))
    return 1;
  if (pva->levinfo)
    {
      int zaxisID = pva->zaxisID;
      size_t nlevs = (size_t)zaxisInqSize(zaxisID);
      diff |= (memcmp(pva->levinfo, pvb->levinfo, sizeof (levinfo_t) * nlevs)
               != 0);
      if (diff)
        return 1;
    }
  size_t natts = a->vars[varIDA].atts.nelems;
  if (natts != b->vars[varIDB].atts.nelems)
    return 1;
  for (size_t attID = 0; attID < natts; ++attID)
    diff |= vlist_att_compare(a, varIDA, b, varIDB, (int)attID);
  if ((diff |= ((pva->ensdata == NULL) ^ (pvb->ensdata == NULL))))
    return 1;
  if (pva->ensdata)
    diff = (memcmp(pva->ensdata, pvb->ensdata, sizeof (*(pva->ensdata)))) != 0;
  return diff;
}



enum {
  vlistvar_nints = 20,
  vlistvar_ndbls = 3,
};

int vlistVarGetPackSize(vlist_t *p, int varID, void *context)
{
  var_t *var = p->vars + varID;
  int varsize = serializeGetSize(vlistvar_nints, DATATYPE_INT, context)
    + serializeGetSize(vlistvar_ndbls, DATATYPE_FLT64, context);
  if (var->name)
    varsize += serializeGetSize((int)strlen(var->name), DATATYPE_TXT, context);
  if (var->longname)
    varsize += serializeGetSize((int)strlen(var->longname), DATATYPE_TXT, context);
  if (var->stdname)
    varsize += serializeGetSize((int)strlen(var->stdname), DATATYPE_TXT, context);
  if (var->units)
    varsize += serializeGetSize((int)strlen(var->units), DATATYPE_TXT, context);
  varsize += serializeGetSize(4 * zaxisInqSize(var->zaxisID),
                              DATATYPE_INT, context);
  varsize += vlistAttsGetSize(p, varID, context);
  return varsize;
}

void vlistVarPack(vlist_t *p, int varID, char * buf, int size, int *position,
                  void *context)
{
  double dtempbuf[vlistvar_ndbls];
  var_t *var = p->vars + varID;
  int tempbuf[vlistvar_nints], namesz, longnamesz, stdnamesz, unitssz;

  tempbuf[0] = var->flag;
  tempbuf[1] = var->gridID;
  tempbuf[2] = var->zaxisID;
  tempbuf[3] = var->tsteptype;
  tempbuf[4] = namesz = var->name?(int)strlen(var->name):0;
  tempbuf[5] = longnamesz = var->longname?(int)strlen(var->longname):0;
  tempbuf[6] = stdnamesz = var->stdname?(int)strlen(var->stdname):0;
  tempbuf[7] = unitssz = var->units?(int)strlen(var->units):0;
  tempbuf[8] = var->datatype;
  tempbuf[9] = var->param;
  tempbuf[10] = var->instID;
  tempbuf[11] = var->modelID;
  tempbuf[12] = var->tableID;
  tempbuf[13] = var->timave;
  tempbuf[14] = var->timaccu;
  tempbuf[15] = var->missvalused;
  tempbuf[16] = var->comptype;
  tempbuf[17] = var->complevel;
  int nlevs = var->levinfo ? zaxisInqSize(var->zaxisID) : 0;
  tempbuf[18] = nlevs;
  tempbuf[19] = var->iorank;
  dtempbuf[0] = var->missval;
  dtempbuf[1] = var->scalefactor;
  dtempbuf[2] = var->addoffset;
  serializePack(tempbuf, vlistvar_nints, DATATYPE_INT,
                buf, size, position, context);
  serializePack(dtempbuf, vlistvar_ndbls, DATATYPE_FLT64,
                buf, size, position, context);
  if (namesz)
    serializePack(var->name, namesz, DATATYPE_TXT, buf, size, position, context);
  if (longnamesz)
    serializePack(var->longname, longnamesz, DATATYPE_TXT,
                  buf, size, position, context);
  if (stdnamesz)
    serializePack(var->stdname, stdnamesz, DATATYPE_TXT,
                  buf, size, position, context);
  if (unitssz)
    serializePack(var->units, unitssz, DATATYPE_TXT,
                  buf, size, position, context);
  if (nlevs)
    {
      int levbuf[nlevs][4];
      for (int levID = 0; levID < nlevs; ++levID)
        {
          levbuf[levID][0] = var->levinfo[levID].flag;
          levbuf[levID][1] = var->levinfo[levID].index;
          levbuf[levID][2] = var->levinfo[levID].mlevelID;
          levbuf[levID][3] = var->levinfo[levID].flevelID;
        }
      serializePack(levbuf, nlevs * 4, DATATYPE_INT,
                    buf, size, position, context);
    }
  vlistAttsPack(p, varID, buf, size, position, context);
}

static inline int
imax(int a, int b)
{
  return a>=b?a:b;
}


void vlistVarUnpack(int vlistID, char * buf, int size, int *position,
		    int originNamespace, void *context)
{
  double dtempbuf[vlistvar_ndbls];
  int tempbuf[vlistvar_nints];
  char *varname = NULL;
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  serializeUnpack(buf, size, position,
                  tempbuf, vlistvar_nints, DATATYPE_INT, context);
  serializeUnpack(buf, size, position,
                  dtempbuf, vlistvar_ndbls, DATATYPE_FLT64, context);

  /* ------------------------------------------- */
  /* NOTE: Tile sets  currently not supported!!! */
  /* ------------------------------------------- */

  int newvar = vlistDefVar ( vlistID,
			 namespaceAdaptKey ( tempbuf[1], originNamespace ),
			 namespaceAdaptKey ( tempbuf[2], originNamespace ),
                             tempbuf[3]);
  if (tempbuf[4] || tempbuf[5] || tempbuf[6] || tempbuf[7])
    varname = (char *)xmalloc((size_t)imax(imax(imax(tempbuf[4],tempbuf[5]),
                                                tempbuf[6]),
                                           tempbuf[7]) + 1);
  if (tempbuf[4])
  {
    serializeUnpack(buf, size, position,
                    varname, tempbuf[4], DATATYPE_TXT, context);
    varname[tempbuf[4]] = '\0';
    vlistDefVarName(vlistID, newvar, varname);
  }
  if (tempbuf[5])
  {
    serializeUnpack(buf, size, position,
                    varname, tempbuf[5], DATATYPE_TXT, context);
    varname[tempbuf[5]] = '\0';
    vlistDefVarLongname(vlistID, newvar, varname);
  }
  if (tempbuf[6])
  {
    serializeUnpack(buf, size, position,
                    varname, tempbuf[6], DATATYPE_TXT, context);
    varname[tempbuf[6]] = '\0';
    vlistDefVarStdname(vlistID, newvar, varname);
  }
  if (tempbuf[7])
  {
    serializeUnpack(buf, size, position,
                    varname, tempbuf[7], DATATYPE_TXT, context);
    varname[tempbuf[7]] = '\0';
    vlistDefVarUnits(vlistID, newvar, varname);
  }
  free(varname);
  vlistDefVarDatatype(vlistID, newvar, tempbuf[8]);
  vlistDefVarInstitut ( vlistID, newvar,
			namespaceAdaptKey ( tempbuf[10], originNamespace ));
  vlistDefVarModel ( vlistID, newvar,
		     namespaceAdaptKey ( tempbuf[11], originNamespace ));
  vlistDefVarTable(vlistID, newvar, tempbuf[12]);
  /* FIXME: changing the table might change the param code */
  vlistDefVarParam(vlistID, newvar, tempbuf[9]);
  vlistDefVarTimave(vlistID, newvar, tempbuf[13]);
  vlistDefVarTimaccu(vlistID, newvar, tempbuf[14]);
  if (tempbuf[15])
    vlistDefVarMissval(vlistID, newvar, dtempbuf[0]);
  vlistDefVarScalefactor(vlistID, newvar, dtempbuf[1]);
  vlistDefVarAddoffset(vlistID, newvar, dtempbuf[2]);
  vlistDefVarCompType(vlistID, newvar, tempbuf[16]);
  vlistDefVarCompLevel(vlistID, newvar, tempbuf[17]);
  int nlevs = tempbuf[18];
  if (nlevs)
    {
      int levbuf[nlevs][4];
      var_t *var = vlistptr->vars + newvar;
      int i, flagSetLev = 0;
      cdiVlistCreateVarLevInfo(vlistptr, newvar);
      serializeUnpack(buf, size, position,
                      levbuf, nlevs * 4, DATATYPE_INT, context);
      for (i = 0; i < nlevs; ++i)
        {
          vlistDefFlag(vlistID, newvar, i, levbuf[i][0]);
          vlistDefIndex(vlistID, newvar, i, levbuf[i][1]);
          // FIXME: these lack an accessor function
          var->levinfo[i].mlevelID = levbuf[i][2];
          var->levinfo[i].flevelID = levbuf[i][3];
          if (levbuf[i][0] == tempbuf[0])
            flagSetLev = i;
        }
      vlistDefFlag(vlistID, newvar, flagSetLev, levbuf[flagSetLev][0]);
    }
  vlistDefVarIOrank(vlistID, newvar, tempbuf[19]);
  vlistAttsUnpack(vlistID, newvar, buf, size, position, context);
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
