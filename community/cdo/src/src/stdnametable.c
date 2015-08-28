#include <stdio.h>
#include <string.h>
#include <assert.h>
#include "stdnametable.h"


typedef struct
{
  int   varid;
  int   echamcode;
  char *name;
  char *stdname;     /* Standard name */
  char *units;       /* Units         */
}
stdnametable_t;


const stdnametable_t stdnametable[] = {
  /* varid                       code    name                standard name                 units */
  { air_pressure,                   1,  "apres",            "air_pressure",               "Pa" },
  { pressure_thickness,             2,  "dpress",           "pressure_thickness",         "Pa" },
  { surface_geopotential,         129,  "geosp",            "surface_geopotential",       "m2 s-2" },
  { geopotential,                 129,  "z",                "geopotential",               "m2 s-2" },
  { air_temperature,              130,  "ta",               "air_temperature",            "K" },
  { specific_humidity,            133,  "hus",              "specific_humidity",          "1" },
  { surface_air_pressure,         134,  "aps",              "surface_air_pressure",       "Pa" },
  { air_pressure_at_sea_level,    151,  "psl",              "air_pressure_at_sea_level",  "Pa" },
  { geopotential_height,          156,  "zh",               "geopotential_height",        "m" },
};


static int stdnametable_idx(int varid)
{
  int idx;
  int num_entries = (int) (sizeof(stdnametable)/sizeof(stdnametable_t));

  for ( idx = 0; idx < num_entries; ++idx )
    if ( stdnametable[idx].varid == varid ) break;

  assert( idx < num_entries );

  return (idx);
}


int var_echamcode(int varid)
{
  return (stdnametable[stdnametable_idx(varid)].echamcode);
}

const char* var_name(int varid)
{
  return (stdnametable[stdnametable_idx(varid)].name);
}

const char* var_stdname(int varid)
{
  return (stdnametable[stdnametable_idx(varid)].stdname);
}

const char* var_units(int varid)
{
  return (stdnametable[stdnametable_idx(varid)].units);
}

int echamcode_from_stdname(const char* stdname)
{
  int code = -1;

  if      ( strcmp(stdname, var_stdname(surface_geopotential))      == 0 ) code = 129;
  else if ( strcmp(stdname, var_stdname(geopotential))              == 0 ) code = 129;
  else if ( strcmp(stdname, var_stdname(air_temperature))           == 0 ) code = 130;
  else if ( strcmp(stdname, var_stdname(specific_humidity))         == 0 ) code = 133;
  else if ( strcmp(stdname, var_stdname(surface_air_pressure))      == 0 ) code = 134;
  else if ( strcmp(stdname, var_stdname(air_pressure_at_sea_level)) == 0 ) code = 151;
  else if ( strcmp(stdname, var_stdname(geopotential_height))       == 0 ) code = 156;

  return (code);
}


void echam_gribcodes(gribcode_t *gribcodes)
{
  gribcodes->geopot  =  129;
  gribcodes->temp    =  130;
  gribcodes->hum     =  133;
  gribcodes->ps      =  134;
  gribcodes->lsp     =  152;
  gribcodes->gheight =  156;
}


void wmo_gribcodes(gribcode_t *gribcodes)
{
  gribcodes->geopot  =   6;
  gribcodes->temp    =  11;
  gribcodes->hum     =   0;
  gribcodes->ps      =   1;
  gribcodes->lsp     =   0;
  gribcodes->gheight =   7;
}
