#ifndef _STDNAMETABLE_H
#define _STDNAMETABLE_H

enum stdnameid {air_pressure,
		pressure_thickness,
                surface_geopotential,
		geopotential,
		air_temperature,
                specific_humidity,
		surface_air_pressure,
		air_pressure_at_sea_level,
		geopotential_height};

int         var_echamcode(int varid);
const char* var_name(int varid);
const char* var_stdname(int varid);
const char* var_units(int varid);

int echamcode_from_stdname(const char* stdname);

typedef struct
{
  int geopot;
  int temp;
  int hum;
  int ps;
  int lsp;
  int gheight;
} gribcode_t;

void echam_gribcodes(gribcode_t *gribcodes);
void wmo_gribcodes(gribcode_t *gribcodes);

#endif
