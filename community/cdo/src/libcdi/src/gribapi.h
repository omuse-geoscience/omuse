#ifndef _GRIBAPI_H
#define _GRIBAPI_H

#ifdef HAVE_LIBGRIB_API
#include <grib_api.h>
#ifndef  _ERROR_H
#include "error.h"
#endif
#endif

#ifndef  _CDI_INT_H
#include "cdi_int.h"
#endif

#define  GRIBAPI_MISSVAL  -9.E33

/* GRIB2 Level Types */
#define  GRIB2_LTYPE_SURFACE               1
#define  GRIB2_LTYPE_CLOUD_BASE            2
#define  GRIB2_LTYPE_CLOUD_TOP             3
#define  GRIB2_LTYPE_ISOTHERM0             4
#define  GRIB2_LTYPE_TOA                   8
#define  GRIB2_LTYPE_SEA_BOTTOM            9
#define  GRIB2_LTYPE_ATMOSPHERE           10
#define  GRIB2_LTYPE_ISOBARIC            100
#define  GRIB2_LTYPE_MEANSEA             101
#define  GRIB2_LTYPE_ALTITUDE            102
#define  GRIB2_LTYPE_HEIGHT              103
#define  GRIB2_LTYPE_SIGMA               104
#define  GRIB2_LTYPE_HYBRID              105
#define  GRIB2_LTYPE_LANDDEPTH           106
#define  GRIB2_LTYPE_ISENTROPIC          107
#define  GRIB2_LTYPE_SNOW                114
#define  GRIB2_LTYPE_REFERENCE           150
#define  GRIB2_LTYPE_SEADEPTH            160  /* Depth Below Sea Level                                 */
#define  GRIB2_LTYPE_LAKE_BOTTOM         162  /* Lake or River Bottom                                  */
#define  GRIB2_LTYPE_SEDIMENT_BOTTOM     163  /* Bottom Of Sediment Layer                              */
#define  GRIB2_LTYPE_SEDIMENT_BOTTOM_TA  164  /* Bottom Of Thermally Active Sediment Layer             */
#define  GRIB2_LTYPE_SEDIMENT_BOTTOM_TW  165  /* Bottom Of Sediment Layer Penetrated By Thermal Wave   */
#define  GRIB2_LTYPE_MIX_LAYER           166  /* Mixing Layer                                          */

/* GRIB2 Data representation type (Grid Type) */
#define  GRIB2_GTYPE_LATLON                0  /*  latitude/longitude                                   */
#define  GRIB2_GTYPE_LATLON_ROT            1  /*  rotated latitude/longitude                           */
#define  GRIB2_GTYPE_LATLON_STR            2  /*  stretched latitude/longitude                         */
#define  GRIB2_GTYPE_LATLON_ROTSTR         3  /*  rotated and stretched latitude/longitude             */
#define  GRIB2_GTYPE_GAUSSIAN             40  /*  gaussian grid                                        */
#define  GRIB2_GTYPE_GAUSSIAN_ROT         41  /*  rotated gaussian grid                                */
#define  GRIB2_GTYPE_GAUSSIAN_STR         42  /*  stretched gaussian grid                              */
#define  GRIB2_GTYPE_GAUSSIAN_ROTSTR      43  /*  rotated and stretched gaussian grid                  */
#define  GRIB2_GTYPE_LCC                  30  /*  Lambert conformal                                    */
#define  GRIB2_GTYPE_SPECTRAL             50  /*  spherical harmonics                                  */
#define  GRIB2_GTYPE_GME                 100  /*  hexagonal GME grid                                   */
#define  GRIB2_GTYPE_UNSTRUCTURED        101  /*  General Unstructured Grid                            */

const char *gribapiLibraryVersionString(void);
void gribContainersNew(stream_t * streamptr);
void gribContainersDelete(stream_t * streamptr);

#ifdef HAVE_LIBGRIB_API
static inline void *gribHandleNew(int editionNumber)
{
  void *gh = (void *)grib_handle_new_from_samples(NULL, (editionNumber == 1) ? "GRIB1" : "GRIB2");

  if ( gh == NULL ) Error("grib_handle_new_from_samples failed!");

  return gh;
}

static inline void gribHandleDelete(void *gh)
{
  grib_handle_delete(gh);
}
#else
#define gribHandleNew(editionNumber) (NULL)
#define gribHandleDelete(gh)
#endif

typedef struct {
  int init;
  void *gribHandle;
}
gribContainer_t;

#endif  /* _GRIBAPI_H */
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
