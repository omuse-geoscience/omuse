#ifndef _GRID_H
#define _GRID_H

#include "cdi.h"

typedef unsigned char mask_t;

typedef struct {
  int     self;
  int     type;                   /* grid type                      */
  int     prec;                   /* grid precision                 */
  int     proj;                   /* grid projection                */
  mask_t *mask;
  mask_t *mask_gme;
  double *xvals;
  double *yvals;
  double *area;
  double *xbounds;
  double *ybounds;
  double  xfirst, yfirst;
  double  xlast, ylast;
  double  xinc, yinc;
  double  lcc_originLon;          /* Lambert Conformal Conic        */
  double  lcc_originLat;
  double  lcc_lonParY;
  double  lcc_lat1;
  double  lcc_lat2;
  double  lcc_xinc;
  double  lcc_yinc;
  int     lcc_projflag;
  int     lcc_scanflag;
  int     lcc_defined;
  double  lcc2_lon_0;             /* Lambert Conformal Conic 2      */
  double  lcc2_lat_0;
  double  lcc2_lat_1;
  double  lcc2_lat_2;
  double  lcc2_a;
  int     lcc2_defined;
  double  laea_lon_0;             /* Lambert Azimuthal Equal Area   */
  double  laea_lat_0;
  double  laea_a;
  int     laea_defined;
  double  xpole, ypole, angle;    /* rotated north pole             */
  int     isCyclic;               /* TRUE for global cyclic grids   */
  int     isRotated;              /* TRUE for rotated grids         */
  int     xdef;                   /* 0: undefined 1:xvals 2:x0+xinc */
  int     ydef;                   /* 0: undefined 1:yvals 2:y0+yinc */
  int     nd, ni, ni2, ni3;       /* parameter for GRID_GME         */
  int     number, position;       /* parameter for GRID_REFERENCE   */
  char   *reference;
  unsigned char uuid[CDI_UUID_SIZE]; /* uuid for grid reference        */
  int     trunc;                  /* parameter for GRID_SPECTEAL    */
  int     nvertex;
  int    *rowlon;
  int     nrowlon;
  int     size;
  int     xsize;                  /* number of values along X */
  int     ysize;                  /* number of values along Y */
  int     np;                     /* number of parallels between a pole and the equator */
  int     locked;
  int     lcomplex;
  int     hasdims;
  char    xname[CDI_MAX_NAME];
  char    yname[CDI_MAX_NAME];
  char    xlongname[CDI_MAX_NAME];
  char    ylongname[CDI_MAX_NAME];
  char    xstdname[CDI_MAX_NAME];
  char    ystdname[CDI_MAX_NAME];
  char    xunits[CDI_MAX_NAME];
  char    yunits[CDI_MAX_NAME];
  char   *name;
}
grid_t;


void grid_init(grid_t *gridptr);
void grid_free(grid_t *gridptr);

unsigned cdiGridCount(void);

const double *gridInqXvalsPtr(int gridID);
const double *gridInqYvalsPtr(int gridID);

const double *gridInqXboundsPtr(int gridID);
const double *gridInqYboundsPtr(int gridID);
const double *gridInqAreaPtr(int gridID);

int gridCompare(int gridID, const grid_t *grid);
int gridGenerate(const grid_t *grid);

void cdiGridGetIndexList(unsigned, int * );

void
gridUnpack(char * unpackBuffer, int unpackBufferSize,
           int * unpackBufferPos, int originNamespace, void *context,
           int force_id);

int  varDefGrid(int vlistID, const grid_t *grid, int mode);

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
