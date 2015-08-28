#ifndef _AFTERBURNER_H
#define _AFTERBURNER_H

/* =============================================== */
/* These include files should be standard on all   */
/* UNIX systems.                                   */
/* =============================================== */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <errno.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <pwd.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <limits.h>

#ifndef _ERROR_H
#  include "error.h"
#endif
#ifndef _DMEMORY_H
#  include "dmemory.h"
#endif

#ifndef TRUE
#define TRUE  1
#endif
#ifndef FALSE
#define FALSE 0
#endif

#define MaxLevel 1024

#define MaxCodes  277

#define S_ECHAM5  1

struct Date
{
   int yr;
   int mo;
   int dy;
   int hr;
   int mn;
};

struct Control
{
  int    Verbose;
 
  int    Mean;
  int    MeanCount0;
  int    MeanCount;
  int    Multi;
  int    Nfiles;
  int    TermCount;

  int    OutputInterval;
  int    EndOfInterval;

  int    AnalysisData; /* 0 = ECHAM Data, 1 = ECMWF Spectral Analyses */
  int    DayIn;        /* day increment of infiles if Multi = TRUE    */
  int    Debug;
  int    Extrapolate;
  int    Szip;

  int    istreamID;
  int    ostreamID;
  int    ostreamID2;
  int    ivlistID;
  int    ovlistID;
  int    ovlistID2;
  int    taxisID;
  int    taxisID2;

  struct Date NextDate;
  struct Date NewDate;
  struct Date OldDate;
  struct Date StartDate;

  int    nvct;
  double *vct;

  int    *vert_index;
  int    *pnmiss;
  double *Orography;
  double *p_of_height;

  int    Type;
  int    unitsel;
  
  int    Fouriers;
  int    Latitudes;
  int    Longitudes;
  int    HalfLevels;
  int    Gaussian;
  int    Spectral;

  int    Truncation;
  int    Waves;

  int    Dim3FC,    Dim3SP,    Dim3GP;
  int    DimFC,     DimGP,     DimSP;
  int    DimSP_half;

  double *poli;
  double *pold;
  double *pdev;
  double *pol2;
  double *pol3;

  double *dv2uv_f1;
  double *dv2uv_f2;

  int    NumCodesRequest;
  
  int    NumLevel;
  int    NumLevelFound;
  int    NumLevelRequest;
  double LevelRequest[MaxLevel];

  double *rcoslat;
  double *coslat;
  double *DerivationFactor;
  double *Field;
};
  
struct Variable
{
  int     needed0;   /* var needed for process  */
  int     needed;    /* var needed for process  */
  int     selected;  /* var selected for output */
  int     detected;  /* var detected in input   */
  int     comp;      /* compute var if selected and not detected */
  int     sfit;
  int     hlev;
  int     plev;
  int     ivarID;
  int     ovarID;    /* 1st variable ID */
  int     ovarID2;   /* 2nd variable ID used for variance */
  int     tableID;
  int     igridID;
  int     ogridID;
  int     izaxisID;
  int     ozaxisID;
  int     nmiss0;
  int     nmiss;
  double  missval;
  double *spectral;
  double *spectral0;
  double *fourier;
  double *hybrid;
  double *hybrid0;
  double *height;
  double *grid;
  double *grid0;
  double *mean;
  double *variance;
  int    *samp;
};

/* FFT */
void fft_set(double *trigs, long *ifax, long n);
void fc2gp(double *restrict trig, long *restrict ifax, double *restrict fc, double *restrict gp, long nlat, long nlon, long nlev, long nfc);
void gp2fc(double *trig, long *ifax, const double *restrict gp, double *restrict fc, long nlat, long nlon, long nlev, long nfc);

/* Convert Spectral Array to new resolution */
void sp2sp(double *arrayIn, int truncIn, double *arrayOut, int truncOut);
void sp2fc(const double *sa, double *fa, const double *poli, long nlev, long nlat, long nfc, long nt);
void fc2sp(double *fa, double *sa, double *poli, int klev, int nlat, int nfc, int nt);

/* Physc */
void dv2ps(const double * restrict div, double * restrict pot, long nlev, long ntr);
void dv2uv(double *d, double *o, double *u, double *v, double *f, double *g, int nt, int nsp, int nlev);
void scaluv(double *fu, double rclat[], int nlat, int lot);
void uv2dv(double *fu, double *fv, double *sd, double *sv,
           double *pol2, double *pol3, int klev, int nlat, int nt);
void geninx(long ntr, double *f, double *g);

#define    LOW_CLOUD   34
#define    MID_CLOUD   35
#define    HIH_CLOUD   36
#define    LOW_WATER   37  /* not used ?   */
#define    MID_WATER   38  /* not used ?   */
#define    HIH_WATER   39  /* not used ?   */
#define    ALL_WATER   40  /* not used ?   */

#define GEOPOTENTIAL  129
#define  TEMPERATURE  130
#define       U_WIND  131
#define       V_WIND  132
#define     HUMIDITY  133
#define           PS  134
#define        OMEGA  135
#define    VORTICITY  138
#define           TS  139
#define       STREAM  148
#define      VELOPOT  149
#define          SLP  151
#define         LNPS  152
#define   DIVERGENCE  155
#define GEOPOTHEIGHT  156
#define    RHUMIDITY  157

#define   SW_BOT_CLF  189  /* not used ?   */
#define   LW_BOT_CLF  190  /* not used ?   */
#define   SW_TOP_CLF  191  /* not used ?   */
#define   LW_TOP_CLF  192  /* not used ?   */
#define  NET_TOP_CLF  193  /* not computed */

#define    WINDSPEED  259
#define       PRECIP  260
#define      NET_TOP  261
#define      NET_BOT  262
#define     NET_HEAT  263
#define    NET_WATER  264
#define       SW_CLF  265
#define       LW_CLF  266
#define      NET_CLF  267
#define       SW_ATM  268
#define       LW_ATM  269
#define      NET_ATM  270
#define  SURF_RUNOFF  271
#define        DPSDX  273
#define        DPSDY  274
#define  FRESH_WATER  275
#define      PS_PROG  276  /* PS for prognostic timestep */
#define   HALF_PRESS  277
#define   FULL_PRESS  278
#define       THETAH  279
#define       THETAF  280

void after_read_vct(const char *vctfile, double **vct, int *nvct);

void after_gp2sp(struct Control *globs, struct Variable *vars, int ccode);
void after_GP2FC(double *gp, double *fc, long nlat, long nlon, long nlev, long nfc);
void after_FC2GP(double *fc, double *gp, long nlat, long nlon, long nlev, long nfc);
void after_FCrh2FCsh(struct Control *globs, struct Variable *vars);
void after_SPuv2SPdv(struct Control *globs, struct Variable *vars);
void after_FCsh2FCrh(struct Control *globs, struct Variable *vars);

void after_EchamCompGP(struct Control *globs, struct Variable *vars);
void after_processPL(struct Control *globs, struct Variable *vars);
void after_processML(struct Control *globs, struct Variable *vars);

void after_AnalysisAddRecord(struct Control *globs, struct Variable *vars, int code, int gridID, int zaxisID, int levelID, int nmiss);
void after_EchamAddRecord(struct Control *globs, struct Variable *vars, int code, int gridID, int zaxisID, int levelID, int nmiss);

void  after_AnalysisDependencies(struct Variable *vars, int ncodes);
void  after_EchamDependencies(struct Variable *vars, int ncodes, int type, int source);

void after_legini_setup(struct Control *globs, struct Variable *vars);

#endif /*  afterburner.h  */
