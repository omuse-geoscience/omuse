#include <stdio.h>
#include <stdlib.h>

#include "cdi.h"

#if defined(HAVE_CONFIG_H)
#include "config.h"
#if defined(CDO)
#include "cdo_int.h"
#include "pstream_write.h"
#else
#define  OPENMP4  201307
#if defined(_OPENMP) && defined(OPENMP4) && _OPENMP >= OPENMP4
#define  HAVE_OPENMP4  1
#endif
#endif
#endif

#include "error.h"
#include "afterburner.h"
#include "constants.h"
#include "compare.h"
#include "after_vertint.h"

int afterDebug = 0;
int labort_after = TRUE;

#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define MIN(a,b) ((a) < (b) ? (a) : (b))


static
char *FieldName(int code, const char *text)
{
  static char name[256];

  sprintf(name, "[%3d].%s", code, text);

  return (name);
}

/* Free array space */
static
void *FreeMemory(void *ptr)
{
  free(ptr);

  return (NULL);
}

static
void FreeSpectral(struct Variable *vars)
{
  int code;

  for ( code = MaxCodes-1; code >= 0; --code )
    if ( vars[code].spectral )
      vars[code].spectral = FreeMemory(vars[code].spectral);
}

static
void FreeFourier(struct Variable *vars)
{
  int code;

  for ( code = 0; code < MaxCodes; code++ )
    if ( vars[code].fourier )
      vars[code].fourier = FreeMemory(vars[code].fourier);
}

static
void FreeHybrid(struct Variable *vars)
{
  int code;

  for ( code = 0; code < MaxCodes; code++ )
    if ( vars[code].hybrid )
      vars[code].hybrid = FreeMemory(vars[code].hybrid);
}

static
void FreeGrid(struct Variable *vars)
{
  int code;

  for ( code = 0; code < MaxCodes; code++ )
    if ( vars[code].grid )
      vars[code].grid = FreeMemory(vars[code].grid);
}

static
void FreeSamp(struct Variable *vars)
{
  int code;

  for ( code = 0; code < MaxCodes; code++ )
    if ( vars[code].samp )
      vars[code].samp = FreeMemory(vars[code].samp);
}

/* alloc_dp - Allocate space for double array */
static
double *alloc_dp(int words, char *array_name)
{
  double *result = NULL;

  if ( words > 0 )
    {
      result = (double *) malloc(words * sizeof(double));

      if ( result == NULL ) SysError(array_name, "No Memory!");
    }

  return(result);
}

/* after_copy_array - Copy array of type double */
static
void after_copy_array(void *destination, void *source, int words)
{
   memcpy(destination, source, words*sizeof(double));
}

/* after_zero_array -  Set array of type double to zero */
static
void after_zero_array(double *field, int words)
{
   memset((char *)field, 0, words*sizeof(double));
}

static
void IniQuaSum(double *dest, const double *restrict src, int len)
{
  int i;

  for ( i = 0; i < len; i++ )
    dest[i] = src[i] * src[i];
}

static 
void AddQuaSum(double *dest, const double *restrict src, int len)
{
  int i;

  for ( i = 0; i < len; i++ )
    dest[i] += src[i] * src[i];
}

static
void VarQuaSum(double *Variance, const double *restrict Sum, int len, int n)
{
  int i;
  double rn1;

  rn1 = 1.0 / (n-1);

  for ( i = 0; i < len; i++ )
    Variance[i] = (Variance[i] - Sum[i] * Sum[i] * n) * rn1;

  for ( i = 0; i < len; i++ )
    if (Variance[i] > 0.0) Variance[i] = sqrt(Variance[i]);
    else                   Variance[i] = 0.0;
}

static
void AddVector(double * restrict dest, const double *restrict src, long len, int *nmiss, double missval)
{
  long i;

  if ( *nmiss > 0 )
    {
      for ( i = 0; i < len; i++ )
	if ( IS_NOT_EQUAL(src[i], missval) )
	  {
	    if ( IS_NOT_EQUAL(dest[i], missval) )
	      dest[i] += src[i];
	    else
	      dest[i] = src[i];
	  }

      *nmiss = 0;
      for ( i = 0; i < len; i++ )
	if ( IS_EQUAL(dest[i], missval) ) *nmiss = *nmiss + 1;

      if ( *nmiss == 0 ) *nmiss = 1;
    }
  else
    {
#if defined (_OPENMP)
#pragma omp parallel for
#endif
      for ( i = 0; i < len; i++ ) dest[i] += src[i];
    }
}

static
void Add2Vectors(double *dest, const double *restrict srcA, const double *restrict srcB, int len)
{
  int i;

  for ( i = 0; i < len; i++ )
    dest[i] = srcA[i] + srcB[i];
}

static
void Sub2Vectors(double *dest, const double *restrict srcA, const double *restrict srcB, int len)
{
  int i;

  for ( i = 0; i < len; i++ )
    dest[i] = srcA[i] - srcB[i];
}

static
void MultVectorScalar(double *dest, const double *restrict src, double factor, int len, int nmiss, double missval)
{
  int i;

  if ( nmiss > 0 )
    {
      for ( i = 0; i < len; i++ )
	{
	  if ( IS_EQUAL(src[i], missval) )
	    dest[i] = missval;
	  else
	    dest[i] = src[i] * factor;
	}
    }
  else
    {
      for ( i = 0; i < len; i++ ) dest[i] = src[i] * factor;
    }
}

static
void DivVectorIvector(double *dest, const double *restrict src, const int *samp, int len, int *nmiss, double missval)
{
  int i;

  *nmiss = 0;

  for ( i = 0; i < len; i++ )
    {
      if ( IS_EQUAL(src[i], missval) || samp[i] == 0 )
	{
	  dest[i] = missval;
	  *nmiss = *nmiss + 1;
	}
      else
	dest[i] = src[i] / samp[i];
    }
}


void after_read_vct(const char *vctfile, double **vct, int *nvct)
{
  int n;
  char line[1024];
  double va, vb;

  FILE *fp = fopen(vctfile, "r");
  if ( fp == NULL ) SysError( "Open failed on %s", vctfile);

  while ( fgets(line, 1023, fp) ) nvct++;

  *nvct *= 2;
  *vct = (double *) malloc(*nvct*sizeof(double));

  rewind(fp);
  for ( int i = 0; i < *nvct/2; i++ )
    {
      fgets(line, 1023, fp);
      sscanf(line, "%d %lg %lg", &n, &va, &vb);
      *vct[i]         = va;
      *vct[i+*nvct/2] = vb;
    }
  fprintf(stdout, "  Reading VCT for %d hybrid levels from file %s\n", *nvct/2-1, vctfile);

  fclose(fp);
}


void after_gp2sp(struct Control *globs, struct Variable *vars, int ccode)
{
  struct Variable *var = &vars[ccode];
  
  if ( var->spectral == NULL )
    {
      if ( var->hybrid == NULL )
	{
	  fprintf(stderr,"%d.hybrid not found\n", ccode);
	  exit(99);
	}
      
      if ( var->fourier == NULL )
	{
	  long fieldSize = globs->DimFC * var->hlev;
	  var->fourier = alloc_dp(fieldSize, "gp2sp.fourier");
	  after_GP2FC(var->hybrid, var->fourier,
		      globs->Latitudes, globs->Longitudes, var->hlev, globs->Fouriers);
	}

      var->spectral = alloc_dp(globs->Dim3SP, "gp2sp.spectral");
      fc2sp(var->fourier, var->spectral,
            globs->pold, var->hlev, globs->Latitudes, globs->Fouriers, globs->Truncation);    
    }
}


void after_GP2FC(double *gp, double *fc, long nlat, long nlon, long nlev, long nfc)
{
  static long ifax[10];
  static double *trig = NULL;

  if ( ifax[9] != nlon )
    {
      if ( trig ) free (trig);
      trig = (double *) malloc(nlon * sizeof(double));
      fft_set (trig, ifax, nlon);
    }

  gp2fc(trig, ifax, gp, fc, nlat, nlon, nlev, nfc);
}


void after_FC2GP(double *fc, double *gp, long nlat, long nlon, long nlev, long nfc)
{
  static long ifax[10];
  static double *trig = NULL;

  if ( ifax[9] != nlon )
    {
      if ( trig ) free (trig);
      trig = (double *) malloc(nlon * sizeof(double));
      fft_set (trig, ifax, nlon);
    }

  fc2gp(trig, ifax, fc, gp, nlat, nlon, nlev, nfc);
}

/* HUMTEST */

static
void sh2rh(int AnalysisData, double *sphum, double *rhum, double *t, int lev,
	   int dimgpout, double *level, double *fullpresshybrid)
{
   int lp,i;
   int lpi,lfp;
   double  es, qsatr;
   double *fullp;
   double  RALPW, RBETW, RGAMW;
   // double  RALPS, RBETS, RGAMS;
   double  RALP , RBET , RGAM ;

   /* ***************************************************** */
   /* Define constants for calculation in presence of water */
   /* ***************************************************** */
   RGAMW = (C_RCW - C_RCPV) / C_RV;
   RBETW = C_RLVTT / C_RV + RGAMW * C_RTT;
   RALPW = log(C_RESTT) + RBETW / C_RTT + RGAMW * log(C_RTT);

   /* ***************************************************** */
   /* Define constants for calculation in presence of  ice  */
   /* ***************************************************** */
   /*
   RGAMS = (C_RCS - C_RCPV) / C_RV;
   RBETS = C_RLSTT / C_RV + RGAMS * C_RTT;
   RALPS = log(C_RESTT) + RBETS / C_RTT + RGAMS * log(C_RTT);
   */
   if (AnalysisData) fullp = level;
   else              fullp = fullpresshybrid;

   /***************************************************/
   /* Diagnostics of saturation water vapour pressure */
   /* over ice makes no sense, therefore ...          */
   /* Hint of Michael Ponater                08.10.97 */
   /***************************************************/
   RGAM = RGAMW; RBET = RBETW; RALP = RALPW;
   for (lp = 0; lp < lev; lp++) {
      for (i = 0; i < dimgpout; i++) {
         lpi = lp*dimgpout + i;
         lfp = (1 - AnalysisData)*lpi + AnalysisData*lp;
	 /*
	 if (t[lpi] < C_RTT) {
	   RGAM = RGAMS; RBET = RBETS; RALP = RALPS;
	 } else {
	   RGAM = RGAMW; RBET = RBETW; RALP = RALPW;
	 }
	 */
         es = (exp(RALP - RBET / t[lpi] - RGAM * log(t[lpi]))) / fullp[lfp];
         // qsat = es / (1. + C_RETV * (1. - es));
	 qsatr = (1. + C_RETV * (1. - es)) / es;
         rhum[lpi] = sphum[lpi] * 100. * qsatr;
      }
   }
}

static
void rh2sh(double *sphum, double *rhum, double *t, int lev, int dimgpout, double *level)
{
   int lp,i;
   int lpi;
   double  es,qsat;
   double  RALPW, RBETW, RGAMW;
   // double  RALPS, RBETS, RGAMS;
   double  RALP , RBET , RGAM ;

/* ***************************************************** */
/* Define constants for calculation in presence of water */
/* ***************************************************** */
   RGAMW = (C_RCW - C_RCPV) / C_RV;
   RBETW = C_RLVTT / C_RV + RGAMW * C_RTT;
   RALPW = log(C_RESTT) + RBETW / C_RTT + RGAMW * log(C_RTT);

/* ***************************************************** */
/* Define constants for calculation in presence of  ice  */
/* ***************************************************** */
   // RGAMS = (C_RCS - C_RCPV) / C_RV;
   // RBETS = C_RLSTT / C_RV + RGAMS * C_RTT;
   // RALPS = log(C_RESTT) + RBETS / C_RTT + RGAMS * log(C_RTT);

/***************************************************/
/* Diagnostics of saturation water vapour pressure */
/* over ice makes no sense, therefore ...          */
/* Hint of Michael Ponater                08.10.97 */
/***************************************************/

   RGAM = RGAMW; RBET = RBETW; RALP = RALPW;
   for (lp = 0; lp < lev; lp++) {
      for (i = 0; i < dimgpout; i++) {
         lpi = lp*dimgpout + i;
/*       if (t[lpi] < C_RTT) { */
/*          RGAM = RGAMS; RBET = RBETS; RALP = RALPS; */
/*       }  else { */
/*          RGAM = RGAMW; RBET = RBETW; RALP = RALPW; */
/*       } */
         es = (exp(RALP - RBET / t[lpi] - RGAM * log(t[lpi]))) / level[lp];
         qsat = es / (1. + C_RETV * (1. - es));
         sphum[lpi] = rhum[lpi] * qsat /  100.;
      }
   }
}


void after_FCrh2FCsh(struct Control *globs, struct Variable *vars)
{
   long fieldSize = globs->DimGP * globs->NumLevelRequest;

   if ( vars[RHUMIDITY].grid == NULL )
     vars[RHUMIDITY].grid = alloc_dp(fieldSize, "vars[RHUMIDITY].grid");
   if ( vars[TEMPERATURE].grid == NULL )
     vars[TEMPERATURE].grid = alloc_dp(fieldSize, "vars[TEMPERATURE].grid");
   if ( vars[HUMIDITY].grid == NULL )
     vars[HUMIDITY].grid = alloc_dp(fieldSize, "vars[HUMIDITY].grid");

   after_FC2GP(vars[RHUMIDITY].fourier, vars[RHUMIDITY].grid,
	       globs->Latitudes, globs->Longitudes, vars[RHUMIDITY].plev, globs->Fouriers);
   after_FC2GP(vars[TEMPERATURE].fourier, vars[TEMPERATURE].grid,
	       globs->Latitudes, globs->Longitudes, vars[TEMPERATURE].plev, globs->Fouriers);

   rh2sh(vars[HUMIDITY].grid, vars[RHUMIDITY].grid, vars[TEMPERATURE].grid, globs->NumLevelRequest,
	 globs->DimGP, globs->LevelRequest);

   after_GP2FC(vars[HUMIDITY].grid, vars[HUMIDITY].fourier,
	       globs->Latitudes, globs->Longitudes, vars[HUMIDITY].plev, globs->Fouriers);

   vars[HUMIDITY].grid    = FreeMemory(vars[HUMIDITY].grid);
   vars[RHUMIDITY].grid   = FreeMemory(vars[RHUMIDITY].grid);
   vars[TEMPERATURE].grid = FreeMemory(vars[TEMPERATURE].grid);
}


void after_SPuv2SPdv(struct Control *globs, struct Variable *vars)
{
   int i;
   long fieldSize;
   double *Div, *DivOut, *Vor, *VorOut;

   Div = DivOut = vars[DIVERGENCE].spectral;
   Vor = VorOut =  vars[VORTICITY].spectral;
   fieldSize = globs->DimFC * globs->NumLevelRequest;

   if ( vars[U_WIND].fourier == NULL )
     vars[U_WIND].fourier = alloc_dp(fieldSize, "vars[U_WIND].fourier");
   if ( vars[V_WIND].fourier == NULL )
     vars[V_WIND].fourier = alloc_dp(fieldSize, "vars[V_WIND].fourier");

   sp2fc(vars[U_WIND].spectral, vars[U_WIND].fourier, globs->poli,
         globs->NumLevelRequest, globs->Latitudes, globs->Fouriers, globs->Truncation);
   sp2fc(vars[V_WIND].spectral, vars[V_WIND].fourier, globs->poli,
         globs->NumLevelRequest, globs->Latitudes, globs->Fouriers, globs->Truncation);
   uv2dv(vars[U_WIND].fourier, vars[V_WIND].fourier, Div, Vor,
         globs->pol2, globs->pol3, globs->NumLevelRequest, globs->Latitudes, globs->Truncation);

   vars[U_WIND].fourier = FreeMemory(vars[U_WIND].fourier);
   vars[V_WIND].fourier = FreeMemory(vars[V_WIND].fourier);

   for (i = 0; i < globs->NumLevelRequest; ++i) {
      sp2sp(Div, globs->Truncation, DivOut, globs->Truncation);
      sp2sp(Vor, globs->Truncation, VorOut, globs->Truncation);
      Div    += globs->DimSP;
      Vor    += globs->DimSP;
      DivOut += globs->DimSP;
      VorOut += globs->DimSP;
   }
}


void after_FCsh2FCrh(struct Control *globs, struct Variable *vars)
{
   long fieldSize;

   fieldSize = globs->DimGP * globs->NumLevelRequest;

   if ( vars[RHUMIDITY].grid == NULL )
     vars[RHUMIDITY].grid = alloc_dp(fieldSize, "vars[RHUMIDITY].grid");
   if ( vars[TEMPERATURE].grid == NULL )
     vars[TEMPERATURE].grid = alloc_dp(fieldSize, "vars[TEMPERATURE].grid");
   if ( vars[HUMIDITY].grid == NULL )
     vars[HUMIDITY].grid = alloc_dp(fieldSize, "vars[HUMIDITY].grid");

   after_FC2GP(vars[HUMIDITY].fourier, vars[HUMIDITY].grid,
	       globs->Latitudes, globs->Longitudes, vars[HUMIDITY].plev, globs->Fouriers);
   after_FC2GP(vars[TEMPERATURE].fourier, vars[TEMPERATURE].grid,
	       globs->Latitudes, globs->Longitudes, vars[TEMPERATURE].plev, globs->Fouriers);

   sh2rh(globs->AnalysisData, vars[HUMIDITY].grid, vars[RHUMIDITY].grid, vars[TEMPERATURE].grid, globs->NumLevelRequest,
	 globs->DimGP, globs->LevelRequest, vars[FULL_PRESS].hybrid);

   after_GP2FC(vars[RHUMIDITY].grid, vars[RHUMIDITY].fourier,
	       globs->Latitudes, globs->Longitudes, vars[RHUMIDITY].plev, globs->Fouriers);

   vars[HUMIDITY].grid    = FreeMemory(vars[HUMIDITY].grid);
   vars[RHUMIDITY].grid   = FreeMemory(vars[RHUMIDITY].grid);
   vars[TEMPERATURE].grid = FreeMemory(vars[TEMPERATURE].grid);
}
/* ENDE HUMTEST */

static
void CheckAnalyses(struct Variable *vars)
{
  for ( int code = 0; code < 272; code++ )
    if ( vars[code].needed     && 
	 code != DIVERGENCE    &&
	 code != VORTICITY     &&
	 code != STREAM        &&
	 code != U_WIND        &&
	 code != HUMIDITY      &&
	 code != VELOPOT       &&
	 code != V_WIND        &&
	 code != RHUMIDITY     &&
	 code != GEOPOTHEIGHT  && 
	 code != PS            &&
         vars[code].spectral == NULL && vars[code].grid == NULL )
      {
	if ( labort_after )
	  Error("Code  %3d not found", code);
	else
	  Warning("Code  %3d not found", code);
      }
}

/*  Process Pressure Level data */
void after_processPL(struct Control *globs, struct Variable *vars)
{
  int code, l;
  long fieldSize;
  int lindex, nlevel;
  int offset;

  globs->MeanCount++;
  globs->TermCount++;

  if ( globs->MeanCount == 1 )
    {
      if ( globs->Debug ) fprintf(stderr, "CheckAnalyses: %d %d\n", globs->TermCount, globs->MeanCount);
      CheckAnalyses(vars);
      globs->StartDate = globs->OldDate;
    }
  if ( globs->TermCount  > 120 ) globs->Debug = 0;

  /* ============================== */
  /* Computations in spectral space */
  /* ============================== */

  if ( vars[TEMPERATURE].needed )
    {
      vars[TEMPERATURE].hlev = 2;
      vars[TEMPERATURE].plev = globs->NumLevelRequest;
      vars[TEMPERATURE].sfit = TRUE;
    }

  if ( vars[GEOPOTHEIGHT].comp )
    {
      vars[GEOPOTHEIGHT].hlev = 2;
      vars[GEOPOTHEIGHT].plev = globs->NumLevelRequest;
      vars[GEOPOTHEIGHT].sfit = TRUE;
    }

  if ( vars[GEOPOTHEIGHT].comp && vars[GEOPOTENTIAL].detected )
    {
      if ( vars[GEOPOTHEIGHT].spectral == NULL )
	vars[GEOPOTHEIGHT].spectral = alloc_dp(globs->DimSP*globs->NumLevelRequest, "GEOPOTHEIGHT.spectral");
      MultVectorScalar(vars[GEOPOTHEIGHT].spectral, vars[GEOPOTENTIAL].spectral,
		       C_RG, globs->DimSP*globs->NumLevelRequest, 0, 0);
      vars[GEOPOTENTIAL].needed = vars[GEOPOTENTIAL].selected;
    }

  if ( globs->Type == 50 && vars[HUMIDITY].needed && vars[HUMIDITY].spectral == NULL )
    {
      vars[HUMIDITY].plev = globs->NumLevelRequest;
      vars[HUMIDITY].sfit = TRUE;
      vars[HUMIDITY].spectral = alloc_dp(globs->DimSP*globs->NumLevelRequest, "vars[HUMIDITY].spectral");
      /*
	SPrh2SPsh();
      */
      vars[RHUMIDITY].needed = vars[RHUMIDITY].selected;
      vars[TEMPERATURE].needed = vars[TEMPERATURE].selected;
    }

  if ( vars[U_WIND].spectral && vars[V_WIND].spectral &&
       (vars[DIVERGENCE].comp || vars[VORTICITY].comp) )
    {
      vars[DIVERGENCE].hlev = vars[VORTICITY].hlev =    2;
      vars[DIVERGENCE].plev = vars[VORTICITY].plev = globs->NumLevelRequest;
      vars[DIVERGENCE].sfit = vars[VORTICITY].sfit = TRUE;
      if ( vars[DIVERGENCE].spectral == NULL )
	vars[DIVERGENCE].spectral = alloc_dp(globs->DimSP*globs->NumLevelRequest, "vars[DIVERGENCE].spectral");
      if ( vars[VORTICITY].spectral == NULL )
	vars[VORTICITY].spectral = alloc_dp(globs->DimSP*globs->NumLevelRequest, "vars[VORTICITY].spectral");
      after_SPuv2SPdv(globs, vars);
    }

  if ( vars[U_WIND].comp || vars[V_WIND].comp )
    {
      vars[U_WIND].hlev = vars[V_WIND].hlev =    2;
      vars[U_WIND].plev = vars[V_WIND].plev = globs->NumLevelRequest;
      vars[U_WIND].sfit = vars[V_WIND].sfit = TRUE;
      if ( vars[U_WIND].spectral == NULL )
	vars[U_WIND].spectral = alloc_dp(globs->DimSP*globs->NumLevelRequest, "vars[U_WIND].spectral");
      if ( vars[V_WIND].spectral == NULL )
	vars[V_WIND].spectral = alloc_dp(globs->DimSP*globs->NumLevelRequest, "vars[V_WIND].spectral");
      dv2uv(vars[DIVERGENCE].spectral, vars[VORTICITY].spectral,
	    vars[U_WIND].spectral, vars[V_WIND].spectral,
	    globs->dv2uv_f1, globs->dv2uv_f2,
	    globs->Truncation, globs->DimSP, globs->NumLevelRequest);
    }

  if ( vars[VELOPOT].comp )
    {
      vars[VELOPOT].hlev =    2;
      vars[VELOPOT].plev = globs->NumLevelRequest;
      vars[VELOPOT].sfit = TRUE;
      if ( vars[VELOPOT].spectral == NULL )
	vars[VELOPOT].spectral = alloc_dp(globs->DimSP*globs->NumLevelRequest, "vars[VELOPOT].spectral");
      dv2ps(vars[DIVERGENCE].spectral, vars[VELOPOT].spectral, globs->NumLevelRequest, globs->Truncation);
    }

  if ( vars[STREAM].comp )
    {
      vars[STREAM].hlev =    2;
      vars[STREAM].plev = globs->NumLevelRequest;
      vars[STREAM].sfit = TRUE;
      if ( vars[STREAM].spectral == NULL )
	vars[STREAM].spectral = alloc_dp(globs->DimSP*globs->NumLevelRequest, "vars[STREAM].spectral");
      dv2ps(vars[VORTICITY].spectral, vars[STREAM].spectral, globs->NumLevelRequest, globs->Truncation);
    }

  /* --------------------------- */
  /*  Output of spectral fields  */
  /* --------------------------- */

  if ( globs->Type == 50 )
    {
      for ( code = 0; code < MaxCodes; code++ )
	if ( vars[code].selected )
	  {
	    if ( ! vars[code].spectral ) Error("Code %d not available on spectral space!", code);
	      
	    nlevel = zaxisInqSize(vars[code].ozaxisID);
	    for ( lindex = 0; lindex < nlevel; lindex++ )
	      {
		offset = lindex*globs->DimSP;
		streamDefRecord(globs->ostreamID, vars[code].ovarID, lindex);
		streamWriteRecord(globs->ostreamID, vars[code].spectral+offset, 0);
	      }	      
	  }

      FreeSpectral(vars);
      return;
    }

  /* =============================== */
  /* Transformation to fourier space */
  /* Computations in fourier space   */
  /* =============================== */

  if ( globs->Type >= 60 )
    {
      for (code = 0; code < MaxCodes; code++)
	if (vars[code].needed && vars[code].spectral)
	  {
	    if (vars[code].fourier == NULL)
	      {
		fieldSize = vars[code].plev * globs->DimFC;
		vars[code].fourier = alloc_dp(fieldSize, FieldName(code,"fourier"));
	      }
	    sp2fc(vars[code].spectral,vars[code].fourier,globs->poli,
		  vars[code].plev,globs->Latitudes,globs->Fouriers,globs->Truncation);
	  }
      if ( vars[U_WIND].needed && vars[U_WIND].fourier )
	scaluv(vars[U_WIND].fourier, globs->rcoslat, globs->Latitudes, globs->Fouriers*globs->NumLevelRequest);
      if ( vars[V_WIND].needed && vars[V_WIND].fourier )
	scaluv(vars[V_WIND].fourier, globs->rcoslat, globs->Latitudes, globs->Fouriers*globs->NumLevelRequest);

      /* HUMTEST */
      if ( globs->Type < 70 && vars[HUMIDITY].needed && vars[HUMIDITY].fourier == NULL )
	{
	  vars[HUMIDITY].plev = globs->NumLevelRequest;
	  vars[HUMIDITY].sfit = TRUE;
	  vars[HUMIDITY].fourier = alloc_dp(globs->DimFC*globs->NumLevelRequest, "vars[HUMIDITY].fourier");

	  after_FCrh2FCsh(globs, vars);

	  vars[RHUMIDITY].needed = vars[RHUMIDITY].selected;
	  vars[TEMPERATURE].needed = vars[TEMPERATURE].selected;
	}

      if ( globs->Type < 70 && vars[RHUMIDITY].needed && vars[RHUMIDITY].fourier == NULL )
	{
	  vars[RHUMIDITY].plev = globs->NumLevelRequest;
	  vars[RHUMIDITY].sfit = TRUE;
	  vars[RHUMIDITY].fourier = alloc_dp(globs->DimFC*globs->NumLevelRequest, "vars[RHUMIDITY].fourier");

	  after_FCsh2FCrh(globs, vars);

	  vars[HUMIDITY].needed = vars[HUMIDITY].selected;
	  vars[TEMPERATURE].needed = vars[TEMPERATURE].selected;
	}
      /* ENDE HUMTEST */
    }

  FreeSpectral(vars);

  /* -------------------------- */
  /*  Output of fourier fields  */
  /* -------------------------- */

  if ( globs->Type == 60 )
    {
      for ( code = 0; code < MaxCodes; code++ )
	if ( vars[code].selected )
	  {
	    if ( ! vars[code].fourier ) Error("Code %d not available on fourier space!", code);

	    nlevel = zaxisInqSize(vars[code].ozaxisID);
	    for ( lindex = 0; lindex < nlevel; lindex++ )
	      {
		offset = lindex*globs->DimFC;
		streamDefRecord(globs->ostreamID, vars[code].ovarID, lindex);
		streamWriteRecord(globs->ostreamID, vars[code].fourier+offset, 0);
	      }	      
	  }

      FreeFourier(vars);
      return;
    }

  /* ----------------------- */
  /*  Output of zonal means  */
  /* ----------------------- */

  if ( globs->Type == 61 )
    {
      for ( code = 0; code < MaxCodes; code++ )
	if ( vars[code].selected )
	  {
	    if ( ! vars[code].fourier ) Error("Code %d not available on zonal mean!", code);

	    nlevel = zaxisInqSize(vars[code].ozaxisID);
	    for ( lindex = 0; lindex < nlevel; lindex++ )
	      {
		offset = lindex*globs->DimFC;
		streamDefRecord(globs->ostreamID, vars[code].ovarID, lindex);
		streamWriteRecord(globs->ostreamID, vars[code].fourier+offset, 0);
	      }	      
	  }

      FreeFourier(vars);
      return;
   }

  /* ============================ */
  /* Transformation to gridpoints */
  /* ============================ */

  if ( vars[PS].comp && vars[LNPS].grid )
    {
      if (vars[PS].grid == NULL) vars[PS].grid = alloc_dp(globs->DimGP, "Ps");
      for (l = 0; l < globs->DimGP; l++) vars[PS].grid[l] = exp(vars[LNPS].grid[l]);
    }

  if ( globs->Type >= 70 )
    {
      for (code = 0; code < MaxCodes; code++)
	if (vars[code].needed && vars[code].fourier)
	  {
	    if (vars[code].grid == NULL)
	      {
		fieldSize = vars[code].plev * globs->DimGP;
		vars[code].grid = alloc_dp(fieldSize,FieldName(code,"grid"));
	      }

	    after_FC2GP(vars[code].fourier,vars[code].grid,
			globs->Latitudes,globs->Longitudes,vars[code].plev,globs->Fouriers);
	  }
    }

   FreeFourier(vars);

  /* HUMTEST */
  /* -------------------------------- */
  /*  Computation in gridpoint space  */
  /* -------------------------------- */

  if ( vars[RHUMIDITY].needed && vars[RHUMIDITY].grid == NULL )
    {
      vars[RHUMIDITY].plev = globs->NumLevelRequest;
      vars[RHUMIDITY].sfit = TRUE;
      vars[RHUMIDITY].grid = alloc_dp(globs->DimGP*globs->NumLevelRequest, "vars[RHUMIDITY].grid");
      sh2rh(globs->AnalysisData, vars[HUMIDITY].grid, vars[RHUMIDITY].grid, vars[TEMPERATURE].grid, globs->NumLevelRequest,
	    globs->DimGP, globs->LevelRequest, vars[FULL_PRESS].hybrid);
      vars[HUMIDITY].needed = vars[HUMIDITY].selected;
      vars[TEMPERATURE].needed = vars[TEMPERATURE].selected;
    }

  if ( vars[HUMIDITY].needed && vars[HUMIDITY].grid == NULL )
    {
      vars[HUMIDITY].plev = globs->NumLevelRequest;
      vars[HUMIDITY].sfit = TRUE;
      vars[HUMIDITY].grid = alloc_dp(globs->DimGP*globs->NumLevelRequest, "vars[HUMIDITY].grid");
      rh2sh(vars[HUMIDITY].grid, vars[RHUMIDITY].grid, vars[TEMPERATURE].grid, globs->NumLevelRequest,
	    globs->DimGP, globs->LevelRequest);
      vars[RHUMIDITY].needed = vars[RHUMIDITY].selected;
      vars[TEMPERATURE].needed = vars[TEMPERATURE].selected;
    }
  /* HUMTEST ENDE */

  /* -------------------------- */
  /*  Computation of Means      */
  /* -------------------------- */

  if ( globs->Mean )
    for ( code = 0; code < MaxCodes; code++ )
      if ( vars[code].needed && vars[code].grid )
	{
	  fieldSize = globs->DimGP * vars[code].plev;
	  if (vars[code].mean == NULL)
	    vars[code].mean = alloc_dp(fieldSize,FieldName(code,"mean"));

	  if (globs->MeanCount == 1)
	    after_copy_array(vars[code].mean, vars[code].grid, fieldSize);
	  else
	    AddVector(vars[code].mean, vars[code].grid, fieldSize,
		      &vars[code].nmiss, vars[code].missval);

	  if ( globs->EndOfInterval )
	    {
	      if ( vars[code].samp == NULL )
		MultVectorScalar(vars[code].mean, vars[code].mean, 1.0/globs->MeanCount, fieldSize,
				 vars[code].nmiss, vars[code].missval);
	      else
		DivVectorIvector(vars[code].mean, vars[code].mean, vars[code].samp, fieldSize,
				 &vars[code].nmiss, vars[code].missval);
	    }
	}

  /* -------------------------- */
  /*  Computation of Variances  */
  /* -------------------------- */

  if ( globs->Mean > 1 )
    for (code = 0; code < MaxCodes; code++)
      if (vars[code].needed && vars[code].mean)
	{
	  fieldSize = globs->DimGP * vars[code].plev;
	  if (vars[code].variance == NULL)
	    vars[code].variance = alloc_dp(fieldSize,FieldName(code,"var"));
	  if (globs->MeanCount == 1)
	    IniQuaSum(vars[code].variance,vars[code].grid,fieldSize);
	  else
	    AddQuaSum(vars[code].variance,vars[code].grid,fieldSize);

	  if (globs->EndOfInterval) VarQuaSum(vars[code].variance,vars[code].mean,fieldSize, globs->MeanCount);
	}

  if ( globs->Mean && !globs->EndOfInterval )
    {
      FreeGrid(vars);
      return;
    }

  /* ---------------------------------------------- */
  /*  Output of pressure level means and variances  */
  /* ---------------------------------------------- */

  if ( globs->Type == 70 && globs->Mean && globs->EndOfInterval )
    {
      for ( code = 0; code < MaxCodes; code++ )
	if ( vars[code].selected )
	  {
	    nlevel = zaxisInqSize(vars[code].ozaxisID);
	    for ( lindex = 0; lindex < nlevel; lindex++ )
	      {
		offset = lindex*globs->DimGP;
		if ( globs->Mean != 2 )
		  {
		    streamDefRecord(globs->ostreamID, vars[code].ovarID, lindex);
		    streamWriteRecord(globs->ostreamID, vars[code].mean+offset, vars[code].nmiss);
		  }
		if ( globs->Mean >= 2 )
		  {
		    streamDefRecord(globs->ostreamID2, vars[code].ovarID2, lindex);
		    streamWriteRecord(globs->ostreamID2, vars[code].variance+offset, vars[code].nmiss);
		  }
	      }
	  }

      FreeSamp(vars);
      FreeGrid(vars);
      return;
    }

  /* -------------------------------- */
  /*  Output of pressure level grids  */
  /* -------------------------------- */

  if ( globs->Type == 70 )
    {
      for ( code = 0; code < MaxCodes; code++ )
	if ( vars[code].selected )
	  {
	    nlevel = zaxisInqSize(vars[code].ozaxisID);
	    for ( lindex = 0; lindex < nlevel; lindex++ )
	      {
		offset = lindex*globs->DimGP;
		streamDefRecord(globs->ostreamID, vars[code].ovarID, lindex);
		streamWriteRecord(globs->ostreamID, vars[code].grid+offset, vars[code].nmiss);
	      }
	  }

      FreeGrid(vars);
      return;
    }
}

static
void theta(double *pthetaf, double *pthetah, double *ph, double *ps,
	   double *tf, double *ts, int levels, int dimgp, int dim3gp)
{
   int h,l;
   double  kappa;
   double *thetah = pthetah;
   double *thetaf = pthetaf;

   kappa = PlanetRD / C_RCPD;

   for (h = 0; h < dimgp; h++) thetah[h] = 0.0;
   thetah += dimgp;
   for (l = 0; l < levels - 1; l++) {
      for (h = 0; h < dimgp; h++) {
         thetah[h] = 0.5 * (tf[h] + tf[h+dimgp]) * pow((ps[h]/ph[h]),kappa);
      }
      ph += dimgp;
      tf += dimgp;
      thetah += dimgp;
   }
   after_copy_array(thetah,ts,dimgp);
   thetah = pthetah;
   for (h = 0; h < dim3gp; h++) {
      thetaf[h] = 0.5 * (thetah[h] + thetah[h+dimgp]);
   }
}

static
void windSpeed(double *uvspeed, double *u, double *v, int dim3gp)
{
  int i;

  for (i = 0; i < dim3gp; i++)
    uvspeed[i] = sqrt(u[i] * u[i] + v[i] * v[i]);
}

static
void Omega(double *omega_in, double *divergence, double *u_wind, double *v_wind,
	   double *halfpress, double *fullpress, double *dpsdx, double *dpsdy,
	   double *vct, int dimgp, int nlev)
{
  int i, j;
  double DeltaHybrid, Cterm, Pterm;
  double *diver, *halfp, *fullp, *uwind, *vwind;
  double *omega = omega_in;

  /* Compute half level part of vertical velocity */

  for ( i = 0; i < dimgp; i++ ) omega[i] = 0.0;

  for ( j = 0; j < nlev; j++ )
    {
      omega = omega_in   + j*dimgp;
      halfp = halfpress  + j*dimgp;
      diver = divergence + j*dimgp;
      uwind = u_wind     + j*dimgp;
      vwind = v_wind     + j*dimgp;

      DeltaHybrid = vct[nlev+j+2] - vct[nlev+j+1];
#if defined (SX)
#pragma vdir nodep
#endif
#if defined (__uxp__)
#pragma loop novrec
#endif
#if defined (_OPENMP)
#pragma omp parallel for
#endif
      for ( i = 0; i < dimgp; i++ )
	{
	  omega[i+dimgp] = omega[i]
	                 - diver[i] * (halfp[i+dimgp] - halfp[i])
	                 - DeltaHybrid * (uwind[i]*dpsdx[i] + vwind[i]*dpsdy[i]);
	}
    }

  /* interpolate to full levels  */

  for ( j = 0; j < nlev; j++ )
    {
      omega = omega_in   + j*dimgp;
#if defined (SX)
#pragma vdir nodep
#endif
#if defined (__uxp__)
#pragma loop novrec
#endif
#if defined (_OPENMP)
#pragma omp parallel for
#endif
      for ( i = 0; i < dimgp; i++ )
	omega[i] = 0.5 * (omega[i] + omega[i+dimgp]);
    }

  /* compute full level part of vertical velocity */

#if defined (_OPENMP)
#pragma omp parallel for default(shared) private(i, omega, halfp, fullp, uwind, vwind, DeltaHybrid, Cterm, Pterm)
#endif
  for ( j = 0; j < nlev; j++ )
    {
      omega = omega_in   + j*dimgp;
      halfp = halfpress  + j*dimgp;
      fullp = fullpress  + j*dimgp;
      uwind = u_wind     + j*dimgp;
      vwind = v_wind     + j*dimgp;

      DeltaHybrid = vct[nlev+j+2] - vct[nlev+j+1];
      if ( fabs(DeltaHybrid) > 0 )
	{
	  Cterm = vct[j+1]*vct[nlev+j+1] - vct[j]*vct[nlev+j+2];
#if defined (__uxp__)
#pragma loop novrec
#endif
	  for ( i = 0; i < dimgp; i++ )
	    {
	      if ( Cterm != 0.0 )
		Pterm = Cterm / (halfp[i+dimgp]-halfp[i]) * log(halfp[i+dimgp]/halfp[i]);
	      else
		Pterm = 0.0;

	      omega[i] += fullp[i] * (uwind[i]*dpsdx[i] + vwind[i]*dpsdy[i])
	 	        / (halfp[i+dimgp] - halfp[i]) * (DeltaHybrid + Pterm);
	    }
	}
    }
}


void MakeGeopotHeight(double *geop, double *gt, double *gq, double *ph, int nhor, int nlev)
{
  int i, j;
  double vtmp;
  double zrg;
  double z2log2;
  double *restrict geopl, *restrict gtl, *restrict gql, *restrict phl;

  z2log2 = 2.0 * log(2.0);
  vtmp   = (C_RV / PlanetRD) - 1.0;
  zrg    = 1.0 / PlanetGrav;

  if ( gq ) /* Humidity is present */
    {
      for ( j = nlev ; j > 1 ; j-- )
	{
	  geopl = geop + nhor*(j-1);
	  gtl   = gt   + nhor*(j-1);
	  gql   = gq   + nhor*(j-1);
	  phl   = ph   + nhor*(j-1);
#if defined (SX)
#pragma vdir nodep
#endif
#if defined (_OPENMP)
#pragma omp parallel for
#endif
	  for ( i = 0; i < nhor; i++ )
	    geopl[i] = geopl[i+nhor] + PlanetRD * gtl[i] * (1.0 + vtmp * gql[i])
	             * log(phl[i+nhor] / phl[i]);
	}

#if defined (SX)
#pragma vdir nodep
#endif
#if defined (_OPENMP)
#pragma omp parallel for
#endif
      for ( i = 0; i < nhor; i++ )
	geop[i] = geop[i+nhor] + PlanetRD * gt[i] * (1.0 + vtmp * gq[i]) * z2log2;
    }
  else    /* No humidity */
    {
      geopl = geop + nhor;
      phl   = ph   + nhor;
      
      for ( j = nlev ; j > 1 ; j-- )
#if defined (SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
        for ( i = nhor * (j-1) ; i < nhor * j ; i++ )
          geop[i] = geopl[i] + PlanetRD * gt[i] * log(phl[i] / ph[i]);

#if defined (SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
      for ( i = 0; i < nhor; i++ )
	geop[i] = geopl[i] + PlanetRD * gt[i] * z2log2;
    }

#if defined (SX)
#pragma vdir nodep
#endif
#if defined (_OPENMP)
#pragma omp parallel for
#endif
  for ( i = 0; i < nhor * (nlev+1); i++ ) geop[i] *= zrg;
}

#define  SCALESLP        (101325.0)

/* ======================================== */
/* LayerWater integral liquid water content */
/* ======================================== */

void LayerWater(double *ww, double *ll, double pmax, double pmin,
		int DimGP, int HalfLevels, double *vct)
{
   int  i,k;
   int  MaxLev, MinLev;
   double pph[MaxLevel];

   for (k = 0; k < HalfLevels; k++)
      pph[k] = vct[k] + vct[k+HalfLevels] * SCALESLP;
   for (k = 0; k < HalfLevels; k++)
      if (pph[k] > pmax) break;
   MaxLev = k - 1;
   for (k = HalfLevels - 1; k >= 0; k--)
      if (pph[k] < pmin) break;
   MinLev = k;

   after_zero_array (ll, DimGP);

   for (k = MaxLev; k <= MinLev; k++) {
     for (i = 0;     i < DimGP;  i++)
        ll[i] += ww[i+k*DimGP] * (pph[k+1] - pph[k]);
   }
   for (i = 0; i < DimGP; i++) ll[i] /= PlanetGrav;
}

/* ================================================= */
/* LayerCloud calculates random overlap cloud cover */
/* ================================================= */

void LayerCloud(double *cc, double *ll, double pmax, double pmin,
		int DimGP, int HalfLevels, double *vct)
{
   int  i, k;
   int  MaxLev, MinLev;
   double pph[MaxLevel];
   double ZEPSEC = 1.0e-12;

   for (k = 0; k < HalfLevels; k++)
      pph[k] = vct[k] + vct[k+HalfLevels] * SCALESLP;
   for (k = 0; k < HalfLevels; k++)
      if (pph[k] > pmax) break;
   MaxLev = k - 1;
   for (k  =  HalfLevels - 1; k >=0; k--)
      if (pph[k] < pmin) break;
   MinLev = k;

   for (i = 0; i < DimGP; i++) ll[i] = 1. - cc[i+MaxLev*DimGP];

   for (k = MaxLev + 1; k <= MinLev; k++) {
     for (i = 0;     i < DimGP;  i++)
         ll[i] *= (1. - MAX(cc[i+(k-1)*DimGP],cc[i+k*DimGP]))
                / (1. - MIN(cc[i+(k-1)*DimGP],1.-ZEPSEC));
   }
   for (i = 0; i < DimGP; i++) ll[i] = 1. - ll[i];
}

/* Grid Point Computations */
void after_EchamCompGP(struct Control *globs, struct Variable *vars)
{
  if ( vars[GEOPOTHEIGHT].comp || vars[SLP].comp || vars[THETAF].needed ||
       vars[HALF_PRESS].needed || vars[RHUMIDITY].comp || vars[OMEGA].comp ||
       globs->Type >= 30 )
    {
      if ( vars[FULL_PRESS].hybrid == NULL )
	vars[FULL_PRESS].hybrid = alloc_dp(globs->Dim3GP, "vars[FULL_PRESS].hybrid");

      vars[HALF_PRESS].hlev = globs->NumLevel+1;
      vars[HALF_PRESS].plev = globs->NumLevelRequest;
      vars[HALF_PRESS].sfit = FALSE;

      if ( vars[HALF_PRESS].hybrid == NULL )
	vars[HALF_PRESS].hybrid = alloc_dp(globs->Dim3GP+globs->DimGP, "vars[HALF_PRESS].hybrid");

      presh(vars[FULL_PRESS].hybrid, vars[HALF_PRESS].hybrid, globs->vct, vars[PS_PROG].hybrid,
	    globs->NumLevel, globs->DimGP);
    }

  if ( globs->unitsel > 2 ) vars[FULL_PRESS].hybrid = FreeMemory(vars[FULL_PRESS].hybrid);

  if ( vars[THETAF].needed )
    {
      vars[THETAF].hlev = globs->NumLevel;
      vars[THETAF].plev = globs->NumLevelRequest;
      vars[THETAF].sfit = TRUE;
      if ( vars[THETAF].hybrid == NULL )
	vars[THETAF].hybrid = alloc_dp(globs->Dim3GP, "vars[THETAF].hybrid");
      if ( vars[THETAH].hybrid == NULL )
	vars[THETAH].hybrid = alloc_dp(globs->Dim3GP, "vars[THETAH].hybrid");
      theta(vars[THETAF].hybrid, vars[THETAH].hybrid, vars[HALF_PRESS].hybrid, vars[PS_PROG].hybrid,
	    vars[TEMPERATURE].hybrid, vars[TS].hybrid, globs->NumLevel, globs->DimGP, globs->Dim3GP);
    }

  if ( vars[GEOPOTHEIGHT].comp )
    {
      vars[GEOPOTHEIGHT].hlev = globs->NumLevel+1;
      vars[GEOPOTHEIGHT].plev = globs->NumLevelRequest;
      vars[GEOPOTHEIGHT].sfit = TRUE;
      vars[GEOPOTHEIGHT].hybrid = alloc_dp(globs->Dim3GP+globs->DimGP, "vars[GEOPOTHEIGHT].hybrid");

      after_copy_array(vars[GEOPOTHEIGHT].hybrid+globs->Dim3GP, globs->Orography, globs->DimGP);
      MakeGeopotHeight(vars[GEOPOTHEIGHT].hybrid, vars[TEMPERATURE].hybrid,
		       vars[HUMIDITY].hybrid, vars[HALF_PRESS].hybrid, globs->DimGP, globs->NumLevel);

      vars[HUMIDITY].needed = vars[HUMIDITY].selected;
    }
  else if ( vars[GEOPOTHEIGHT].hybrid && vars[GEOPOTHEIGHT].hlev == globs->NumLevel )
    {
      vars[GEOPOTHEIGHT].hlev = globs->NumLevel+1;
      vars[GEOPOTHEIGHT].sfit = TRUE;
      vars[GEOPOTHEIGHT].hybrid = realloc(vars[GEOPOTHEIGHT].hybrid, (globs->Dim3GP+globs->DimGP)*sizeof(double));
      after_copy_array(vars[GEOPOTHEIGHT].hybrid+globs->Dim3GP, globs->Orography, globs->DimGP);
      for ( int i = 0; i < globs->DimGP; i++ ) vars[GEOPOTHEIGHT].hybrid[globs->Dim3GP+i] /= PlanetGrav;
    }

  if ( vars[DPSDX].needed || vars[DPSDY].needed )
    for ( int l = 0; l < globs->DimGP; l++ )
      {
	vars[DPSDX].hybrid[l] *= vars[PS_PROG].hybrid[l];
	vars[DPSDY].hybrid[l] *= vars[PS_PROG].hybrid[l];
      }

  if ( vars[OMEGA].comp )
    {
      vars[OMEGA].hlev = globs->NumLevel+1;
      vars[OMEGA].plev = globs->NumLevelRequest;
      vars[OMEGA].sfit = TRUE;
      vars[OMEGA].hybrid = alloc_dp(globs->Dim3GP+globs->DimGP, "OMEGA.hybrid");

      Omega(vars[OMEGA].hybrid, vars[DIVERGENCE].hybrid, vars[U_WIND].hybrid, vars[V_WIND].hybrid,
	    vars[HALF_PRESS].hybrid, vars[FULL_PRESS].hybrid, vars[DPSDX].hybrid, vars[DPSDY].hybrid,
	    globs->vct, globs->DimGP, globs->NumLevel);
      
      vars[DPSDX].needed = vars[DPSDX].selected;
      vars[DPSDY].needed = vars[DPSDY].selected;
    }

  if ( vars[WINDSPEED].comp )
    {
      vars[WINDSPEED].hlev = globs->NumLevel;
      vars[WINDSPEED].plev = globs->NumLevelRequest;
      vars[WINDSPEED].sfit = TRUE;
      vars[WINDSPEED].hybrid = alloc_dp(globs->Dim3GP, "vars[WINDSPEED].hybrid");

      windSpeed(vars[WINDSPEED].hybrid, vars[U_WIND].hybrid, vars[V_WIND].hybrid, globs->Dim3GP);
    }

  if ( vars[RHUMIDITY].comp )
    {
      vars[RHUMIDITY].hlev = globs->NumLevel;
      vars[RHUMIDITY].plev = globs->NumLevelRequest;
      vars[RHUMIDITY].sfit = FALSE;
      vars[RHUMIDITY].hybrid = alloc_dp(globs->Dim3GP, "vars[RHUMIDITY].hybrid");

      sh2rh(globs->AnalysisData, vars[HUMIDITY].hybrid, vars[RHUMIDITY].hybrid, vars[TEMPERATURE].hybrid, globs->NumLevel,
	    globs->DimGP, globs->LevelRequest, vars[FULL_PRESS].hybrid);

      vars[TEMPERATURE].needed = vars[TEMPERATURE].selected;
      vars[HUMIDITY].needed    = vars[HUMIDITY].selected;
    }

  if ( vars[PS].comp )
    {
      vars[PS].hlev = 1;
      vars[PS].plev = 1;
      vars[PS].sfit = TRUE; /* ??? */
      vars[PS].hybrid = alloc_dp(globs->DimGP, "vars[PS].hybrid");
      for ( int l = 0; l < globs->DimGP; l++ ) vars[PS].hybrid[l] = exp(vars[LNPS].hybrid[l]);
    }

  if ( vars[SLP].comp )
    {
      vars[SLP].hlev = 1;
      vars[SLP].plev = 1;
      vars[SLP].sfit = TRUE;
      vars[SLP].hybrid = alloc_dp(globs->DimGP, "vars[SLP].hybrid");

      extra_P(vars[SLP].hybrid, vars[HALF_PRESS].hybrid + globs->Dim3GP,
	      vars[FULL_PRESS].hybrid + globs->Dim3GP - globs->DimGP , globs->Orography,
	      vars[TEMPERATURE].hybrid + globs->Dim3GP - globs->DimGP , globs->DimGP);
      vars[TEMPERATURE].needed = vars[TEMPERATURE].selected || vars[GEOPOTHEIGHT].selected;
    }

  if ( vars[PRECIP].comp )
    {
      vars[PRECIP].hlev = vars[PRECIP].plev = 1;
      vars[PRECIP].sfit = FALSE;
      vars[PRECIP].hybrid = alloc_dp(globs->DimGP, "PRECIP.hybrid");
      Add2Vectors(vars[PRECIP].hybrid, vars[142].hybrid, vars[143].hybrid, globs->DimGP);
    }

  if ( vars[NET_TOP].comp )
    {
      vars[NET_TOP].hlev = vars[NET_TOP].plev = 1;
      vars[NET_TOP].sfit = FALSE;
      vars[NET_TOP].hybrid = alloc_dp(globs->DimGP, "NET_TOP.hybrid");
      Add2Vectors(vars[NET_TOP].hybrid, vars[178].hybrid, vars[179].hybrid, globs->DimGP);
    }

  if ( vars[NET_BOT].comp )
    {
      vars[NET_BOT].hlev = vars[NET_BOT].plev = 1;
      vars[NET_BOT].sfit = FALSE;
      vars[NET_BOT].hybrid = alloc_dp(globs->DimGP, "NET_BOT.hybrid");
      Add2Vectors(vars[NET_BOT].hybrid, vars[176].hybrid, vars[177].hybrid, globs->DimGP);
    }

  if ( vars[NET_HEAT].comp )
    {
      vars[NET_HEAT].hlev = vars[NET_HEAT].plev = 1;
      vars[NET_HEAT].sfit = FALSE;
      vars[NET_HEAT].hybrid = alloc_dp(globs->DimGP, "NET_HEAT.hybrid");
      /*
      if ( Source == S_ECHAM5 )
	{
	  MultVectorScalar(vars[NET_HEAT].hybrid, vars[218].hybrid, (-3.345e5), globs->DimGP,
	                   vars[218].nmiss, vars[218].missval);
	  Add2Vectors(vars[NET_HEAT].hybrid, vars[NET_HEAT].hybrid, vars[176].hybrid, globs->DimGP);
	  Add2Vectors(vars[NET_HEAT].hybrid, vars[NET_HEAT].hybrid, vars[177].hybrid, globs->DimGP);
	  Add2Vectors(vars[NET_HEAT].hybrid, vars[NET_HEAT].hybrid, vars[146].hybrid, globs->DimGP);
	  Add2Vectors(vars[NET_HEAT].hybrid, vars[NET_HEAT].hybrid, vars[147].hybrid, globs->DimGP);
	  Add2Vectors(vars[NET_HEAT].hybrid, vars[NET_HEAT].hybrid, vars[206].hybrid, globs->DimGP);
	  Sub2Vectors(vars[NET_HEAT].hybrid, vars[NET_HEAT].hybrid, vars[208].hybrid, globs->DimGP);
	  Sub2Vectors(vars[NET_HEAT].hybrid, vars[NET_HEAT].hybrid, vars[209].hybrid, globs->DimGP);
	}
      else
      */
	{
	  MultVectorScalar(vars[NET_HEAT].hybrid, vars[218].hybrid, C_TIMES_RHOH2O, globs->DimGP,
			   vars[218].nmiss, vars[218].missval);
	  Add2Vectors(vars[NET_HEAT].hybrid, vars[NET_HEAT].hybrid, vars[176].hybrid, globs->DimGP);
	  Add2Vectors(vars[NET_HEAT].hybrid, vars[NET_HEAT].hybrid, vars[177].hybrid, globs->DimGP);
	  Add2Vectors(vars[NET_HEAT].hybrid, vars[NET_HEAT].hybrid, vars[146].hybrid, globs->DimGP);
	  Add2Vectors(vars[NET_HEAT].hybrid, vars[NET_HEAT].hybrid, vars[147].hybrid, globs->DimGP);
	  Sub2Vectors(vars[NET_HEAT].hybrid, vars[NET_HEAT].hybrid, vars[220].hybrid, globs->DimGP);
	}
    }

  if ( vars[NET_WATER].comp )
    {
      vars[NET_WATER].hlev = vars[NET_WATER].plev = 1;
      vars[NET_WATER].sfit = FALSE;
      vars[NET_WATER].hybrid = alloc_dp(globs->DimGP, "NET_WATER.hybrid");
      Sub2Vectors(vars[NET_WATER].hybrid, vars[182].hybrid ,      vars[160].hybrid, globs->DimGP);
      Add2Vectors(vars[NET_WATER].hybrid, vars[NET_WATER].hybrid, vars[142].hybrid, globs->DimGP);
      Add2Vectors(vars[NET_WATER].hybrid, vars[NET_WATER].hybrid, vars[143].hybrid, globs->DimGP);
    }

  if ( vars[LOW_WATER].comp)
    {
      vars[LOW_WATER].hlev = vars[LOW_WATER].plev = 1;
      vars[LOW_WATER].sfit = FALSE;
      vars[LOW_WATER].hybrid = alloc_dp(globs->DimGP, "vars[LOW_WATER].hybrid");
      LayerWater(vars[222].hybrid, vars[LOW_WATER].hybrid, 75000.,101300., globs->DimGP, globs->HalfLevels, globs->vct);
    }

  if ( vars[MID_WATER].comp)
    {
      vars[MID_WATER].hlev = vars[MID_WATER].plev = 1;
      vars[MID_WATER].sfit = FALSE;
      vars[MID_WATER].hybrid = alloc_dp(globs->DimGP, "vars[MID_WATER].hybrid");
      LayerWater(vars[222].hybrid, vars[MID_WATER].hybrid, 46000., 73000., globs->DimGP, globs->HalfLevels, globs->vct);
    }

  if ( vars[HIH_WATER].comp )
    {
      vars[HIH_WATER].hlev = vars[HIH_WATER].plev = 1;
      vars[HIH_WATER].sfit = FALSE;
      vars[HIH_WATER].hybrid = alloc_dp(globs->DimGP, "vars[HIH_WATER].hybrid");
      LayerWater(vars[222].hybrid, vars[HIH_WATER].hybrid,  5000., 44000., globs->DimGP, globs->HalfLevels, globs->vct);
    }

  if ( vars[ALL_WATER].comp )
    {
      vars[ALL_WATER].hlev = vars[ALL_WATER].plev = 1;
      vars[ALL_WATER].sfit = FALSE;
      vars[ALL_WATER].hybrid = alloc_dp(globs->DimGP, "vars[ALL_WATER].hybrid");
      LayerWater(vars[222].hybrid, vars[ALL_WATER].hybrid,  5000.,101300., globs->DimGP, globs->HalfLevels, globs->vct);
    }

  if ( vars[LOW_CLOUD].comp )
    {
      vars[LOW_CLOUD].hlev = vars[LOW_CLOUD].plev = 1;
      vars[LOW_CLOUD].sfit = FALSE;
      vars[LOW_CLOUD].hybrid = alloc_dp(globs->DimGP, "vars[LOW_CLOUD].hybrid");
      LayerCloud(vars[223].hybrid, vars[LOW_CLOUD].hybrid, 75000.,101300., globs->DimGP, globs->HalfLevels, globs->vct);
    }

  if ( vars[MID_CLOUD].comp )
    {
      vars[MID_CLOUD].hlev = vars[MID_CLOUD].plev = 1;
      vars[MID_CLOUD].sfit = FALSE;
      vars[MID_CLOUD].hybrid = alloc_dp(globs->DimGP, "vars[MID_CLOUD].hybrid");
      LayerCloud(vars[223].hybrid, vars[MID_CLOUD].hybrid, 46000., 73000., globs->DimGP, globs->HalfLevels, globs->vct);
    }

  if ( vars[HIH_CLOUD].comp )
    {
      vars[HIH_CLOUD].hlev = vars[HIH_CLOUD].plev = 1;
      vars[HIH_CLOUD].sfit = FALSE;
      vars[HIH_CLOUD].hybrid = alloc_dp(globs->DimGP, "vars[HIH_CLOUD].hybrid");
      LayerCloud(vars[223].hybrid, vars[HIH_CLOUD].hybrid,  5000., 44000., globs->DimGP, globs->HalfLevels, globs->vct);
    }

  if ( vars[SW_CLF].comp )
    {
      vars[SW_CLF].hlev = vars[SW_CLF].plev = 1;
      vars[SW_CLF].sfit = FALSE;
      vars[SW_CLF].hybrid = alloc_dp(globs->DimGP, "SW_CLF.hybrid");
      Sub2Vectors(vars[SW_CLF].hybrid, vars[178].hybrid, vars[224].hybrid, globs->DimGP);
    }

  if ( vars[SW_BOT_CLF].comp )
    {
      vars[SW_BOT_CLF].hlev = vars[SW_BOT_CLF].plev = 1;
      vars[SW_BOT_CLF].sfit = FALSE;
      vars[SW_BOT_CLF].hybrid = alloc_dp(globs->DimGP, "vars[SW_BOT_CLF].hybrid");
      Sub2Vectors(vars[SW_BOT_CLF].hybrid, vars[176].hybrid, vars[185].hybrid, globs->DimGP);
    }

  if ( vars[SW_TOP_CLF].comp )
    {
      vars[SW_TOP_CLF].hlev = vars[SW_TOP_CLF].plev = 1;
      vars[SW_TOP_CLF].sfit = FALSE;
      vars[SW_TOP_CLF].hybrid = alloc_dp(globs->DimGP, "vars[SW_TOP_CLF].hybrid");
      Sub2Vectors(vars[SW_TOP_CLF].hybrid, vars[178].hybrid, vars[187].hybrid, globs->DimGP);
    }

  if ( vars[LW_CLF].comp )
    {
      vars[LW_CLF].hlev = vars[LW_CLF].plev = 1;
      vars[LW_CLF].sfit = FALSE;
      vars[LW_CLF].hybrid = alloc_dp(globs->DimGP, "LW_CLF.hybrid");
      Sub2Vectors(vars[LW_CLF].hybrid, vars[179].hybrid, vars[225].hybrid, globs->DimGP);
    }

  if ( vars[LW_BOT_CLF].comp )
    {
      vars[LW_BOT_CLF].hlev = vars[LW_BOT_CLF].plev = 1;
      vars[LW_BOT_CLF].sfit = FALSE;
      vars[LW_BOT_CLF].hybrid = alloc_dp(globs->DimGP, "vars[LW_BOT_CLF].hybrid");
      Sub2Vectors(vars[LW_BOT_CLF].hybrid, vars[177].hybrid, vars[186].hybrid, globs->DimGP);
    }

  if ( vars[LW_TOP_CLF].comp )
    {
      vars[LW_TOP_CLF].hlev = vars[LW_TOP_CLF].plev = 1;
      vars[LW_TOP_CLF].sfit = FALSE;
      vars[LW_TOP_CLF].hybrid = alloc_dp(globs->DimGP, "vars[LW_TOP_CLF].hybrid");
      Sub2Vectors(vars[LW_TOP_CLF].hybrid, vars[179].hybrid, vars[188].hybrid, globs->DimGP);
    }

  if ( vars[NET_CLF].comp )
    {
      vars[NET_CLF].hlev = vars[NET_CLF].plev = 1;
      vars[NET_CLF].sfit = FALSE;
      vars[NET_CLF].hybrid = alloc_dp(globs->DimGP, "NET_CLF.hybrid");
      Add2Vectors(vars[NET_CLF].hybrid, vars[178].hybrid, vars[179].hybrid, globs->DimGP);
      Sub2Vectors(vars[NET_CLF].hybrid, vars[NET_CLF].hybrid, vars[224].hybrid, globs->DimGP);
      Sub2Vectors(vars[NET_CLF].hybrid, vars[NET_CLF].hybrid, vars[225].hybrid, globs->DimGP);
    }

  if ( vars[SW_ATM].comp )
    {
      vars[SW_ATM].hlev = vars[SW_ATM].plev = 1;
      vars[SW_ATM].sfit = FALSE;
      vars[SW_ATM].hybrid = alloc_dp(globs->DimGP, "vars[SW_ATM].hybrid");
      Sub2Vectors(vars[SW_ATM].hybrid, vars[178].hybrid, vars[176].hybrid,globs->DimGP);
    }

  if ( vars[LW_ATM].comp )
    {
      vars[LW_ATM].hlev = vars[LW_ATM].plev = 1;
      vars[LW_ATM].sfit = FALSE;
      vars[LW_ATM].hybrid = alloc_dp(globs->DimGP, "vars[LW_ATM].hybrid");
      Sub2Vectors(vars[LW_ATM].hybrid, vars[179].hybrid, vars[177].hybrid,globs->DimGP);
    }

  if ( vars[NET_ATM].comp )
    {
      vars[NET_ATM].hlev = vars[NET_ATM].plev = 1;
      vars[NET_ATM].sfit = FALSE;
      vars[NET_ATM].hybrid = alloc_dp(globs->DimGP, "vars[NET_ATM].hybrid");
      Add2Vectors(vars[NET_ATM].hybrid, vars[178].hybrid, vars[179].hybrid,globs->DimGP);
      Sub2Vectors(vars[NET_ATM].hybrid, vars[NET_ATM].hybrid, vars[176].hybrid,globs->DimGP);
      Sub2Vectors(vars[NET_ATM].hybrid, vars[NET_ATM].hybrid, vars[177].hybrid,globs->DimGP);
    }

  if ( vars[SURF_RUNOFF].comp )
    {
      vars[SURF_RUNOFF].hlev = vars[SURF_RUNOFF].plev = 1;
      vars[SURF_RUNOFF].sfit = FALSE;
      vars[SURF_RUNOFF].hybrid = alloc_dp(globs->DimGP, "vars[SURF_RUNOFF].hybrid");
      Sub2Vectors(vars[SURF_RUNOFF].hybrid, vars[182].hybrid, vars[221].hybrid,globs->DimGP);
      Add2Vectors(vars[SURF_RUNOFF].hybrid, vars[SURF_RUNOFF].hybrid, vars[142].hybrid,globs->DimGP);
      Add2Vectors(vars[SURF_RUNOFF].hybrid, vars[SURF_RUNOFF].hybrid, vars[143].hybrid,globs->DimGP);
    }

  if ( vars[FRESH_WATER].comp )
    {
      vars[FRESH_WATER].hlev = vars[FRESH_WATER].plev = 1;
      vars[FRESH_WATER].sfit = FALSE;
      vars[FRESH_WATER].hybrid = alloc_dp(globs->DimGP, "vars[FRESH_WATER].hybrid");
      Add2Vectors(vars[FRESH_WATER].hybrid, vars[142].hybrid, vars[143].hybrid, globs->DimGP);
      Add2Vectors(vars[FRESH_WATER].hybrid, vars[FRESH_WATER].hybrid, vars[182].hybrid, globs->DimGP);
    }
}

static
void CheckContent(struct Variable *vars, int timestep)
{
  for ( int code = 0; code < 272; code++ )
    {
      /*  if ( code == GEOPOTENTIAL ) continue; */
      if ( code ==          SLP ) continue;
      if ( code == GEOPOTHEIGHT ) continue;
      if ( code ==       STREAM ) continue;
      if ( code ==      VELOPOT ) continue;
      if ( code ==       U_WIND ) continue;
      if ( code ==       V_WIND ) continue;
      if ( code ==        OMEGA ) continue;
      if ( code ==    RHUMIDITY ) continue;
      if ( code ==    LOW_CLOUD ) continue;
      if ( code ==    MID_CLOUD ) continue;
      if ( code ==    HIH_CLOUD ) continue;
      if ( code ==           PS ) continue;
      if ( code ==     HUMIDITY )
	{
	  if ( vars[code].needed && !vars[code].selected &&
	       vars[code].spectral == NULL &&
	       vars[code].hybrid   == NULL )
	    {
	      Warning( "No humidity in data file, set to zero !");
	      vars[code].needed = FALSE;
	    }
	}
      else
	{
	  if ( vars[code].needed && !vars[code].comp &&
	       vars[code].spectral == NULL &&
	       vars[code].hybrid   == NULL )
	    {
	      if ( labort_after )
		Error( "Code  %3d not found at timestep %d!", code, timestep);
	      else
		Warning( "Code  %3d not found at timestep %d!", code, timestep);
	    }
	}
    }
  /*
  if ( NumLevelRequest > 0 )
    {
      vars[HALF_PRESS].needed = 1;
      vars[FULL_PRESS].needed = 1;
    }

  code = HALF_PRESS;
  if ( vars[code].needed && !vars[code].comp &&
       vars[code].spectral == NULL && vars[code].hybrid == NULL )
    Error( "Hybrid model level not found!");

  code = FULL_PRESS;
  if ( vars[code].needed && !vars[code].comp &&
       vars[code].spectral == NULL && vars[code].hybrid == NULL )
    Error( "Hybrid model level not found!");
  */
}

static
void Derivate(double field[], double derilam[], int levels,
	      int Waves, int Latitudes, double DerivationFactor[])
{
   int l, n, lev;
   int i;

   i = 0;
   for (lev = 0; lev < levels; lev++)
   for (n = 0; n < Waves    ; n++) {
     for (l = 0; l < Latitudes; l++) {
       derilam[i] = -n * field[i+Latitudes] * DerivationFactor[l];
       i++;
     }
     for (l = 0; l < Latitudes; l++) {
       derilam[i] =  n * field[i-Latitudes] * DerivationFactor[l];
       i++;
     }
   }
}

/* Process Model Level data */
void after_processML(struct Control *globs, struct Variable *vars)
{
  int code,l,i;
  long fieldSize = 0;
  int lindex, nlevel;
  int offset;
  int leveltype;
  int nmiss;
  double *pressureLevel = NULL;

  globs->MeanCount++;
  globs->TermCount++;

  CheckContent(vars, globs->TermCount);

  if ( globs->MeanCount == 1 )
    {
      if ( globs->Debug ) Message( "TermCount = %d MeanCount = %d", globs->TermCount, globs->MeanCount);
      globs->StartDate = globs->OldDate;
    }

  if ( globs->TermCount  > 120 ) globs->Debug = 0;

  /* ============================== */
  /* Computations in spectral space */
  /* ============================== */

  if ( vars[U_WIND].comp || vars[V_WIND].comp )
    {
      vars[U_WIND].hlev = vars[V_WIND].hlev = vars[DIVERGENCE].hlev;
      vars[U_WIND].plev = vars[V_WIND].plev = vars[DIVERGENCE].plev;
      vars[U_WIND].sfit = vars[V_WIND].sfit = TRUE;
      vars[U_WIND].spectral = alloc_dp(globs->Dim3SP, "vars[U_WIND].spectral");
      vars[V_WIND].spectral = alloc_dp(globs->Dim3SP, "vars[V_WIND].spectral");

      if ( vars[DIVERGENCE].spectral == NULL ) after_gp2sp(globs, vars, DIVERGENCE);
      if ( vars[VORTICITY].spectral == NULL )  after_gp2sp(globs, vars, VORTICITY);

      dv2uv(vars[DIVERGENCE].spectral, vars[VORTICITY].spectral,
	    vars[U_WIND].spectral, vars[V_WIND].spectral,
	    globs->dv2uv_f1, globs->dv2uv_f2,
	    globs->Truncation, globs->DimSP, vars[DIVERGENCE].hlev);
    }

   if ( vars[VELOPOT].comp && globs->Type < 30 )
     {
       vars[VELOPOT].hlev = vars[DIVERGENCE].hlev;
       vars[VELOPOT].plev = vars[DIVERGENCE].plev;
       vars[VELOPOT].spectral = alloc_dp(globs->Dim3SP, "vars[VELOPOT].spectral");

       if ( vars[DIVERGENCE].spectral == NULL ) after_gp2sp(globs, vars, DIVERGENCE);

       dv2ps(vars[DIVERGENCE].spectral, vars[VELOPOT].spectral, vars[DIVERGENCE].hlev, globs->Truncation);
     }

   if ( vars[STREAM].comp && globs->Type < 30 )
     {
       vars[STREAM].hlev = vars[VORTICITY].hlev;
       vars[STREAM].plev = vars[VORTICITY].plev;
       vars[STREAM].spectral = alloc_dp(globs->Dim3SP, "vars[STREAM].spectral");

       if ( vars[VORTICITY].spectral == NULL ) after_gp2sp(globs, vars, VORTICITY);

       dv2ps(vars[VORTICITY].spectral, vars[STREAM].spectral,
	     vars[VORTICITY].hlev, globs->Truncation);
     }

   if ( vars[VORTICITY].spectral && !vars[VORTICITY].needed )
        vars[VORTICITY].spectral = FreeMemory(vars[VORTICITY].spectral);

   if ( vars[DIVERGENCE].spectral && !vars[DIVERGENCE].needed )
        vars[DIVERGENCE].spectral = FreeMemory(vars[DIVERGENCE].spectral);

  /* --------------------------- */
  /*  Output of spectral fields  */
  /* --------------------------- */

  if ( globs->Type == 0 )
    {
      for ( code = 0; code < MaxCodes; code++ )
	if ( vars[code].selected )
	  {
	    if ( ! vars[code].spectral ) Error("Code %d not available on spectral space!", code);

	    nlevel = zaxisInqSize(vars[code].ozaxisID);
	    for ( lindex = 0; lindex < nlevel; lindex++ )
	      {
		offset = lindex*globs->DimSP;
		streamDefRecord(globs->ostreamID, vars[code].ovarID, lindex);
		streamWriteRecord(globs->ostreamID, vars[code].spectral+offset, 0);
	      }
	  }

      FreeSpectral(vars);
      return;
    }

  /* ------------------------------- */
  /*  Computations in fourier space  */
  /* ------------------------------- */

  if ( globs->Type >= 10 )
    {
      for ( code = 0; code < MaxCodes; code++ )
	if ( vars[code].needed && vars[code].spectral )
	  {
	    if ( vars[code].fourier == NULL )
	      {
		fieldSize = vars[code].hlev * globs->DimFC;
		vars[code].fourier = alloc_dp(fieldSize,
						FieldName(code,"fourier"));
		sp2fc(vars[code].spectral, vars[code].fourier, globs->poli,
		      vars[code].hlev, globs->Latitudes, globs->Fouriers, globs->Truncation);
	      }
	    if ( code != LNPS )
	      vars[code].spectral = FreeMemory(vars[code].spectral);
	  }

      /*    if (globs->Type < 60) globs->poli = FreeMemory(globs->poli); */
      /*    if (globs->Type < 50) globs->pol2 = FreeMemory(globs->pol2); */
      /*    if (globs->Type < 50) globs->pol3 = FreeMemory(globs->pol3); */

      if ( vars[U_WIND].needed && vars[U_WIND].fourier )
	scaluv(vars[U_WIND].fourier, globs->rcoslat, globs->Latitudes, globs->Fouriers*globs->NumLevel);
      if ( vars[V_WIND].needed && vars[V_WIND].fourier )
	scaluv(vars[V_WIND].fourier, globs->rcoslat, globs->Latitudes, globs->Fouriers*globs->NumLevel);

      if ( vars[DPSDX].needed )
	{
	  vars[DPSDX].hlev = 1;
	  vars[DPSDX].plev = 1;
	  vars[DPSDX].sfit = FALSE;
	  vars[DPSDX].fourier = alloc_dp(globs->DimFC, "vars[DPSDX].fourier");
	  if ( vars[LNPS].fourier == NULL )  after_gp2sp(globs, vars, LNPS);
	  Derivate(vars[LNPS].fourier, vars[DPSDX].fourier, 1, globs->Waves, globs->Latitudes, globs->DerivationFactor);
	}
      if ( vars[DPSDY].needed )
	{
	  vars[DPSDY].hlev = 1;
	  vars[DPSDY].plev = 1;
	  vars[DPSDY].sfit = FALSE;
	  vars[DPSDY].fourier = alloc_dp(globs->DimFC, "vars[DPSDY].fourier");
	  if ( vars[LNPS].spectral == NULL )  after_gp2sp(globs, vars, LNPS);
	  sp2fc(vars[LNPS].spectral, vars[DPSDY].fourier, globs->pdev,
		vars[DPSDY].hlev, globs->Latitudes, globs->Fouriers, globs->Truncation);
	}
    }

  FreeSpectral(vars);

  /* -------------------------- */
  /*  Output of fourier fields  */
  /* -------------------------- */

  if ( globs->Type == 10 )
    {
      for ( code = 0; code < MaxCodes; code++ )
	if ( vars[code].selected )
	  {
	    if ( ! vars[code].fourier ) Error("Code %d not available on fourier space!", code);
	    
	    nlevel = zaxisInqSize(vars[code].ozaxisID);
	    for ( lindex = 0; lindex < nlevel; lindex++ )
	      {
		offset = lindex*globs->DimFC;
		streamDefRecord(globs->ostreamID, vars[code].ovarID, lindex);
		streamWriteRecord(globs->ostreamID, vars[code].fourier+offset, 0);
	      }
	  }

      FreeFourier(vars);
      return;
    }

  /* ----------------------- */
  /*  Output of zonal means  */
  /* ----------------------- */

  if ( globs->Type == 11 )
    {
      for ( code = 0; code < MaxCodes; code++ )
	if ( vars[code].selected )
	  {
	    if ( ! vars[code].fourier ) Error("Code %d not available on zonal mean!", code);

	    nlevel = zaxisInqSize(vars[code].ozaxisID);
	    for ( lindex = 0; lindex < nlevel; lindex++ )
	      {
		offset = lindex*globs->DimFC;
		streamDefRecord(globs->ostreamID, vars[code].ovarID, lindex);
		streamWriteRecord(globs->ostreamID, vars[code].fourier+offset, 0);
	      }
	  }

      FreeFourier(vars);
      return;
    }

  /* ------------------------------ */
  /*  Transformation to gridpoints  */
  /* ------------------------------ */

  if ( globs->Type >= 20 )
    {
      for ( code = 0; code < MaxCodes; code++ )
	if ( vars[code].needed && vars[code].fourier )
	  {
	    if ( vars[code].hybrid == NULL )
	      {
		fieldSize = globs->DimGP * vars[code].hlev;
		vars[code].hybrid = alloc_dp(fieldSize,
					       FieldName(code,"hybrid"));
		after_FC2GP(vars[code].fourier,vars[code].hybrid,
			    globs->Latitudes,globs->Longitudes,vars[code].hlev,globs->Fouriers);
	      }
	    vars[code].fourier = FreeMemory(vars[code].fourier);
	  }

      if ( vars[PS_PROG].comp && vars[PS_PROG].hybrid == NULL )
	{
	  vars[PS_PROG].hybrid = alloc_dp(globs->DimGP, "PS_PROG");
	  if ( vars[LNPS].hybrid )
	    {
	      for (l = 0; l < globs->DimGP; l++) vars[PS_PROG].hybrid[l] = exp(vars[LNPS].hybrid[l]);
	    }
	  else if ( vars[PS].hybrid )
	    {
	      Warning("log surface pressure (code 152) not found - using surface pressure (code 134)!");
	      after_copy_array(vars[PS_PROG].hybrid, vars[PS].hybrid, globs->DimGP);
	    }
	  else
	    {
	      Error("surface pressure not found!");
	    }
	}
      vars[LNPS].needed = vars[LNPS].selected;
      
      if ( globs->Orography == NULL )
	{
	  globs-> Orography = alloc_dp(globs->DimGP , "Orography");
	  if ( vars[GEOPOTENTIAL].hybrid )
	    after_copy_array(globs->Orography, vars[GEOPOTENTIAL].hybrid, globs->DimGP);
	  else
	    {
	      if ( vars[GEOPOTENTIAL].selected || globs->Type >= 30 )
		{
		  Warning("Orography not found - using zero orography!");
		  after_zero_array(globs->Orography, globs->DimGP);
		}
	    }
	}
      vars[GEOPOTENTIAL].needed = vars[GEOPOTENTIAL].selected;

      after_EchamCompGP(globs, vars);
    }

  FreeFourier(vars);

  if ( globs->Type == 20 )
    {
      /* ----------------------- */
      /*  Means on hybrid grids  */
      /* ----------------------- */

      if ( globs->Mean )
	{
	  for ( code = 0; code < MaxCodes; code++ )
	    {
	      if ( vars[code].selected && vars[code].hybrid )
		{
		  fieldSize = globs->DimGP * vars[code].hlev;

		  if ( vars[code].mean == NULL )
		    vars[code].mean = alloc_dp(fieldSize, FieldName(code, "mean"));

		  if ( globs->Mean > 1 && vars[code].variance == NULL )
		    vars[code].variance = alloc_dp(fieldSize, FieldName(code, "variance"));

		  if ( globs->MeanCount == 1 )
		    {
		      after_copy_array(vars[code].mean, vars[code].hybrid, fieldSize);
		      if ( globs->Mean > 1 )
			IniQuaSum(vars[code].variance, vars[code].hybrid, fieldSize);
		    }
		  else
		    {
		      AddVector(vars[code].mean, vars[code].hybrid, fieldSize,
				&vars[code].nmiss, vars[code].missval);
		      if ( globs->Mean > 1 )
			AddQuaSum(vars[code].variance, vars[code].hybrid, fieldSize);
		    }

		  if ( globs->EndOfInterval )
		    {
		      if ( vars[code].samp == NULL )
			MultVectorScalar(vars[code].hybrid, vars[code].mean, 1.0/globs->MeanCount, fieldSize,
					 vars[code].nmiss, vars[code].missval);
		      else
			DivVectorIvector(vars[code].hybrid, vars[code].mean, vars[code].samp, fieldSize,
					 &vars[code].nmiss, vars[code].missval);

		      if ( globs->Mean > 1 )
			VarQuaSum(vars[code].variance, vars[code].hybrid, fieldSize, globs->MeanCount);
		    }
		}
	    }
	}

      /* ---------------------------- */
      /* Output of hybrid level grids */
      /* ---------------------------- */

      if ( globs->Mean == 0 || globs->EndOfInterval )
	{
	  for ( code = 0; code < MaxCodes; code++ )
	    if ( vars[code].selected )
	      {
		if ( vars[code].hybrid == NULL )
		  Error("Internal problem. Code %d not allocated!", code);

		nlevel = zaxisInqSize(vars[code].ozaxisID);
		for ( lindex = 0; lindex < nlevel; lindex++ )
		  {
		    offset = lindex*globs->DimGP;
		    if ( globs->Mean != 2 )
		      {
			streamDefRecord(globs->ostreamID, vars[code].ovarID, lindex);
			streamWriteRecord(globs->ostreamID, vars[code].hybrid+offset, vars[code].nmiss);
		      }
		    if ( globs->Mean >= 2 )
		      {
			streamDefRecord(globs->ostreamID2, vars[code].ovarID2, lindex);
			streamWriteRecord(globs->ostreamID2, vars[code].variance+offset, vars[code].nmiss);
		      }
		  }
	      }

	  FreeSamp(vars);
	}

      FreeHybrid(vars);
      return;
    }

  /* -------------------------------------- */
  /* Vertical interpolation / extrapolation */
  /* -------------------------------------- */

  if ( globs->Type >= 30 )
    {
      if ( globs->vert_index == NULL )
	globs->vert_index = (int *) malloc(globs->NumLevelRequest*globs->DimGP*sizeof(int));

      if ( globs->unitsel )
	{
	  if ( globs->p_of_height == NULL ) globs->p_of_height = alloc_dp(globs->NumLevelRequest,"p_of_height");
	  height2pressure(globs->p_of_height, globs->LevelRequest, globs->NumLevelRequest);
	  pressureLevel = globs->p_of_height;
	}
      else
	{
	  pressureLevel = globs->LevelRequest;
	}

      genind(globs->vert_index, pressureLevel, vars[FULL_PRESS].hybrid, globs->DimGP, globs->NumLevelRequest, globs->NumLevel);

      nmiss = 0;
      if ( ! globs->Extrapolate )
	{
	  if ( globs->pnmiss == NULL ) globs->pnmiss = (int *) malloc(globs->NumLevelRequest*sizeof(int));  
	  genindmiss(globs->vert_index, pressureLevel, globs->DimGP, globs->NumLevelRequest, vars[PS_PROG].hybrid, globs->pnmiss);
	  for ( i = 0; i < globs->NumLevelRequest; i++ ) nmiss += globs->pnmiss[i];
	}

      for ( code = 0; code < MaxCodes; code++ )
	if ( vars[code].needed && vars[code].hybrid )
	  {
	    leveltype = zaxisInqType(vars[code].izaxisID);
	    if ( vars[code].hlev == 1 || leveltype != ZAXIS_HYBRID || (vars[code].hlev < globs->NumLevel) )
	      {
		if ( vars[code].grid ) FreeMemory(vars[code].grid);
		vars[code].grid = vars[code].hybrid;
		vars[code].hybrid = NULL;
	      }
	    else
	      {
		if ( vars[code].grid == NULL )
		  {
		    fieldSize = globs->DimGP * globs->NumLevelRequest;
		    vars[code].grid = alloc_dp(fieldSize, FieldName(code, "grid"));
		  }

		if ( code == TEMPERATURE )
		  {
		    interp_T(globs->Orography, vars[TEMPERATURE].hybrid, vars[TEMPERATURE].grid,
			     vars[FULL_PRESS].hybrid, vars[HALF_PRESS].hybrid, globs->vert_index,
			     pressureLevel, globs->NumLevelRequest, globs->DimGP, globs->NumLevel, vars[code].missval);
		  }
		else if ( code == GEOPOTHEIGHT )
		  {
		    if ( vars[TEMPERATURE].hybrid == NULL ) Error("Code  130 not found!");
		    interp_Z(globs->Orography, vars[GEOPOTHEIGHT].hybrid, vars[GEOPOTHEIGHT].grid,
			     vars[FULL_PRESS].hybrid, vars[HALF_PRESS].hybrid, globs->vert_index,
			     vars[TEMPERATURE].hybrid,
			     pressureLevel, globs->NumLevelRequest, globs->DimGP, globs->NumLevel, vars[code].missval);
		  }
		else
		  {
		    interp_X(vars[code].hybrid, vars[code].grid, vars[FULL_PRESS].hybrid,
			     globs->vert_index, pressureLevel, globs->NumLevelRequest, globs->DimGP, globs->NumLevel, vars[code].missval);
		  }

		if ( ! globs->Extrapolate ) vars[code].nmiss = nmiss;

		if ( code != TEMPERATURE )
		  vars[code].hybrid = FreeMemory(vars[code].hybrid);
	      }
	  }
    }

  vars[TEMPERATURE].needed = vars[TEMPERATURE].selected;
  FreeHybrid(vars);
  if ( vars[HALF_PRESS].hybrid )
    vars[HALF_PRESS].hybrid = FreeMemory(vars[HALF_PRESS].hybrid);

  /* -------------------------------- */
  /*  Output of pressure level grids  */
  /* -------------------------------- */

  if ( globs->Type == 30 && globs->Mean == 0 )
    {
      for ( code = 0; code < MaxCodes; code++ )
	if ( vars[code].selected && vars[code].grid )
	  {
	    nlevel = zaxisInqSize(vars[code].ozaxisID);
	    for ( lindex = 0; lindex < nlevel; lindex++ )
	      {
		offset = lindex*globs->DimGP;
		if ( globs->Mean != 2 )
		  {
		    streamDefRecord(globs->ostreamID, vars[code].ovarID, lindex);
		    streamWriteRecord(globs->ostreamID, vars[code].grid+offset, vars[code].nmiss);
		  }
	      }
	  }
      
      FreeGrid(vars);
      return;
    }

  /* ---------------------- */
  /*  Computation of Means  */
  /* ---------------------- */

  if ( globs->Type >= 30 && globs->Mean )
    for ( code = 0; code < MaxCodes; code++ )
      if ( vars[code].needed && vars[code].grid )
	{
	  fieldSize = globs->DimGP * vars[code].plev;

	  if ( vars[code].mean == NULL )
	    vars[code].mean = alloc_dp(fieldSize, FieldName(code, "mean"));

	  if ( globs->MeanCount == 1 )
	    after_copy_array(vars[code].mean, vars[code].grid, fieldSize);
	  else
	    AddVector(vars[code].mean, vars[code].grid, fieldSize,
		      &vars[code].nmiss, vars[code].missval);

	  if ( globs->EndOfInterval )
	    {
	      if ( vars[code].samp == NULL )
		MultVectorScalar(vars[code].mean, vars[code].mean, 1.0/globs->MeanCount, fieldSize,
				 vars[code].nmiss, vars[code].missval);
	      else
		DivVectorIvector(vars[code].mean, vars[code].mean, vars[code].samp, fieldSize,
				 &vars[code].nmiss, vars[code].missval);
	    }
	}

  /* -------------------------- */
  /*  Computation of Variances  */
  /* -------------------------- */

  if ( globs->Type >= 30 && globs->Mean > 1 )
    for ( code = 0; code < MaxCodes; code++ )
      if ( vars[code].needed && vars[code].mean )
	{
	  fieldSize = globs->DimGP * vars[code].plev;

	  if (vars[code].variance == NULL)
	    vars[code].variance = alloc_dp(fieldSize, FieldName(code, "var"));

	  if ( globs->MeanCount == 1 )
	    IniQuaSum(vars[code].variance, vars[code].grid, fieldSize);
	  else
	    AddQuaSum(vars[code].variance, vars[code].grid, fieldSize);

	  if ( globs->EndOfInterval )
	    VarQuaSum(vars[code].variance, vars[code].mean, fieldSize, globs->MeanCount);
	}

   if ( globs->Mean && !globs->EndOfInterval )
     {
       FreeGrid(vars);
       return;
     }

  /* --------------------------------------------- */
  /*  Output of pressure level means and variances */
  /* --------------------------------------------- */

  if ( globs->Type == 30 && globs->Mean )
    {
      for ( code = 0; code < MaxCodes; code++ )
	if ( vars[code].selected && vars[code].mean )
	  {
	    nlevel = zaxisInqSize(vars[code].ozaxisID);
	    for ( lindex = 0; lindex < nlevel; lindex++ )
	      {
		offset = lindex*globs->DimGP;
		if ( globs->Mean != 2 )
		  {
		    streamDefRecord(globs->ostreamID, vars[code].ovarID, lindex);
		    streamWriteRecord(globs->ostreamID, vars[code].mean+offset, vars[code].nmiss);
		  }
		if ( globs->Mean >= 2 )
		  {
		    streamDefRecord(globs->ostreamID2, vars[code].ovarID2, lindex);
		    streamWriteRecord(globs->ostreamID2, vars[code].variance+offset, vars[code].nmiss);
		  }
	      }
	  }

      FreeSamp(vars);
      FreeGrid(vars);
      return;
    }


  /* ------------------ */
  /*  Free mean fields  */
  /* ------------------ */

  if ( globs->Type >= 40 && globs->Mean )
    for ( code = 0; code < MaxCodes; code++ )
      if ( vars[code].mean )
	{
	  if ( vars[code].variance ) vars[code].variance = FreeMemory(vars[code].variance);
	  if ( vars[code].grid )     vars[code].grid     = FreeMemory(vars[code].grid);
	  vars[code].grid = vars[code].mean;
	  vars[code].mean = NULL;
	}

  /* --------------------------------- */
  /*  Transformation to fourier space  */
  /* --------------------------------- */

  if ( globs->Type >= 40 )
    {
      for ( code = 0; code < MaxCodes; code++ )
	if ( vars[code].needed && vars[code].grid && (vars[code].sfit || globs->Type < 70) )
	  {
	    if ( vars[code].nmiss > 0 )
	      Error("Missing values for code %d unsupported with TYPE > 30!", code);

	    if ( vars[code].fourier == NULL )
	      {
		fieldSize = globs->DimFC * vars[code].plev;
		vars[code].fourier = alloc_dp(fieldSize,FieldName(code,"fourier"));
	      }

	    after_GP2FC(vars[code].grid, vars[code].fourier,
			globs->Latitudes, globs->Longitudes, vars[code].plev, globs->Fouriers);

	    if ( vars[code].grid && (vars[code].sfit || globs->Type < 70) )
	      vars[code].grid = FreeMemory(vars[code].grid);
	  }
    }

  for (code = 0; code < MaxCodes; code++)
    if (vars[code].grid && (vars[code].sfit || globs->Type < 70))
      vars[code].grid = FreeMemory(vars[code].grid);

  /* -------------------------- */
  /*  Output of fourier fields  */
  /* -------------------------- */

  if ( globs->Type == 40 )
    {
      for ( code = 0; code < MaxCodes; code++ )
	if ( vars[code].selected )
	  {
	    if ( ! vars[code].fourier ) Error("Code %d not available on fourier space!", code);
	    
	    nlevel = zaxisInqSize(vars[code].ozaxisID);
	    for ( lindex = 0; lindex < nlevel; lindex++ )
	      {
		offset = lindex*globs->DimFC;
		streamDefRecord(globs->ostreamID, vars[code].ovarID, lindex);
		streamWriteRecord(globs->ostreamID, vars[code].fourier+offset, vars[code].nmiss);
	      }
	  }

      FreeFourier(vars);
      return;
    }

  /* --------------------- */
  /* Output of zonal means */
  /* --------------------- */

  if ( globs->Type == 41 )
    {
      for ( code = 0; code < MaxCodes; code++ )
	if ( vars[code].selected )
	  {
	    if ( ! vars[code].fourier ) Error("Code %d not available on zonal mean!", code);

	    nlevel = zaxisInqSize(vars[code].ozaxisID);
	    for ( lindex = 0; lindex < nlevel; lindex++ )
	      {
		offset = lindex*globs->DimFC;
		streamDefRecord(globs->ostreamID, vars[code].ovarID, lindex);
		streamWriteRecord(globs->ostreamID, vars[code].fourier+offset, vars[code].nmiss);
	      }
	  }

      FreeFourier(vars);
      return;
    }

  /* ---------------------------------- */
  /*  Transformation to spectral space  */
  /* ---------------------------------- */

  if ( globs->Type >= 50 )
    {
      if ( vars[U_WIND].needed && vars[U_WIND].fourier )
	scaluv(vars[U_WIND].fourier, globs->coslat, globs->Latitudes, globs->Fouriers*globs->NumLevelRequest);
      if ( vars[V_WIND].needed && vars[V_WIND].fourier )
	scaluv(vars[V_WIND].fourier, globs->coslat, globs->Latitudes, globs->Fouriers*globs->NumLevelRequest);

      for ( code = 0; code < MaxCodes; code++ )
	if ( vars[code].needed && vars[code].fourier )
	  {
	    if ( vars[code].spectral == NULL )
	      {
		fieldSize = vars[code].plev * globs->DimSP;
		vars[code].spectral = alloc_dp(fieldSize, FieldName(code,"spectral"));
	      }

	    fc2sp(vars[code].fourier,vars[code].spectral,
		  globs->pold,vars[code].plev,globs->Latitudes,globs->Fouriers,globs->Truncation);
	  }

      if ( vars[DIVERGENCE].needed || vars[VORTICITY].needed ||
	   vars[VELOPOT].needed    || vars[STREAM].needed )
	{
	  if ( vars[DIVERGENCE].spectral == NULL )
	    vars[DIVERGENCE].spectral = alloc_dp(globs->DimSP*globs->NumLevelRequest,
						   "vars[DIVERGENCE].spectral");
	  if ( vars[VORTICITY].spectral == NULL )
	    vars[VORTICITY].spectral = alloc_dp(globs->DimSP*globs->NumLevelRequest,
						  "vars[VORTICITY].spectral");
	  if ( (vars[U_WIND].fourier == 0 || vars[V_WIND].fourier == 0) && globs->NumLevelRequest )
	    Error("uwind or vwind missing!");
	  uv2dv(vars[U_WIND].fourier, vars[V_WIND].fourier,
		vars[DIVERGENCE].spectral, vars[VORTICITY].spectral,
		globs->pol2, globs->pol3, globs->NumLevelRequest, globs->Latitudes, globs->Truncation);
	}

      if ( vars[VELOPOT].needed )
	{
	  vars[VELOPOT].hlev = vars[DIVERGENCE].hlev;
	  vars[VELOPOT].plev = vars[DIVERGENCE].plev;
	  vars[VELOPOT].sfit = TRUE;
	  if ( vars[VELOPOT].spectral == NULL )
	    vars[VELOPOT].spectral = alloc_dp(globs->DimSP*globs->NumLevelRequest, "vars[VELOPOT].spectral");
	  dv2ps(vars[DIVERGENCE].spectral, vars[VELOPOT].spectral, globs->NumLevelRequest, globs->Truncation);
	}

      if ( vars[STREAM].needed )
	{
	  vars[STREAM].hlev = vars[VORTICITY].hlev;
	  vars[STREAM].plev = vars[VORTICITY].plev;
	  vars[STREAM].sfit = TRUE;
	  if ( vars[STREAM].spectral == NULL )
	    vars[STREAM].spectral = alloc_dp(globs->DimSP*globs->NumLevelRequest, "vars[STREAM].spectral");
	  dv2ps(vars[VORTICITY].spectral, vars[STREAM].spectral, globs->NumLevelRequest, globs->Truncation);
	}
    }

  for ( code = 0; code < MaxCodes; code++ )
    if ( vars[code].fourier && (vars[code].sfit || globs->Type < 61) )
      vars[code].fourier = FreeMemory(vars[code].fourier);

  /* --------------------------- */
  /*  Output of spectral fields  */
  /* --------------------------- */

  if ( globs->Type == 50 )
    {
      for ( code = 0; code < MaxCodes; code++ )
	if ( vars[code].selected && vars[code].spectral )
	  {
	    nlevel = zaxisInqSize(vars[code].ozaxisID);
	    for ( lindex = 0; lindex < nlevel; lindex++ )
	      {
		offset = lindex*globs->DimSP;
		streamDefRecord(globs->ostreamID, vars[code].ovarID, lindex);
		streamWriteRecord(globs->ostreamID, vars[code].spectral+offset, 0);
	      }	      
	  }

      FreeSpectral(vars);
      return;
    }

  /* -------------------------------*/
  /*  Computations in fourier space */
  /* -------------------------------*/

  if ( globs->Type >= 60 )
    {
      for (code = 0; code < MaxCodes; code++)
      if (vars[code].needed && vars[code].spectral)
	{
	  if (vars[code].fourier == NULL)
	    {
	      fieldSize = vars[code].plev * globs->DimFC;
	      vars[code].fourier = alloc_dp(fieldSize, FieldName(code,"fourier"));
	    }
	  sp2fc(vars[code].spectral,vars[code].fourier,globs->poli,
		vars[code].plev,globs->Latitudes,globs->Fouriers,globs->Truncation);
	}
      if ( vars[U_WIND].needed && vars[U_WIND].fourier )
	scaluv(vars[U_WIND].fourier, globs->rcoslat, globs->Latitudes, globs->Fouriers*globs->NumLevelRequest);
      if ( vars[V_WIND].needed && vars[V_WIND].fourier )
	scaluv(vars[V_WIND].fourier, globs->rcoslat, globs->Latitudes, globs->Fouriers*globs->NumLevelRequest);
    }

  FreeSpectral(vars);

  /* -------------------------- */
  /*  Output of fourier fields  */
  /* -------------------------- */

  if ( globs->Type == 60 )
    {
      for ( code = 0; code < MaxCodes; code++ )
	if ( vars[code].selected )
	  {
	    if ( ! vars[code].fourier ) Error("Code %d not available on fourier space!", code);

	    nlevel = zaxisInqSize(vars[code].ozaxisID);
	    for ( lindex = 0; lindex < nlevel; lindex++ )
	      {
		offset = lindex*globs->DimFC;
		streamDefRecord(globs->ostreamID, vars[code].ovarID, lindex);
		streamWriteRecord(globs->ostreamID, vars[code].fourier+offset, 0);
	      }	      
	  }

      FreeFourier(vars);
      return;
    }

  /* ----------------------- */
  /*  Output of zonal means  */
  /* ----------------------- */

  if ( globs->Type == 61 )
    {
      for ( code = 0; code < MaxCodes; code++ )
	if ( vars[code].selected )
	  {
	    if ( ! vars[code].fourier ) Error("Code %d not available on zonal mean!", code);

	    nlevel = zaxisInqSize(vars[code].ozaxisID);
	    for ( lindex = 0; lindex < nlevel; lindex++ )
	      {
		offset = lindex*globs->DimFC;
		streamDefRecord(globs->ostreamID, vars[code].ovarID, lindex);
		streamWriteRecord(globs->ostreamID, vars[code].fourier+offset, vars[code].nmiss);
	      }	      
	  }

      FreeFourier(vars);
      return;
    }

  /* ------------------------------ */
  /*  Transformation to gridpoints  */
  /* ------------------------------ */

  if ( globs->Type >= 70 )
    {
      for ( code = 0; code < MaxCodes; code++ )
	if ( vars[code].needed && vars[code].fourier )
	  {
	    fieldSize = vars[code].plev * globs->DimGP;
	    if ( vars[code].grid == NULL )
	      vars[code].grid = alloc_dp(fieldSize, FieldName(code, "grid"));

	    after_FC2GP(vars[code].fourier, vars[code].grid,
			globs->Latitudes, globs->Longitudes, vars[code].plev, globs->Fouriers);
	  }
    }

  FreeFourier(vars);

  /* -------------------------------- */
  /*  Output of pressure level grids  */
  /* -------------------------------- */

  if ( globs->Type == 70 )
    {
      for ( code = 0; code < MaxCodes; code++ )
	if ( vars[code].selected )
	  {
	    nlevel = zaxisInqSize(vars[code].ozaxisID);
	    for ( lindex = 0; lindex < nlevel; lindex++ )
	      {
		offset = lindex*globs->DimGP;
		streamDefRecord(globs->ostreamID, vars[code].ovarID, lindex);
		streamWriteRecord(globs->ostreamID, vars[code].grid+offset, vars[code].nmiss);
	      }
	  }

      FreeSamp(vars);
      FreeGrid(vars);
      return;
    }
}

void after_AnalysisAddRecord(struct Control *globs, struct Variable *vars, int code, int gridID, int zaxisID, int levelID, int nmiss)
{
  long fieldSize;
  int truncation;
  int dataSize;
  int dataOffset;
  int nlevel;
  int gridSize;
  int gridtype, leveltype;

  gridtype   = gridInqType(gridID);
  leveltype  = zaxisInqType(zaxisID);
  nlevel     = zaxisInqSize(zaxisID);
  gridSize   = gridInqSize(gridID);
  dataSize   = gridSize*nlevel;
  dataOffset = gridSize*levelID;

  vars[code].nmiss0 += nmiss;

  if ( gridtype == GRID_SPECTRAL )
    {
      vars[code].sfit = TRUE;
      vars[code].hlev = globs->NumLevelRequest;
      vars[code].plev = globs->NumLevelRequest;
      if ( nlevel > 1 && leveltype == ZAXIS_PRESSURE )
	{
	  if ( code != U_WIND && code != V_WIND )
	    {
	      fieldSize = globs->Dim3SP;
	      if (vars[code].spectral0 == NULL)
		vars[code].spectral0 = alloc_dp(fieldSize, FieldName(code,"spectral"));
	      truncation = gridInqTrunc(gridID);
	      sp2sp(globs->Field, truncation, vars[code].spectral0+levelID*globs->DimSP, globs->Truncation);
	    }
	  else
	    {
	      fieldSize = globs->Dim3SP;
	      if (vars[code].spectral0 == NULL)
		vars[code].spectral0 = alloc_dp(fieldSize, FieldName(code,"spectral"));
	      after_copy_array(vars[code].spectral0+levelID*globs->DimSP,globs->Field,globs->DimSP);
	    }
	}
      else
	Error("Only pressure level data supported for spectral data!");
    }
  else
    {
      if ( nlevel > 1 && leveltype == ZAXIS_PRESSURE )
	{
	  vars[code].sfit = TRUE;
	  fieldSize = globs->Dim3GP;
	  vars[code].hlev = globs->NumLevelRequest;
	  vars[code].plev = globs->NumLevelRequest;
	  if ( vars[code].grid0 == NULL )
	    vars[code].grid0 = alloc_dp(fieldSize, FieldName(code,"grid0"));
	  after_copy_array(vars[code].grid0+levelID*globs->DimGP, globs->Field, globs->DimGP);
	}
      else
	{
	  vars[code].sfit = FALSE;
	  fieldSize = globs->DimGP;
	  vars[code].hlev = 1;
	  vars[code].plev = 1;
	  if ( vars[code].grid0 == NULL )
	    vars[code].grid0 = alloc_dp(fieldSize, FieldName(code,"grid0"));
	  after_copy_array(vars[code].grid0, globs->Field, globs->DimGP);
	}

      if ( globs->Mean > 0 && (nmiss > 0 || vars[code].samp) )
	{
	  int i;
	  if ( vars[code].samp == NULL )
	    {
	      vars[code].samp = (int *) malloc(dataSize*sizeof(int));
	      for ( i = 0; i < dataSize; i++ ) vars[code].samp[i] = globs->MeanCount0;
	    }

	  for ( i = 0; i < gridSize; i++ )
	    if ( IS_NOT_EQUAL(globs->Field[i], vars[code].missval) ) vars[code].samp[i+dataOffset]++;
	}
    }
}


void after_EchamAddRecord(struct Control *globs, struct Variable *vars, int code, int gridID, int zaxisID, int levelID, int nmiss)
{
  int dataSize;
  int dataOffset;
  int nlevel;
  int gridSize;
  int gridtype, leveltype;

  gridtype   = gridInqType(gridID);
  leveltype  = zaxisInqType(zaxisID);
  nlevel     = zaxisInqSize(zaxisID);
  gridSize   = gridInqSize(gridID);
  dataSize   = gridSize*nlevel;
  dataOffset = gridSize*levelID;

  vars[code].nmiss0 += nmiss;

  if ( gridtype == GRID_SPECTRAL )
    {
      /* ---------------------------------------------------------- */
      /* Found spectral field ! If needed, allocate memory and copy */
      /* ---------------------------------------------------------- */
      vars[code].sfit = TRUE;
      vars[code].hlev = nlevel;
      vars[code].plev = 1;
      if ( nlevel > 1 && leveltype == ZAXIS_HYBRID )
	vars[code].plev = globs->NumLevelRequest;

      if ( gridInqTrunc(gridID) != globs->Truncation )
	Error("Resolution error. Required %d - Found %d", globs->Truncation, gridInqTrunc(gridID));

      if ( vars[code].spectral0 == NULL )
	vars[code].spectral0 = alloc_dp(dataSize, FieldName(code, "spectral0"));

      after_copy_array(vars[code].spectral0+dataOffset, globs->Field, gridSize);
    }
  else
    {
      vars[code].hlev = nlevel;
      vars[code].plev = nlevel;
      vars[code].sfit = FALSE;
      if ( nlevel > 1 && leveltype == ZAXIS_HYBRID && nlevel == globs->NumLevel )
	{
	  vars[code].plev = globs->NumLevelRequest;
	  vars[code].sfit = TRUE;
	}

      if ( globs->Mean > 0 && (nmiss > 0 || vars[code].samp) )
	{
	  int i;
	  if ( vars[code].samp == NULL )
	    {
	      vars[code].samp = (int *) malloc(dataSize*sizeof(int));
	      for ( i = 0; i < dataSize; i++ ) vars[code].samp[i] = globs->MeanCount0;
	    }

	  for ( i = 0; i < gridSize; i++ )
	    if ( IS_NOT_EQUAL(globs->Field[i], vars[code].missval) ) vars[code].samp[i+dataOffset]++;
	}

      if ( vars[code].hybrid0 == NULL )
	vars[code].hybrid0 = alloc_dp(dataSize, FieldName(code, "hybrid0"));

      after_copy_array(vars[code].hybrid0+dataOffset, globs->Field, gridSize);
    }
}

static
void MakeDependencies(struct Variable *vars, int varcode, int depcode)
{
  if ( vars[varcode].needed && ! vars[varcode].detected )
    {
      vars[depcode].needed = TRUE;
      vars[varcode].comp   = TRUE;

      if ( afterDebug )
	fprintf(stderr, "Needed code %d to compute code %d\n", depcode, varcode);

      if ( vars[depcode].ivarID == -1 )
	{
	  if ( depcode == U_WIND )
	    {
	      MakeDependencies(vars, U_WIND, DIVERGENCE);
	      MakeDependencies(vars, U_WIND, VORTICITY);
	    }
	  if ( depcode == V_WIND )
	    {
	      MakeDependencies(vars, V_WIND, DIVERGENCE);
	      MakeDependencies(vars, V_WIND, VORTICITY);
	    }
	}

      if ( vars[varcode].ivarID == -1 )
	{
	  if ( vars[depcode].ivarID == -1 )
	    {
	      Error("code %d undefined, needed to compute code %d", depcode, varcode);
	    }
	  else
	    {
	      vars[varcode].ivarID   = vars[depcode].ivarID;
	      vars[varcode].igridID  = vars[depcode].igridID;
	      vars[varcode].ogridID  = vars[depcode].ogridID;
	      vars[varcode].izaxisID = vars[depcode].izaxisID;
	      vars[varcode].ozaxisID = vars[depcode].ozaxisID;
	    }
	}
    }
}

static
void CheckDependencies(struct Variable *vars, int analysisdata)
{
  MakeDependencies(vars, VELOPOT, U_WIND);
  MakeDependencies(vars, VELOPOT, V_WIND);
  MakeDependencies(vars, VELOPOT, VORTICITY);
  MakeDependencies(vars, VELOPOT, DIVERGENCE);

  MakeDependencies(vars, STREAM, U_WIND);
  MakeDependencies(vars, STREAM, V_WIND);
  MakeDependencies(vars, STREAM, VORTICITY);
  MakeDependencies(vars, STREAM, DIVERGENCE);

  MakeDependencies(vars, VORTICITY, U_WIND);
  MakeDependencies(vars, VORTICITY, V_WIND);

  MakeDependencies(vars, DIVERGENCE, U_WIND);
  MakeDependencies(vars, DIVERGENCE, V_WIND);

  MakeDependencies(vars, U_WIND, VORTICITY);
  MakeDependencies(vars, U_WIND, DIVERGENCE);
  MakeDependencies(vars, U_WIND, V_WIND);

  MakeDependencies(vars, V_WIND, VORTICITY);
  MakeDependencies(vars, V_WIND, DIVERGENCE);
  MakeDependencies(vars, V_WIND, U_WIND);

  MakeDependencies(vars, WINDSPEED, U_WIND);
  MakeDependencies(vars, WINDSPEED, V_WIND);

  if ( analysisdata )
    {
      MakeDependencies(vars, RHUMIDITY, HUMIDITY);
      MakeDependencies(vars, RHUMIDITY, TEMPERATURE);
      MakeDependencies(vars, HUMIDITY, RHUMIDITY);
      MakeDependencies(vars, HUMIDITY, TEMPERATURE);
      MakeDependencies(vars, GEOPOTHEIGHT, GEOPOTENTIAL);
    }
  else
    {
      MakeDependencies(vars, THETAF, TEMPERATURE);
      MakeDependencies(vars, SLP, TEMPERATURE);
    }

  MakeDependencies(vars, SW_ATM, 176);
  MakeDependencies(vars, SW_ATM, 178);
  MakeDependencies(vars, LW_ATM, 177);
  MakeDependencies(vars, LW_ATM, 179);
  MakeDependencies(vars, NET_ATM, 176);
  MakeDependencies(vars, NET_ATM, 177);
  MakeDependencies(vars, NET_ATM, 178);
  MakeDependencies(vars, NET_ATM, 179);

  MakeDependencies(vars, SW_BOT_CLF, 176);
  MakeDependencies(vars, SW_BOT_CLF, 185);
  MakeDependencies(vars, LW_BOT_CLF, 177);
  MakeDependencies(vars, LW_BOT_CLF, 186);
  MakeDependencies(vars, SW_TOP_CLF, 178);
  MakeDependencies(vars, SW_TOP_CLF, 187);
  MakeDependencies(vars, LW_TOP_CLF, 179);
  MakeDependencies(vars, LW_TOP_CLF, 188);
  MakeDependencies(vars, NET_TOP_CLF, 178);
  MakeDependencies(vars, NET_TOP_CLF, 179);
  MakeDependencies(vars, NET_TOP_CLF, 187);
  MakeDependencies(vars, NET_TOP_CLF, 188);

  MakeDependencies(vars, ALL_WATER, 222);
  MakeDependencies(vars, LOW_WATER, 222);
  MakeDependencies(vars, MID_WATER, 222);
  MakeDependencies(vars, HIH_WATER, 222);

  MakeDependencies(vars, LOW_CLOUD, 223);
  MakeDependencies(vars, MID_CLOUD, 223);
  MakeDependencies(vars, HIH_CLOUD, 223);

  if ( vars[LOW_CLOUD].comp || vars[MID_CLOUD].comp || vars[HIH_CLOUD].comp )
    {
      static int zaxisID = -999;
      if ( zaxisID == -999 ) zaxisID = zaxisCreate(ZAXIS_SURFACE, 1);

      vars[LOW_CLOUD].izaxisID = zaxisID;
      vars[LOW_CLOUD].ozaxisID = zaxisID;
      vars[MID_CLOUD].izaxisID = zaxisID;
      vars[MID_CLOUD].ozaxisID = zaxisID;
      vars[HIH_CLOUD].izaxisID = zaxisID;
      vars[HIH_CLOUD].ozaxisID = zaxisID;
    }
}

void after_AnalysisDependencies(struct Variable *vars, int ncodes)
{
   int code;

   for ( code = 0; code < ncodes; code++ )
     vars[code].needed = vars[code].selected;

   MakeDependencies(vars, PS, LNPS);

   CheckDependencies(vars, 1);
}

void after_EchamDependencies(struct Variable *vars, int ncodes, int type, int source)
{
  int code;

  for ( code = 0; code < ncodes; code++ )
    vars[code].needed = vars[code].selected;

  for ( code = 0; code < ncodes; code++ )
    if ( vars[code].detected == FALSE ) vars[code].ivarID = -1;

  if ( type >= 50 )
    {
      vars[U_WIND].needed |= vars[DIVERGENCE].needed;
      vars[V_WIND].needed |= vars[DIVERGENCE].needed;
      vars[U_WIND].needed |=  vars[VORTICITY].needed;
      vars[V_WIND].needed |=  vars[VORTICITY].needed;
      vars[U_WIND].needed |=    vars[VELOPOT].needed;
      vars[V_WIND].needed |=    vars[VELOPOT].needed;
      vars[U_WIND].needed |=     vars[STREAM].needed;
      vars[V_WIND].needed |=     vars[STREAM].needed;
    }

  if ( type >= 30 ) vars[LNPS].needed = TRUE;

  if ( type >= 20 )
    {
      MakeDependencies(vars, THETAF, LNPS);
      MakeDependencies(vars, SLP, LNPS);
      MakeDependencies(vars, SLP, GEOPOTENTIAL);
      /*
      MakeDependencies(vars, SLP, HALF_PRESS);
      MakeDependencies(vars, SLP, FULL_PRESS);
      */
      MakeDependencies(vars, RHUMIDITY, TEMPERATURE);
      MakeDependencies(vars, RHUMIDITY, HUMIDITY);
      MakeDependencies(vars, RHUMIDITY, LNPS);
      MakeDependencies(vars, GEOPOTHEIGHT, TEMPERATURE);
      MakeDependencies(vars, GEOPOTHEIGHT, HUMIDITY);
      MakeDependencies(vars, GEOPOTHEIGHT, LNPS);
      MakeDependencies(vars, OMEGA, DIVERGENCE);
      MakeDependencies(vars, OMEGA, U_WIND);
      MakeDependencies(vars, OMEGA, V_WIND);
      MakeDependencies(vars, OMEGA, LNPS);
      MakeDependencies(vars, OMEGA, DPSDX);
      MakeDependencies(vars, OMEGA, DPSDY);
    }

   MakeDependencies(vars, DPSDX, LNPS);
   MakeDependencies(vars, HALF_PRESS, LNPS);
   MakeDependencies(vars, PS, LNPS);

   MakeDependencies(vars, SURF_RUNOFF, 142);
   MakeDependencies(vars, SURF_RUNOFF, 143);
   MakeDependencies(vars, SURF_RUNOFF, 182);
   MakeDependencies(vars, SURF_RUNOFF, 221);  /* snow depth change */

   MakeDependencies(vars, THETAF, TS);

   MakeDependencies(vars, FRESH_WATER, 142);
   MakeDependencies(vars, FRESH_WATER, 143);
   MakeDependencies(vars, FRESH_WATER, 182);

   MakeDependencies(vars, PRECIP, 142);
   MakeDependencies(vars, PRECIP, 143);

   if ( source != S_ECHAM5 )
     {
       MakeDependencies(vars, NET_WATER, 142);
       MakeDependencies(vars, NET_WATER, 143);
       MakeDependencies(vars, NET_WATER, 160);
       MakeDependencies(vars, NET_WATER, 182);

       MakeDependencies(vars, NET_TOP, 178);
       MakeDependencies(vars, NET_TOP, 179);

       MakeDependencies(vars, NET_BOT, 176);
       MakeDependencies(vars, NET_BOT, 177);

       MakeDependencies(vars, NET_HEAT, 146);
       MakeDependencies(vars, NET_HEAT, 147);
       MakeDependencies(vars, NET_HEAT, 176);
       MakeDependencies(vars, NET_HEAT, 177);
       MakeDependencies(vars, NET_HEAT, 218);
       /*
	 if ( source == S_ECHAM5 )
	 {
	 MakeDependencies(vars, NET_HEAT, 206);
	 MakeDependencies(vars, NET_HEAT, 208);
	 MakeDependencies(vars, NET_HEAT, 209);
	 }
	 else
       */
       {
	 MakeDependencies(vars, NET_HEAT, 220);
       }
     }

   MakeDependencies(vars, SW_CLF, 178);
   MakeDependencies(vars, SW_CLF, 224);

   MakeDependencies(vars, LW_CLF, 179);
   MakeDependencies(vars, LW_CLF, 225);

   MakeDependencies(vars, NET_CLF, 178);
   MakeDependencies(vars, NET_CLF, 179);
   MakeDependencies(vars, NET_CLF, 224);
   MakeDependencies(vars, NET_CLF, 225);

   if ( vars[DPSDX].needed || vars[DPSDY].needed || 
	vars[GEOPOTHEIGHT].comp || vars[SLP].comp || vars[THETAF].needed ||
        vars[HALF_PRESS].needed || vars[RHUMIDITY].comp || vars[OMEGA].comp ||
        type >= 30 )
     vars[PS_PROG].comp = TRUE;

   CheckDependencies(vars, 0);
}

void after_legini_full(int ntr, int nlat, double *restrict poli, double *restrict pold, double *restrict pdev,
		       double *restrict pol2, double *restrict pol3, double *restrict coslat);
void after_legini(int ntr, int nlat, double *restrict poli, double *restrict pold, double *restrict coslat);

void after_legini_setup(struct Control *globs, struct Variable *vars)
{
  int ntr   = globs->Truncation;
  int nlat  = globs->Latitudes;
  int dimsp = (ntr + 1) * (ntr + 2);
  int pdim  = dimsp / 2 * nlat;

  globs->poli = (double *) malloc(pdim*sizeof(double));
  
  if ( ! globs->AnalysisData )
    {
      if ( globs->Type >= 20 )  globs->pold = (double *) malloc(pdim*sizeof(double));
      if ( vars[DPSDY].needed ) globs->pdev = (double *) malloc(pdim*sizeof(double));
    }
  
  if ( (vars[DIVERGENCE].needed || vars[VORTICITY].needed ||
	vars[VELOPOT].needed || vars[STREAM].needed ) && globs->Type > 20 )
    {
      globs->pol2 = (double *) malloc(pdim*sizeof(double));
      globs->pol3 = (double *) malloc(pdim*sizeof(double));
    }

  if ( globs->AnalysisData && (globs->Type == 70) && globs->Gaussian && !globs->Spectral )
    {
      if ( globs->poli ) { free(globs->poli); globs->poli = NULL;}
      if ( globs->pol2 ) { free(globs->pol2); globs->pol2 = NULL;}
      if ( globs->pol3 ) { free(globs->pol3); globs->pol3 = NULL;}
      return;
    }

  int flag = 1;
  //if ( globs->poli && globs->pold && globs->pdev == NULL && globs->pol2 == NULL && globs->pol3 == NULL ) flag = 0;

  if ( flag )
    after_legini_full(ntr, nlat, globs->poli, globs->pold, globs->pdev,
		      globs->pol2, globs->pol3, globs->coslat);
  else
    after_legini(ntr, nlat, globs->poli, globs->pold, globs->coslat);

  for ( int jgl = 0; jgl < nlat; ++jgl )
    globs->rcoslat[jgl] = 1.0 / globs->coslat[jgl];
  
  for ( int jgl = 0; jgl < nlat; ++jgl )
    globs->DerivationFactor[jgl] = globs->rcoslat[jgl] / PlanetRadius;
}
