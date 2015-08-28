#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "constants.h"

#define  OPENMP4  201307
#if defined(_OPENMP) && defined(OPENMP4) && _OPENMP >= OPENMP4
#define  HAVE_OPENMP4  1
#endif

void gaussaw(double *pa, double *pw, int nlat);

static
void jspleg1(double *pleg, double plat, int ktrunc, double *work)
{
  /*
     jspleg1 - Routine to calculate legendre functions

     Purpose
     --------

     This routine calculates the legendre functions for one latitude.
     (but not their derivatives)


     Interface
     ----------

     jspleg1( pleg, plat, ktrunc)


     Input parameters
     ----------------

     plat      - Latitude in radians
     ktrunc    - Spectral truncation


     Output parameters
     -----------------

     pleg      - Array of legendre functions for one latitude.
                 The array must be at least (KTRUNC+1)*(KTRUNC+4)/2 
                 words long.

     Method
     ------

     Recurrence relation with explicit relations for P(m,m) and 
     P(m,m+1)


     AUTHOR
     ------

     J.D.Chambers         ECMWF        9 November 1993


     Modifications
     -------------

     None

  */
  int itout1, i1m, ilm, jm, jcn, im2;
  double zsin, zcos, zf1m, zre1, zf2m, znsqr, ze1, ze2;
  double zjmsqr;
  double *zhlp1, *zhlp2, *zhlp3;


  /* Initialization */

  itout1 = ktrunc+1;
  /*  zsin   = sin(plat); */
  zsin   = plat;
  zcos   = sqrt(1.-zsin*zsin);

  zhlp1 = work;
  zhlp2 = work + itout1;
  zhlp3 = work + itout1 + itout1;

  /*  Step 1.        M = 0, N = 0 and N = 1 */

  ilm     = 1;
  pleg[0] = 1.0;
  zf1m    = sqrt(3.0);
  pleg[1] = zf1m*zsin;

  /*  Step 2.       Sum for M = 0 to T (T = truncation) */

  for ( jm = 1; jm < itout1; jm++ )
    {
      zhlp1[jm] = sqrt(2.*jm+3.);
      zhlp2[jm] = 1./sqrt(2.*jm);
    }

  zhlp1[0] = sqrt(3.);

  for ( jm = 0; jm < itout1; jm++ )
    {
      i1m  = jm - 1;
      zre1 = zhlp1[jm];
      ze1  = 1./zre1;

      /*   Step 3.       M > 0 only */

      if ( i1m != -1 )
	{
          zf2m = zf1m*zcos*zhlp2[jm];
          zf1m = zf2m*zre1;

	  /*  Step 4.       N = M and N = M+1 */

          ilm       = ilm+1;
          pleg[ilm] = zf2m;
          ilm       = ilm+1;
          pleg[ilm] = zf1m*zsin;

	  /* When output truncation is reached, return to calling program */

          if ( jm == (itout1-1) ) break;
	}

      /*  Step 5.       Sum for N = M+2 to T+1 */

      zjmsqr = jm*jm;
      im2 = i1m+2;

      for ( jcn = im2; jcn < itout1; jcn++ )
	{
          znsqr      = (jcn + 1)*(jcn + 1);
	  zhlp3[jcn] = sqrt((4.*znsqr-1.)/(znsqr-zjmsqr));
	}

      for ( jcn = im2; jcn < itout1; jcn++ )
	{
          ze2        = zhlp3[jcn];
          ilm        = ilm+1;
          pleg[ilm]  = ze2*(zsin*pleg[ilm-1]-ze1*pleg[ilm-2]);
          ze1        = 1./ze2;
	}
    }
}


/* ============================================= */
/* phcs - Compute values of Legendre polynomials */
/*        and their meridional derivatives       */
/* ============================================= */
static
void phcs(double *pnm, double *hnm, int waves, double pmu,
	  double *ztemp1, double *ztemp2)
{
  int twowaves;

  int jk, jn, jm;

  double jnmjk;
  double zcos2;
  double lat;
  double zan;
  double zsinpar;
  double zcospar;
  double zsqp;
  double zcosfak;
  double zsinfak;
  double zq;
  double zwm2;
  double zw;
  double zwq;
  double zq2m1;
  double zwm2q2;
  double z2q2;
  double zcnm;
  double zdnm;
  double zenm;

  twowaves  = waves << 1;

  zcos2     = sqrt(1.0 - pmu*pmu);
  lat       = acos(pmu);
  zan       = 1.0;

  ztemp1[0] = 0.5;

  for ( jn = 1; jn < twowaves; jn++ )
    {
      zsqp    = 1.0 / sqrt((double)(jn + jn*jn));
      zan    *= sqrt(1.0 - 1.0/(4*jn*jn));

      zcospar = cos(lat * jn);
      zsinpar = sin(lat * jn) * jn * zsqp;
      zcosfak = 1.0;

      for ( jk = 2; jk < jn; jk += 2 )
	{
	  jnmjk = jn - jk;
	  zcosfak *= (jk-1.0) * (jn+jnmjk+2.0) / (jk * (jn+jnmjk+1.0));
	  zsinfak  = zcosfak * (jnmjk) * zsqp;
	  zcospar += zcosfak * cos(lat * jnmjk);
	  zsinpar += zsinfak * sin(lat * jnmjk);
	}

      /*  code for jk == jn */

      if ((jn & 1) == 0)
	{
	  zcosfak *= (double)((jn-1) * (jn+2)) / (double)(jn * (jn+1));
	  zcospar += zcosfak * 0.5;
	}
      ztemp1[jn  ] = zan * zcospar;
      ztemp2[jn-1] = zan * zsinpar;
    }

  memcpy(pnm, ztemp1, waves*sizeof(double));
  pnm += waves;
  memcpy(pnm, ztemp2, waves*sizeof(double));
  pnm += waves;

  hnm[0] = 0.0;
  for (jn = 1; jn < waves; jn++)
    hnm[jn] = jn * (pmu * ztemp1[jn] - sqrt((jn+jn+1.0) / (jn+jn-1.0)) * ztemp1[jn-1]);

  hnm += waves;

  hnm[0] = pmu * ztemp2[0];

  for (jn = 1; jn < waves; jn++)
    hnm[jn] = (jn+1) * pmu * ztemp2[jn]
            - sqrt(((jn+jn+3.0)*((jn+1)*(jn+1)-1.0)) / (jn+jn+1.0)) * ztemp2[jn-1];
	    
  hnm += waves;

  for (jm = 2; jm < waves; jm++)
    {
      pnm[0] = sqrt(1.0 + 1.0 / (jm+jm)) * zcos2 * ztemp2[0];
      hnm[0] = jm * pmu * pnm[0];
#if defined(CRAY)
#pragma _CRI novector
#endif
#if defined(__uxp__)
#pragma loop scalar
#endif
      for (jn = 1; jn < twowaves-jm; jn++)
	{
          zq      = jm + jm + jn - 1;
          zwm2    = zq + jn;
          zw      = zwm2 + 2;
          zwq     = zw*zq;
          zq2m1   = zq*zq - 1.;
          zwm2q2  = zwm2*zq2m1;
          z2q2    = zq2m1*2;
          zcnm    = sqrt((zwq*(zq-2.))/(zwm2q2-z2q2));
          zdnm    = sqrt((zwq*(jn+1.))/zwm2q2);
          zenm    = sqrt(zw * jn /((zq+1.0) * zwm2));
          pnm[jn] = zcnm * ztemp1[jn] - pmu
                  * (zdnm * ztemp1[jn+1] - zenm * pnm[jn-1]);
          hnm[jn] = (jm + jn) * pmu * pnm[jn]
                  - sqrt(zw * jn * (zq+1) / zwm2) * pnm[jn-1];
	}
      memcpy(ztemp1, ztemp2, twowaves*sizeof(double));
      memcpy(ztemp2, pnm   , twowaves*sizeof(double));
      pnm += waves;
      hnm += waves;
    }
}


void after_legini_full(int ntr, int nlat, double *restrict poli, double *restrict pold, double *restrict pdev,
		       double *restrict pol2, double *restrict pol3, double *restrict coslat)
{
  int jgl, jm, jn;
  int jsp;
  double gmusq;

  int waves =  ntr + 1;
  int dimsp = (ntr + 1) * (ntr + 2);

  double *gmu = (double*) malloc(nlat * sizeof(double));
  double *gwt = (double*) malloc(nlat * sizeof(double));

  gaussaw(gmu, gwt, nlat);

#if ! defined(_OPENMP)
  double *pnm    = (double*) malloc(dimsp * sizeof(double));
  double *hnm    = (double*) malloc(dimsp * sizeof(double));
  double *ztemp1 = (double*) malloc((waves<<1) * sizeof(double));
  double *ztemp2 = (double*) malloc((waves<<1) * sizeof(double));
#endif

#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(jm, jn, jsp, gmusq)
#endif
  for ( jgl = 0; jgl < nlat; jgl++ )
    {
#if defined(_OPENMP)
      double *pnm    = (double*) malloc(dimsp * sizeof(double));
      double *hnm    = (double*) malloc(dimsp * sizeof(double));
      double *ztemp1 = (double*) malloc((waves<<1) * sizeof(double));
      double *ztemp2 = (double*) malloc((waves<<1) * sizeof(double));
#endif
      gmusq = 1.0 - gmu[jgl]*gmu[jgl];
      coslat[jgl] = sqrt(gmusq);

      phcs(pnm, hnm, waves, gmu[jgl], ztemp1, ztemp2);

      double zgwt = gwt[jgl];
      double zrafgmusqr = 1./(PlanetRadius * gmusq);
      double zradsqrtgmusqr = 1./(-PlanetRadius * sqrt(gmusq));

      const int lpold = pold != NULL;
      const int lpdev = pdev != NULL;
      const int lpol2 = pol2 != NULL;
      const int lpol3 = pol3 != NULL;

      jsp = jgl;
      for ( jm = 0; jm < waves; jm++ )
	for ( jn = 0; jn < waves - jm; jn++ )
	  {
	               poli[jsp] = pnm[jm*waves+jn] * 2.0;
	    if (lpold) pold[jsp] = pnm[jm*waves+jn] * zgwt;
	    if (lpdev) pdev[jsp] = hnm[jm*waves+jn] * 2.0 * zradsqrtgmusqr;
	    if (lpol2) pol2[jsp] = hnm[jm*waves+jn] * zgwt * zrafgmusqr;
	    if (lpol3) pol3[jsp] = pnm[jm*waves+jn] * zgwt * jm * zrafgmusqr;
	    jsp += nlat;
	  }
      
#if defined(_OPENMP)
      free(ztemp2);
      free(ztemp1);
      free(pnm);
      free(hnm);
#endif
    }

#if ! defined(_OPENMP)
  free(ztemp2);
  free(ztemp1);
  free(pnm);
  free(hnm);
#endif
  free(gwt);
  free(gmu);
}


void after_legini(int ntr, int nlat, double *restrict poli, double *restrict pold, double *restrict coslat)
{
  int jgl, jm, jn, is;
  int isp, latn, lats;

  int waves  =  ntr + 1;
  int dimpnm = (ntr + 1)*(ntr + 4)/2;

  double *gmu  = (double*) malloc(nlat * sizeof(double));
  double *gwt  = (double*) malloc(nlat * sizeof(double));
  double *pnm  = (double*) malloc(dimpnm * sizeof(double));
  double *work = (double*) malloc(3*waves * sizeof(double));

  gaussaw(gmu, gwt, nlat);
  for ( jgl = 0; jgl < nlat; ++jgl ) gwt[jgl] *= 0.5;

  for ( jgl = 0; jgl < nlat; ++jgl )
    coslat[jgl] = sqrt(1.0 - gmu[jgl]*gmu[jgl]);

  for ( jgl = 0; jgl < nlat/2; jgl++ )
    {
      double zgwt = gwt[jgl];

      jspleg1(pnm, gmu[jgl], ntr, work);

      latn = jgl;
      isp = 0;
      for ( jm = 0; jm < waves; jm++ )
	{
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
	  for ( jn = 0; jn < waves - jm; jn++ )
	    {
	      is = (jn+1)%2 * 2 - 1;
	      lats = latn - jgl + nlat - jgl - 1;
	      poli[latn] = pnm[isp];
	      pold[latn] = pnm[isp] * zgwt;
	      poli[lats] = pnm[isp] * is;
	      pold[lats] = pnm[isp] * zgwt * is;
	      latn += nlat;
	      isp++;
	    }
	  isp++;
	}
    }

  free(work);
  free(pnm);
  free(gwt);
  free(gmu);
}


/* to slow for nec, 2.0 instead of 2.3 GFlops ( vector length too small ) */
void sp2fctest(const double *sa, double *fa, const double *poli, int nlev, int nlat, int nfc, int nt)
{
  int lev, jm, jn, latn, lats, is;
  double sar, sai;
  double saris, saiis;
  double pval;
  double *restrict far, *restrict fai;

  long nsp2 = (nt+1)*(nt+2);

  for ( lev = 0; lev < nlev; lev++ )
    {
      const double *restrict pol = poli;
      const double *restrict sal = sa + lev*nsp2;
      double *fal = fa + lev*nfc*nlat;
      memset(fal, 0, nfc*nlat*sizeof(double));

      for ( jm = 0; jm <= nt; jm++ )
	{
	  for ( jn = 0; jn <= nt - jm; jn++ )
	    {
	      is = (jn+1)%2 * 2 - 1;
	      sar = *sal++;
	      sai = *sal++;
              saris = sar*is;
              saiis = sai*is;
	      far = fal;
	      fai = fal + nlat;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
	      for ( latn = 0; latn < nlat/2; latn++ )
		{
		  lats = nlat - latn - 1;
		  pval = pol[latn];
		  far[latn] += pval * sar;
		  fai[latn] += pval * sai;
		  far[lats] += pval * saris;
		  fai[lats] += pval * saiis;
		}
	      pol += nlat;
	    }
	  fal += 2 * nlat;
	}
    }
}


void sp2fc(const double *sa, double *fa, const double *poli, long nlev, long nlat, long nfc, long nt)
{
  long nsp2 = (nt+1)*(nt+2);

#if defined(_OPENMP)
#pragma omp parallel for default(shared)
#endif
  for ( long lev = 0; lev < nlev; lev++ )
    {
      const double *restrict pol = poli;
      const double *restrict sal = sa + lev*nsp2;
      double *fal = fa + lev*nfc*nlat;
      memset(fal, 0, nfc*nlat*sizeof(double));

      double *restrict far, *restrict fai;
      double sar, sai;
      long jmm, jfc, lat;

      for ( jmm = 0; jmm <= nt; jmm++ )
	{
	  for ( jfc = jmm; jfc <= nt; jfc++ )
	    {
	      sar = *sal++;
	      sai = *sal++;
	      far = fal;
	      fai = fal + nlat;
              /* unaligned loop start
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
              */
	      for ( lat = 0; lat < nlat; lat++ )
		{
		  far[lat] += pol[lat] * sar;
		  fai[lat] += pol[lat] * sai;
		}
	      pol += nlat;
	    }
	  fal += 2 * nlat;
	}
    }
}


void fc2sp(double *fa, double *sa, double *poli, int nlev, int nlat, int nfc, int nt)
{
  int nsp2 = (nt+1)*(nt+2);

#if defined(_OPENMP)
#pragma omp parallel for default(shared)
#endif
  for ( int lev = 0; lev < nlev; lev++ )
    {
      const double *restrict pol = poli;
      double *fal = fa + lev*nfc*nlat;
      double *sal = sa + lev*nsp2;

      const double *restrict far;
      const double *restrict fai;
      double sar, sai;
      int jmm, jfc, lat;

      for ( jmm = 0; jmm <= nt; jmm++ )
	{
	  for ( jfc = jmm; jfc <= nt; jfc++ )
	    {
	      far = fal;
	      fai = fal + nlat;
	      sar = 0.0;
	      sai = 0.0;
#if defined(HAVE_OPENMP4)
#pragma omp simd reduction(+:sar) reduction(+:sai)
#endif
	      for ( lat = 0; lat < nlat; lat++ )
		{
		  sar += pol[lat] * far[lat];
		  sai += pol[lat] * fai[lat];
		}
	      *sal++ = sar;
	      *sal++ = sai;
	      pol += nlat;
	    }
	  fal += 2 * nlat;
	}
    }
}

/* ======================================== */
/* Convert Spectral Array to new truncation */
/* ======================================== */

void sp2sp(double *arrayIn, int truncIn, double *arrayOut, int truncOut)
{
  int n, m;

  if ( truncOut <= truncIn )
    {
      for ( n = 0; n <= truncOut; n++ )
	{
	  for ( m = n; m <= truncOut; m++ )
	    {
	      *arrayOut++ = *arrayIn++ ;
	      *arrayOut++ = *arrayIn++ ;
	    }
	  arrayIn += 2 * (truncIn-truncOut);
	}
    }
  else
    {
      for ( n = 0; n <= truncIn; n++ )
	{
	  for ( m = n; m <= truncIn; m++ )
	    {
	      *arrayOut++ = *arrayIn++ ;
	      *arrayOut++ = *arrayIn++ ;
	    }
	  for ( m = truncIn+1; m <= truncOut; ++m )
	    {
	      *arrayOut++ = 0.0;
	      *arrayOut++ = 0.0;
	    }
	}
      for ( n = truncIn+1; n <= truncOut; ++n )
	for ( m = n; m <= truncOut; ++m )
	  {
	    *arrayOut++ = 0.0;
	    *arrayOut++ = 0.0;
	  }
    }
}

/* ======================================== */
/* Cut spectral wave numbers                */
/* ======================================== */

void spcut(double *arrayIn, double *arrayOut, int trunc, int *waves)
{
  int n, m;

  for ( n = 0; n <= trunc; n++ )
    {
      for ( m = n; m <= trunc; m++ )
	{
	  if ( waves[m] )
	    {
	      *arrayOut++ = *arrayIn++;
	      *arrayOut++ = *arrayIn++;
	    }
	  else
	    {
	      *arrayOut++ = 0.0;
	      *arrayOut++ = 0.0;
	      arrayIn++;
	      arrayIn++;
	    }
	}
    }
}
