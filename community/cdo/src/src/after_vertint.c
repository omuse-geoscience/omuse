#include <stdio.h>
#include <stdlib.h> /* exit */
#include <string.h>
#include <math.h>
#include <assert.h>

#include "compare.h"
#include "constants.h"
#include "after_vertint.h"

#define  SCALEHEIGHT     (-7000.)
#define  SCALESLP        (101325.0)

int Mars = 0;


void height2pressure(double * restrict phlev, const double * restrict hlev, long nphlev)
{
  long k;
  double exp_arg;
  double height;

  for ( k = 0; k < nphlev; k++ )
    {
      height  = hlev[k];
      /*
	unitsel == 1 : hlev[k] is given in meters
	unitsel == 2 : hlev[k] is given in kilometers
	height2pressure needs meters (MKSC-standard)
      */

      exp_arg = height / SCALEHEIGHT;

      phlev[k] = SCALESLP * exp(exp_arg);
    }
}


void pressure2height(double * restrict hlev, const double * restrict plev, long nphlev)
{
  long  k;

  for ( k = 0; k < nphlev; k++ )
    {
      hlev[k] = log(plev[k]/SCALESLP)*SCALEHEIGHT;
    }
}


void presh(double * restrict fullp, double * halfp, const double *restrict vct,
	   const double *restrict ps, long nhlev, long ngp)
{
  long i, lh;
  double zp, ze;
  double *halfpres = halfp;

  if ( ps == NULL )
    {
      fprintf(stderr, "ps undefined!\n");
      exit(EXIT_FAILURE);
    }

  for ( lh = 0; lh < nhlev; lh++ )
    {
      zp = vct[lh];
      ze = vct[lh+nhlev+1];

      for ( i = 0; i < ngp; i++ ) halfpres[i] = zp + ze * ps[i];

      halfpres += ngp;
    }
  memcpy(halfpres, ps, ngp*sizeof(double));

  if ( fullp )
    {
      halfpres = halfp;
      for ( i = 0; i < ngp*nhlev; i++ )
	fullp[i] = 0.5 * (halfpres[i] + halfpres[i+ngp]);
    }

} /* presh */


void genind(int *nx, const double * restrict plev, const double * restrict fullp, long ngp, long nplev, long nhlev)
{
  long  i, lp, lh;
  int *nxl;
  double pres;

  memset(nx, 0, ngp*nplev*sizeof(int));

#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(i, lh, pres, nxl)
#endif
  for ( lp = 0; lp < nplev; lp++ )
    {
      pres = plev[lp];
      nxl  = nx + lp*ngp;
      for ( lh = 0; lh < nhlev; lh++ )
	for ( i = 0; i < ngp ; i++ )
	   {
	     if ( pres > fullp[lh*ngp+i] ) nxl[i] = lh;
	   }
    }

}  /* genind */


void genindmiss(int *nx, const double * restrict plev, int ngp, int nplev, const double * restrict ps_prog, int * restrict pnmiss)
{
  long i, lp;
  int *nxl;
  double pres;

#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(i, pres, nxl)
#endif
  for ( lp = 0; lp < nplev; lp++ )
    {
      pnmiss[lp] = 0;
      pres = plev[lp];
      nxl  = nx + lp*ngp;
      for ( i = 0; i < ngp; i++ )
	{
	  if ( pres > ps_prog[i] )
	    {
	      nxl[i] = -1;
	      pnmiss[lp]++;
	    }
	}
    }

}  /* genindmiss */


void extra_P(double * restrict slp, const double * restrict halfp, const double * restrict fullp,
	     const double * restrict geop, const double * restrict temp, long ngp)
{
  double alpha, tstar, tmsl, zprt, zprtal;
  double zrg;
  double zlapse = 0.0065;
  long j;

  zrg = 1.0 / PlanetGrav;

  for ( j = 0; j < ngp; ++j )
    {
      if ( geop[j] < 0.0001 && geop[j] > -0.0001 ) slp[j] = halfp[j];
      else
	{
	  alpha = PlanetRD * zlapse * zrg;
	  tstar = (1.0 + alpha * (halfp[j]/fullp[j] - 1.0)) * temp[j];

	  if ( tstar < 255.0 ) tstar = 0.5 * (255.0 + tstar);

	  tmsl = tstar + zlapse * zrg * geop[j];
	  if ( tmsl > 290.5 && tstar > 290.5 )
	    {
	      tstar = 0.5 * (290.5 + tstar);
	      tmsl  = tstar;
	    }

	  if ( tmsl-tstar < 0.000001 && tstar-tmsl < 0.000001 )
	    alpha = 0.0;
	  else if ( geop[j] > 0.0001 || geop[j] < -0.0001 )
	    alpha = PlanetRD * (tmsl-tstar) / geop[j];

	  zprt   = geop[j] / (PlanetRD * tstar);
	  zprtal = zprt * alpha;
	  slp[j] = halfp[j] * exp(zprt*(1.0-zprtal*(0.5-zprtal/3.0)));
	}
    }

}  /* extrap */


static 
double extra_T(double pres, double halfp, double fullp, double geop, double temp)
{
  double tstar, ztsz, z1, ztmsl, zalph, peval, zhts, zalp;
  double zrg;
  double zlapse = 0.0065;

  zrg   = 1.0 / PlanetGrav;
  tstar = (1.0 + zlapse * PlanetRD * zrg * (halfp/fullp - 1.0)) * temp;
  ztsz  = tstar;
  z1    = tstar + zlapse * zrg * geop;

  if ( tstar < 255.0 ) tstar = 0.5 * (255.0 + tstar);

  ztmsl = tstar + zlapse * zrg * geop;

  if ( ztmsl > 290.5 && tstar > 290.5 )
    {
      tstar = 0.5 * (290.5 + tstar);
      ztmsl = tstar;
    }

  if ( ztmsl > 290.5 && tstar <= 290.5 ) ztmsl=290.5;

  zalph = PlanetRD*zlapse*zrg;

  if ( ztmsl-tstar < 0.000001 && tstar-ztmsl < 0.000001 ) zalph=0.0;

  if ( (ztmsl-tstar > 0.000001 || tstar-ztmsl > 0.000001 ) &&
       (geop > 0.0001 || geop < -0.0001) )
    zalph = PlanetRD*(ztmsl-tstar)/geop;

  if ( pres <= halfp )
    peval = ((halfp-pres)*temp+ (pres-fullp)*tstar)/ (halfp-fullp);
  else
    {
      ztmsl = z1;
      tstar = ztsz;
      zhts  = geop * zrg;

      if ( zhts > 2000. && z1 > 298. )
	{
	  ztmsl = 298.;
	  if ( zhts < 2500. ) ztmsl = 0.002*((2500.-zhts)*z1+(zhts-2000.)*ztmsl);
	}

      if ( (ztmsl-tstar) < 0.000001 )
	zalph = 0.;
      else if (geop > 0.0001 || geop < -0.0001)
	zalph = PlanetRD*(ztmsl-tstar)/geop;
      else
	zalph = PlanetRD*zlapse*zrg;

      zalp  = zalph*log(pres/halfp);
      peval = tstar*(1.0+zalp*(1.0+zalp*(0.5+0.16666666667*zalp)));
    }

  return peval;

}  /* extra_T */


static 
double extra_Z(double pres, double halfp, double fullp, double geop, double temp)
{
  double alpha, tstar, tmsl, zalp, zalpal;
  double zrg;
  double zlapse = 0.0065;
  double ztlim = 290.5;

  zrg   = 1.0 / PlanetGrav;
  alpha = PlanetRD * zlapse * zrg;
  tstar = (1.0 + alpha * (halfp/fullp - 1.0)) * temp;

  if ( tstar < 255.0 ) tstar = 0.5 * (255.0 + tstar);

  tmsl = tstar + zlapse * zrg * geop;

  if ( tmsl > ztlim && tstar > ztlim )
    {
      tstar = 0.5 * (ztlim + tstar);
      tmsl  = tstar;
    }

  if ( tmsl > ztlim && tstar <= ztlim ) tmsl = ztlim;

  if ( tmsl-tstar < 0.000001 && tstar-tmsl < 0.000001 )
    alpha = 0.0;
  else if ( geop > 0.0001 || geop < -0.0001 )
    alpha = PlanetRD * (tmsl-tstar) / geop;

  zalp   = log(pres/halfp);
  zalpal = zalp * alpha;

  return ((geop - PlanetRD*tstar*zalp*(1.0 + zalpal*(0.5 + zalpal/6.0)))*zrg);
}  /* extra_Z */


void interp_X(const double * restrict gt, double *pt, const double * restrict hyb_press, const int *nx,
	      const double * restrict plev, long nplev, long ngp, long nhlev, double missval)
{
  long lp, i;
  long nl, nh;

#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(i, nl, nh)
#endif
  for ( lp = 0; lp < nplev; lp++ )
    {
      double pres = plev[lp];
      const int *nxl  = nx + lp*ngp;
      double *ptl  = pt + lp*ngp;
      for ( i = 0; i < ngp; i++ )
	{
	  if ( nxl[i] == -1 )
	    ptl[i] = missval;
	  else
	    {
	      nl = nxl[i] * ngp + i;
	      nh = nl + ngp;
	      if ( nh >= ngp*nhlev )
		ptl[i] =  gt[nl];
	      else
		ptl[i] =  gt[nl] + (pres-hyb_press[nl])
		       * (gt[nh] - gt[nl])
                       / (hyb_press[nh] - hyb_press[nl]);
	    }
	}
    }
}  /* interp_X */


void interp_T(const double * restrict geop, const double * restrict gt, double *pt, const double * restrict fullp,
	      const double * restrict halfp, const int *nx, const double * restrict plev, long nplev, long ngp,
	      long nhlev, double missval)
{
  long lp, i;
  long nl, nh;

#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(i, nl, nh)
#endif
  for ( lp = 0; lp < nplev; lp++ )
    {
      double pres = plev[lp];
      const int *nxl  = nx + lp*ngp;
      double *ptl  = pt + lp*ngp;
#if defined(CRAY)
#pragma _CRI inline extra_T
#endif
      for ( i = 0; i < ngp; i++ )
	{
	  nl = nxl[i];
	  if ( nl < 0 )
	    ptl[i] = missval;
	  else
	    {
	      if ( nl > nhlev-2 )
		{
		  if ( Mars )
		    ptl[i] = gt[(nhlev-1)*ngp+i];
		  else
#if defined(SX)
#pragma cdir inline
#endif
		    ptl[i] = extra_T(pres, halfp[nhlev*ngp+i],
				     fullp[(nhlev-1)*ngp+i], geop[i],
				     gt[(nhlev-1)*ngp+i]);
		}
	      else
		{
		  nh = nl + 1;
		  ptl[i] =  gt[nl*ngp+i] + (pres-fullp[nl*ngp+i])
                         * (gt[nh*ngp+i] - gt[nl*ngp+i])
                         / (fullp[nh*ngp+i] - fullp[nl*ngp+i]);
		}
	    }
	}
    }
}  /* interp_T */


void interp_Z(const double * restrict geop, const double * restrict gz, double *pz, const double * restrict fullp,
	      const double * restrict halfp, const int *nx, const double * restrict gt, const double * restrict plev,
	      long nplev, long ngp, long nhlev, double missval)
{
  long lp, i;
  long nl, nh;

  assert(geop != NULL);
  assert(gz != NULL);
  assert(pz != NULL);
  assert(fullp != NULL);
  assert(halfp != NULL);

#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(i, nl, nh)
#endif
  for ( lp = 0; lp < nplev; lp++ )
    {
      double pres = plev[lp];
      const int *nxl  = nx + lp*ngp;
      double *pzl  = pz + lp*ngp;
#if defined(CRAY)
#pragma _CRI inline extra_Z
#endif
      for ( i = 0; i < ngp; i++ )
	{
	  nl = nxl[i];
	  if ( nl < 0 )
	    pzl[i] = missval;
	  else
	    {
	      if ( pres > halfp[(nl+1)*ngp+i] ) nl++;

	      if ( nl > nhlev-1 )
		{
		  if ( Mars )
		    pzl[i] = gt[(nhlev-1)*ngp+i];
		  else
#if defined(SX)
#pragma cdir inline
#endif
		    pzl[i] = extra_Z(pres, halfp[nhlev*ngp+i],
				     fullp[(nhlev-1)*ngp+i], geop[i],
				     gt[(nhlev-1)*ngp+i]);
		}
	      else
		{
		  nh = nl + 1;
		  pzl[i] =  gz[nl*ngp+i] + (pres-halfp[nl*ngp+i])
		         * (gz[nh*ngp+i] - gz[nl*ngp+i])
                         / (halfp[nh*ngp+i] - halfp[nl*ngp+i]);
		}
	    }
	}
    }
}  /* interp_Z */


/*
 * 3d vertical interpolation routine (see vert_interp_lev() in src/Intlevel.c)
 */
void vert_interp_lev3d(int gridsize, double missval, double *vardata1, double *vardata2,
		       int nlev2, int *lev_idx1, int *lev_idx2, double *lev_wgt1, double *lev_wgt2)
{
  int i, ilev;
  int idx1, idx2;
  int offset;
  double wgt1, wgt2;
  double w1, w2;
  double var1L1, var1L2, *var2;

  for ( ilev = 0; ilev < nlev2; ilev++ )
    {
      offset = ilev*gridsize;
      var2 = vardata2 + offset;

      for ( i = 0; i < gridsize; i++ )
	{
          idx1 = lev_idx1[offset+i];
          idx2 = lev_idx2[offset+i];
          wgt1 = lev_wgt1[offset+i];
          wgt2 = lev_wgt2[offset+i];

          /* upper/lower values from input field */
          var1L1 = *(vardata1+idx1);
          var1L2 = *(vardata1+idx2);

          /* if (cdoVerbose) printf("i:%d level %d: idx1=%d idx2=%d (offset+i:%d) wgt1=%g wgt2=%g var1L1:%g var1L2:%g ",
           *                         i,       ilev, idx1,   idx2,    offset+i,    wgt1,   wgt2,   var1L1,   var1L2);
           */
	  w1 = wgt1;
	  w2 = wgt2;
	  if ( DBL_IS_EQUAL(var1L1, missval) ) w1 = 0;
	  if ( DBL_IS_EQUAL(var1L2, missval) ) w2 = 0;

	  if ( IS_EQUAL(w1, 0) && IS_EQUAL(w2, 0) )
	    {
	      var2[i] = missval;
	    }
	  else if ( IS_EQUAL(w1, 0) )
	    {
	      if ( w2 >= 0.5 )
		var2[i] = var1L2;
	      else
		var2[i] = missval;	      
	    }
	  else if ( IS_EQUAL(w2, 0) )
	    {
	      if ( w1 >= 0.5 )
		var2[i] = var1L1;
	      else
		var2[i] = missval;	      
	    }
	  else
	    {
	      var2[i] = var1L1*w1 + var1L2*w2;
	    }
	}
    }
}

#if defined(CDO)
#include "util.h"
/*
 * Create weights for the 3d vertical coordinate
 *
 * The resulting index sets lev_idx1 and lev_idx2 contain absolute numbers,i.e.
 * wrt. the given gridsize. They can directly be used to read values from 3d
 * data fields.
 *
 * 3d version of vert_gen_weights() (src/Intlevel.c)
 */
void vert_gen_weights3d(int expol, int nlev1, int gridsize, double *lev1, int nlev2, double *lev2,
			int *lev_idx1, int *lev_idx2, double *lev_wgt1, double *lev_wgt2)
{
  int i,i1, i2;
  int    idx1 = 0, idx2 = 0;
  double val1, val2 = 0;

  for ( i = 0; i < gridsize; i++ )
    {
      for ( i2 = 0; i2 < nlev2; i2++ )
        {
          /* Because 2 levels were added to the source vertical coordinate (one on
           * top, one at the bottom), its loop starts at 1 */
          for ( i1 = 1; i1 < nlev1; i1++ )
            {
              if ( lev1[(i1-1)*gridsize+i] < lev1[i1*gridsize+i] )
                {
                  idx1 = (i1-1)*gridsize+i;
                  idx2 = i1*gridsize+i;
                }
              else
                {
                  idx1 = i1*gridsize+i;
                  idx2 = (i1-1)*gridsize+i;
                }
              val1 = lev1[idx1];
              val2 = lev1[idx2];

              if ( lev2[i2*gridsize+i] > val1 && lev2[i2*gridsize+i] <= val2 ) break;
            }

          if ( i1 == nlev1 ) 
            {
              if ( expol )
                cdoAbort("Level %g at index %d not found! Use extrapolation", lev2[i2*gridsize],i2);
              else
                cdoAbort("Level %g at index %d not found!");
            }

          if ( i1-1 == 0 ) /* destination levels ios not covert by the first two input z levels */
            {
              lev_idx1[i2*gridsize+i] = gridsize+i;
              lev_idx2[i2*gridsize+i] = gridsize+i;
              lev_wgt1[i2*gridsize+i] = 0;
              if ( expol || IS_EQUAL(lev2[i2*gridsize+i], val2) )
                lev_wgt2[i2*gridsize+i] = 1;
              else
                lev_wgt2[i2*gridsize+i] = 0;
            }
          else if ( i1 == nlev1-1 ) /* destination level is beyond the last value of the input z field */
            {
              lev_idx1[i2*gridsize+i] = (nlev1-2)*gridsize+i;
              lev_idx2[i2*gridsize+i] = (nlev1-2)*gridsize+i;
              if ( expol || IS_EQUAL(lev2[i2*gridsize+i], val2) )
                lev_wgt1[i2*gridsize+i] = 1;
              else
                lev_wgt1[i2*gridsize+i] = 0;
              lev_wgt2[i2*gridsize+i] = 0;
            }
          else /* target z values has two bounday values in input z field */
            {
              lev_idx1[i2*gridsize+i] = idx1;
              lev_idx2[i2*gridsize+i] = idx2;
              lev_wgt1[i2*gridsize+i] = (lev1[idx2]        - lev2[i2*gridsize+i]) / (lev1[idx2] - lev1[idx1]);
              lev_wgt2[i2*gridsize+i] = (lev2[i2*gridsize+i] - lev1[idx1])        / (lev1[idx2] - lev1[idx1]);

            }
  /*         if (cdoVerbose)
   *         {
   *           printf("i:%d i2:%d\ti2*gridsize+i:%d\tlev2[i2*gridsize+i]:%g\tidx1:%d\tidx2:%d\tlev1[idx1]:%g\tlev1[idx2]:%g\t",
   *                   i, i2, i2*gridsize+i,         lev2[i2*gridsize+i],    idx1,    idx2,    lev1[idx1],    lev1[idx2]);
   *           printf("\tlev_wgt1:%g\tlev_wgt2:%g\n", lev_wgt1[i2*gridsize+i], lev_wgt2[i2*gridsize+i]);
   *         }
   */
          /* backshift of the indices because of the two additional levels in input vertical coordinate */
          lev_idx1[i2*gridsize+i] -= gridsize;
          lev_idx2[i2*gridsize+i] -= gridsize;

        }
    }
}


/*
 * Create weights for the 1d vertical coordinate from a 3d vertical coordinate
 *
 * The resulting index sets lev_idx1 and lev_idx2 contain absolute numbers,i.e.
 * wrt. the given gridsize. They can directly be used to read values from 3d
 * data fields.
 *
 * 3d1d version of vert_gen_weights() (src/Intlevel.c)
 */
void vert_gen_weights3d1d(int expol, int nlev1, int gridsize, double *lev1, int nlev2, double *lev2,
			  int *lev_idx1, int *lev_idx2, double *lev_wgt1, double *lev_wgt2)
{
  int i,i1, i2;
  int    idx1 = 0, idx2 = 0;
  double val1, val2 = 0;

  for ( i = 0; i < gridsize; i++ )
    {
      for ( i2 = 0; i2 < nlev2; i2++ )
        {
          /* Because 2 levels were added to the source vertical coordinate (one on
           * top, one at the bottom), its loop starts at 1 */
          for ( i1 = 1; i1 < nlev1; i1++ )
            {
              if ( lev1[(i1-1)*gridsize+i] < lev1[i1*gridsize+i] )
                {
                  idx1 = (i1-1)*gridsize+i;
                  idx2 = i1*gridsize+i;
                }
              else
                {
                  idx1 = i1*gridsize+i;
                  idx2 = (i1-1)*gridsize+i;
                }
              val1 = lev1[idx1];
              val2 = lev1[idx2];

              if ( lev2[i2] > val1 && lev2[i2] <= val2 ) break;
            }

          if ( i1 == nlev1 ) 
            {
              if ( expol )
                cdoAbort("Level %g at index %d not found! Use extrapolation", lev2[i2],i2);
              else
                cdoAbort("Level %g at index %d not found!");
            }

          if ( i1-1 == 0 ) /* destination levels ios not covert by the first two input z levels */
            {
              lev_idx1[i2*gridsize+i] = gridsize+i;
              lev_idx2[i2*gridsize+i] = gridsize+i;
              lev_wgt1[i2*gridsize+i] = 0;
              if ( expol || IS_EQUAL(lev2[i2], val2) )
                lev_wgt2[i2*gridsize+i] = 1;
              else
                lev_wgt2[i2*gridsize+i] = 0;
            }
          else if ( i1 == nlev1-1 ) /* destination level is beyond the last value of the input z field */
            {
              lev_idx1[i2*gridsize+i] = (nlev1-2)*gridsize+i;
              lev_idx2[i2*gridsize+i] = (nlev1-2)*gridsize+i;
              if ( expol || IS_EQUAL(lev2[i2], val2) )
                lev_wgt1[i2*gridsize+i] = 1;
              else
                lev_wgt1[i2*gridsize+i] = 0;
              lev_wgt2[i2*gridsize+i] = 0;
            }
          else /* target z values has two bounday values in input z field */
            {
              lev_idx1[i2*gridsize+i] = idx1;
              lev_idx2[i2*gridsize+i] = idx2;
              lev_wgt1[i2*gridsize+i] = (lev1[idx2]  - lev2[i2]) / (lev1[idx2] - lev1[idx1]);
              lev_wgt2[i2*gridsize+i] = (lev2[i2] - lev1[idx1])  / (lev1[idx2] - lev1[idx1]);

            }
  /*         if (cdoVerbose)
   *         {
   *           printf("i:%d i2:%d\ti2*gridsize+i:%d\tlev2[i2]:%g\tidx1:%d\tidx2:%d\tlev1[idx1]:%g\tlev1[idx2]:%g\t",
   *                   i, i2, i2*gridsize+i,         lev2[i2],    idx1,    idx2,    lev1[idx1],    lev1[idx2]);
   *           printf("\tlev_wgt1:%g\tlev_wgt2:%g\n", lev_wgt1[i2*gridsize+i], lev_wgt2[i2*gridsize+i]);
   *         }
   */
          /* backshift of the indices because of the two additional levels in input vertical coordinate */
          lev_idx1[i2*gridsize+i] -= gridsize;
          lev_idx2[i2*gridsize+i] -= gridsize;

        }
    }
}
#endif
