#include <math.h>

#include "constants.h"


#define  SQUARE_RADIUS   (-PlanetRadius * PlanetRadius)

void dv2ps(const double *restrict div, double *restrict pot, long nlev, long ntr)
{
  long l, m, n;
  double fact;

  for ( l = 0; l <  nlev; l++ )
    {
      /* m == 0 */
      *pot++ = 0.0;
      *pot++ = 0.0;
      div += 2;

      for ( n = 1; n <= ntr; n++ )
	{
	  fact = SQUARE_RADIUS / (n*n + n);
	  *pot++ = *div++ * fact;
	  *pot++ = *div++ * fact;
	}

      /* m >= 0 */
      for ( m = 1; m <= ntr; m++ )
	for ( n = m; n <= ntr; n++ )
	  {
	    fact = SQUARE_RADIUS / (n*n + n);
	    *pot++ = *div++ * fact;
	    *pot++ = *div++ * fact;
	  }
      }
}


void dv2uv(double *d, double *o, double *u, double *v, double *f, double *g,
	   int nt, int nsp, int nlev)
{
  /* d(nsp,nlev), o(nsp,nlev)     ! divergence, vorticity        */
  /* u(nsp,nlev), v(nsp,nlev)     ! zonal wind, meridional wind  */
  /* f(nsp/2)   , g(nsp/2)        ! factor tables                */

  int l, m, n;
  int i;

  for ( l = 0; l < nlev; l++ )
    {
      i = 0;

      for ( m = 0; m < nt-1; m++ )
	{
	  /*********/
	  /* n = m */
	  /*********/

	  if ( m == 0 )
	    {
	      *u++ = -g[i+1] * o[2*(i+1)  ];
	      *u++ = -g[i+1] * o[2*(i+1)+1];
	      *v++ =  g[i+1] * d[2*(i+1)  ];
	      *v++ =  g[i+1] * d[2*(i+1)+1];
	    }
	  else
	    {
	      *u++ = -f[i] * d[2*i+1] - g[i+1] * o[2*(i+1)  ];
	      *u++ =  f[i] * d[2*i  ] - g[i+1] * o[2*(i+1)+1];
	      *v++ = -f[i] * o[2*i+1] + g[i+1] * d[2*(i+1)  ];
	      *v++ =  f[i] * o[2*i  ] + g[i+1] * d[2*(i+1)+1];
	    }
	  ++i;

	  /****************/
	  /* m < n < nt-1 */
	  /****************/

	  for ( n = m+1; n < nt-1; n++ )
	    {
	      *u++ =  g[i] * o[2*(i-1)  ] - f[i] * d[2*i+1] - g[i+1] * o[2*(i+1)  ];
	      *u++ =  g[i] * o[2*(i-1)+1] + f[i] * d[2*i  ] - g[i+1] * o[2*(i+1)+1];
	      *v++ = -g[i] * d[2*(i-1)  ] - f[i] * o[2*i+1] + g[i+1] * d[2*(i+1)  ];
	      *v++ = -g[i] * d[2*(i-1)+1] + f[i] * o[2*i  ] + g[i+1] * d[2*(i+1)+1];
	      ++i;
	    }

	  /************/
	  /* n = nt-1 */
	  /************/

	  *u++ =  g[i] * o[2*(i-1)  ] - f[i] * d[2*i+1];
	  *u++ =  g[i] * o[2*(i-1)+1] + f[i] * d[2*i  ];
	  *v++ = -g[i] * d[2*(i-1)  ] - f[i] * o[2*i+1];
	  *v++ = -g[i] * d[2*(i-1)+1] + f[i] * o[2*i  ];
	  ++i;

	  /**********/
	  /* n = nt */
	  /**********/

	  *u++ =  g[i] * o[2*(i-1)  ];
	  *u++ =  g[i] * o[2*(i-1)+1];
	  *v++ = -g[i] * d[2*(i-1)  ];
	  *v++ = -g[i] * d[2*(i-1)+1];
	  ++i;
	}

      /***************************/
      /* m = nt-1  and  n = nt-1 */
      /***************************/

      *u++ = -f[i] * d[2*i+1];
      *u++ =  f[i] * d[2*i  ];
      *v++ = -f[i] * o[2*i+1];
      *v++ =  f[i] * o[2*i  ];
      ++i;

      /*************************/
      /* m = nt-1  and  n = nt */
      /*************************/

      *u++ =  g[i] * o[2*(i-1)  ];
      *u++ =  g[i] * o[2*(i-1)+1];
      *v++ = -g[i] * d[2*(i-1)  ];
      *v++ = -g[i] * d[2*(i-1)+1];
      ++i;

      /***********************/
      /* m = nt  and  n = nt */
      /***********************/

      *u++ = 0.0;
      *u++ = 0.0;
      *v++ = 0.0;
      *v++ = 0.0;

      d += nsp;
      o += nsp;
    }
}

/*
void scaluv(double *fu, double rclat[], int nlat, int lot)
{
  int l, lat;
  double *ful;

#if defined (_OPENMP)
#pragma tomp parallel for default(shared) private(lat, ful)
#endif
  for ( l = 0; l < lot; l++ )
    {
      ful = fu + l*nlat;
      for ( lat = 0; lat < nlat; lat++ )
	{
	  ful[lat] = rclat[lat];
	}
    }
}
*/


void scaluv(double *fu, double *rclat, int nlat, int lot)
{
  int l,lat;

  for (l = 0; l < lot; l++)
    for (lat = 0; lat < nlat; lat++)
      {
        *fu *= rclat[lat];
        fu++;
      }
}


void uv2dv(double *fu, double *fv, double *sd, double *sv,
           double *pol2, double *pol3, int klev, int nlat, int nt)
{
  int lev, jmm, jfc, lat, nfc, nsp2;
  double dir, dii, vor, voi;
  double *ufr, *ufi, *vfr, *vfi;
  double *ful, *fvl, *sdl, *svl;
  double *po2, *po3;

  nsp2 = (nt+1)*(nt+2);
  nfc  = (nt+1)*2;

#if defined (_OPENMP)
#pragma omp parallel for default(shared) private(jmm, jfc, lat, po2, po3, ful, fvl, sdl, svl, ufr, ufi, vfr, vfi, dir, dii, vor, voi)
#endif
  for ( lev = 0; lev < klev; lev++ )
    {
      po2 = pol2;
      po3 = pol3;
      ful = fu + lev*nfc*nlat;
      fvl = fv + lev*nfc*nlat;
      sdl = sd + lev*nsp2;
      svl = sv + lev*nsp2;
      for ( jmm = 0; jmm <= nt; jmm++ )
	{
	  for  ( jfc = jmm; jfc <= nt; jfc++ )
	    {
	      ufr = ful;
	      ufi = ful + nlat;
	      vfr = fvl;
	      vfi = fvl + nlat;
	      dir = 0.0;
	      dii = 0.0;
	      vor = 0.0;
	      voi = 0.0;
	      for ( lat = 0; lat < nlat; lat++ )
		{
		  dir += vfr[lat] * po2[lat] - ufi[lat] * po3[lat];
		  dii += vfi[lat] * po2[lat] + ufr[lat] * po3[lat];
		  vor -= ufr[lat] * po2[lat] + vfi[lat] * po3[lat];
		  voi -= ufi[lat] * po2[lat] - vfr[lat] * po3[lat];
		}
	      *sdl++ = dir;
	      *sdl++ = dii;
	      *svl++ = vor;
	      *svl++ = voi;
	      po2 += nlat;
	      po3 += nlat;
	    }
	  ful += 2 * nlat;
	  fvl += 2 * nlat;
	}
    }
}


void geninx(long ntr, double *f, double *g)
{
  long m2,n2;
  long m, n ;

  for ( m = 0; m <= ntr; m++ )
    {
      m2 = m * m;
      for ( n = m; n <= ntr; n++ )
	{
	  n2 = n * n;
	  if ( n )
	    {
	      *g++ = -PlanetRadius / n * sqrt((double)(n2-m2)/(double)(4*n2-1));
	      *f++ = -PlanetRadius * m / (double)(n2+n);
	    }
	  else
	    {
	      *g++ = 0.0;
	      *f++ = 0.0;
	    }
	}
    }
}
