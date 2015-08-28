/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2015 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

/*
   This module contains the following operators:

      Specinfo specinfo  Spectral information
*/

#include <ctype.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "pstream.h"


#define NTR2NSP(ntr)          ((ntr+1)*(ntr+2))
#define NSP2NTR(nsp)          ((int) ((((sqrt((double)(4*nsp+1)))-3)/2)))
#define NGP2NLEVEL(ngp)        ((int) (log10(((double)ngp)/80.)/log10(4.)))
#define NGP_ICON(nrooti,nlevel) ((int) (20*nrooti*nrooti*ipow(4, nlevel)))
/*#define NGP_GME(ni)           ((ni+1)*(ni+1)*10)*/
#define NGP_GME(ni)           (2+ni*ni*10)
#define NGP2NI(ngp)           ((int) sqrt((double)ngp/10.) - 1)


static void fac(int nlonin, int *nlonout, int *ierr)
{
  int n2, n3, n5;
  int m;
  
  n2 = 0;
  n3 = 0;
  n5 = 0;

  m = nlonin;

  while (m%2 == 0)
    {
      m = m/2;
      n2++;
    }
  while (m%3 == 0)
    {
      m = m/3;
      n3++;
    }
  while (m%5 == 0)
    {
      m = m/5;
      n5++;
    }

  if (m == 1) {
    *nlonout = nlonin;
    *ierr = 0;
  } else {
    *nlonout =  nlonin+1;
    *ierr = 1;
  }

  return;
}


static int compnlon(int nlat)
{
  int nlon, n;

  nlon = 2 * nlat;

  /* check that FFT works with nlon */
  while ( 1 )
    {
      n = nlon;
      if    ( n % 8 == 0 )  { n /= 8; }
      while ( n % 6 == 0 )  { n /= 6; }
      while ( n % 5 == 0 )  { n /= 5; }
      while ( n % 4 == 0 )  { n /= 4; }
      while ( n % 3 == 0 )  { n /= 3; }
      if    ( n % 2 == 0 )  { n /= 2; }

      if ( n <= 8 ) break;

      nlon = nlon + 4;

      if ( nlon > 9999 )
	{
	  nlon = 2 * nlat;
	  fprintf(stderr, "FFT does not work with len %d!\n", nlon);
	  break;
	}
    }

  return (nlon);
}


static int nlat2nlon(int nlat)
{
  int nlon, m, ierr;

  if ( nlat == 0 )
    cdoAbort("nlat = 0!");

  nlon = 2*nlat;

  fac(nlon, &m, &ierr);
  /* adjust till fft is possible */
  while (ierr != 0) 
    {
      nlon = m;
      /* correct here nlon so that nlat keeps always even */
      while (nlon%4 != 0) nlon++;
      fac(nlon, &m, &ierr);
    }

  return (nlon);
}


int ngp2ntr(int ngp)
{
  int ntr   = (int)lround(sqrt(0.25+ngp)-1.5);
  int nlonl = compnlon(ntr2nlat_linear(ntr));
  int nlatl = nlonl/2;

  ntr = (2*nlatl-1)/2;

  return (ntr);
}


int ipow(int i1, int i2)
{
  int i;
  int i3 = 1;

  for ( i = 0; i < i2; ++i ) i3 *= i1;

  return (i3);
}


void lookup_ni(int nsp, int *nroot, int *ni)
{
  int tbl2[12], tbl3[12], tbl5[12];
  int i, d, d2, n2, d3, n3, d5, n5 ;
    
  for ( i = 0; i < 12; ++i )
    {
      tbl2[i] = 10*2*2*ipow(4, (i+1))+2;
      tbl3[i] = 10*3*3*ipow(4, (i+1))+2;
      tbl5[i] = 10*5*5*ipow(4, (i+1))+2;
    }

  for ( i = 0; i < 12; ++i ) if (tbl2[i] >= nsp) break;
  n2 = i;
  d2 = tbl2[n2]-nsp;

  for ( i = 0; i < 12; ++i ) if (tbl3[i] >= nsp) break;
  n3 = i;
  d3 = tbl3[n3]-nsp;

  for ( i = 0; i < 12; ++i ) if (tbl5[i] >= nsp) break;
  n5 = i;
  d5 = tbl5[n5]-nsp;

  d = d2;
  if ( d3 < d ) d = d3;
  if ( d5 < d ) d = d5;

  if ( d == d2 )
    {
      *nroot = 2;
      *ni = 2*ipow(2, n2+1);
    }
  else if ( d == d3 )
    {
      *nroot = 3;
      *ni = 3*ipow(2, n3+1);
    }
  else if ( d == d5 )
    {
      *nroot = 5;
      *ni = 5*ipow(2, n5+1);
    }
}


void lookup_rl(int nsp, int *nroot, int *nlevel)
{
  int tbl2[12], tbl3[12], tbl5[12];
  int i, d, d2, n2, d3, n3, d5, n5 ;
    
  for ( i = 0; i < 12; ++i )
    {
      tbl2[i] = 20*2*2*ipow(4, (i+1));
      tbl3[i] = 20*3*3*ipow(4, (i+1));
      tbl5[i] = 20*5*5*ipow(4, (i+1));
    }

  for ( i = 0; i < 12; ++i ) if (tbl2[i] >= nsp) break;
  n2 = i;
  d2 = tbl2[n2]-nsp;

  for ( i = 0; i < 12; ++i ) if (tbl3[i] >= nsp) break;
  n3 = i;
  d3 = tbl3[n3]-nsp;

  for ( i = 0; i < 12; ++i ) if (tbl5[i] >= nsp) break;
  n5 = i;
  d5 = tbl5[n5]-nsp;

  d = d2;
  if ( d3 < d ) d = d3;
  if ( d5 < d ) d = d5;

  if ( d == d2 )
    {
      *nroot = 2;
      *nlevel = n2+1;
    }
  else if ( d == d3 )
    {
      *nroot = 3;
      *nlevel = n3+1;
    }
  else if ( d == d5 )
    {
      *nroot = 5;
      *nlevel = n5+1;
    }
}


void *Specinfo(void *argument)
{
  char arg[128], *parg;
  int len, i, nout1 = 0, nout2 = 0;
  int ntr1 = 0, nsp1 = 0, nlat1 = 0, nlon1 = 0, ngp1 = 0, ni1 = 0, ngp_gme1 = 0;
  int ntr2 = 0, nsp2 = 0, nlat2 = 0, nlon2 = 0, ngp2 = 0, ni2 = 0, ngp_gme2 = 0;
  int nlevel1 = 0, nlevel2 = 0, ngp_icon1 = 0, ngp_icon2 = 0;
  int nrootg1 = 0, nrooti1 = 0, nrootg2 = 0, nrooti2 = 0;

  cdoInitialize(argument);

  operatorInputArg("Txx, TLxx, NLON=xx, NLAT=xx, NIxx or ICONRyyLxx");

  len = strlen(operatorArgv()[0]);

  if ( (len+1) >= 128 ) cdoAbort("Parameter string too large!");

  for ( i = 0; i < len; i++ ) arg[i] = toupper(operatorArgv()[0][i]);
  arg[len] = 0;

  if ( arg[0] == 'T' && arg[1] == 'L' )
    {
      parg = &arg[2];
      if ( *parg == '=' ) parg++;
      if ( ! isdigit((int) *parg) ) cdoAbort("Wrong parameter: %s", arg);
      ntr2   = atoi(parg);
      nsp2   = NTR2NSP(ntr2);
      nlat2  = ntr2nlat_linear(ntr2);
      nlon2  = compnlon(nlat2);
      ngp2   = nlon2*nlat2;

      lookup_ni(nsp2, &nrootg2, &ni2);
      lookup_rl(nsp2, &nrooti2, &nlevel2);

      nout1  = FALSE;
      nout2  = TRUE;
    }
  else if ( arg[0] == 'T' )
    {
      parg = &arg[1];
      if ( *parg == '=' ) parg++;
      if ( ! isdigit((int) *parg) ) cdoAbort("Wrong parameter: %s", arg);
      ntr1   = atoi(parg);
      nsp1   = NTR2NSP(ntr1);
      nlat1  = ntr2nlat(ntr1);
      nlon1  = compnlon(nlat1);
      ngp1   = nlon1*nlat1;

      lookup_ni(nsp1, &nrootg1, &ni1);
      lookup_rl(nsp1, &nrooti1, &nlevel1);

      nout1  = TRUE;
      nout2  = FALSE;
    }
  else if ( arg[0] == 'N' && arg[1] == 'I' )
    {
      parg = &arg[2];
      if ( *parg == '=' ) parg++;
      if ( ! isdigit((int) *parg) ) cdoAbort("Wrong parameter: %s", arg);
      ni1    = atoi(parg);
      ni2    = ni1;
      ngp_gme1 = NGP_GME(ni1);
      ngp_gme2 = NGP_GME(ni2);

      ntr1 = ngp2ntr(ngp_gme1);

      nsp1 = NTR2NSP(ntr1);
      ntr1 = NSP2NTR(nsp1);
      ntr2 = ntr1;

      nlat1 = ntr2nlat(ntr1);
      nlon1 = compnlon(nlat1);
      nlat1 = nlon1/2;

      nlat2 = ntr2nlat_linear(ntr2);
      nlon2 = compnlon(nlat2);
      nlat2 = nlon2/2;

      /* lookup_ni(nsp1, &nrootg1, &ni1); */
      lookup_rl(nsp1, &nrooti1, &nlevel1);

      nrootg2 = nrootg1;
      ni2 = ni1;
      nrooti2 = nrooti1;
      nlevel2 = nlevel1;

      nout1  = TRUE;
      nout2  = TRUE;
    }
  else if ( arg[0] == 'N' && arg[1] == 'L' && arg[2] == 'O' && arg[3] == 'N' )
    {
      parg = &arg[4];
      if ( *parg == '=' ) parg++;
      if ( ! isdigit((int) *parg) ) cdoAbort("Wrong parameter: %s", arg);
      nlon1  = atoi(parg);
      nlon2  = nlon1;
      nlat1  = nlon1 / 2;
      nlat2  = nlon2 / 2;
      nlon1  = nlat2nlon(nlat1);
      nlon2  = nlat2nlon(nlat2);
      nlat1  = nlon1 / 2;
      nlat2  = nlon2 / 2;
      ntr1   = (nlat1*2-1)/3;
      ntr2   = (nlat2*2-1)/2;
      ngp1   = nlon1*nlat1;
      ngp2   = nlon2*nlat2;

      nsp1   = NTR2NSP(ntr1);
      nsp2   = NTR2NSP(ntr2);

      lookup_ni(nsp1, &nrootg1, &ni1);
      lookup_rl(nsp1, &nrooti1, &nlevel1);

      lookup_ni(nsp2, &nrootg2, &ni2);
      lookup_rl(nsp2, &nrooti2, &nlevel2);

      nout1  = TRUE;
      nout2  = TRUE;
    }
  else if ( arg[0] == 'N' && arg[1] == 'L' && arg[2] == 'A' && arg[3] == 'T' )
    {
      parg = &arg[4];
      if ( *parg == '=' ) parg++;
      if ( ! isdigit((int) *parg) ) cdoAbort("Wrong parameter: %s", arg);
      nlat1  = atoi(parg);
      nlat2  = nlat1;
      nlon1  = nlat2nlon(nlat1);
      nlon2  = nlat2nlon(nlat2);
      nlat1  = nlon1 / 2;
      nlat2  = nlon2 / 2;
      ntr1   = (nlat1*2-1)/3;
      ntr2   = (nlat2*2-1)/2;
      ngp1   = nlon1*nlat1;
      ngp2   = nlon2*nlat2;

      nsp1   = NTR2NSP(ntr1);
      nsp2   = NTR2NSP(ntr2);

      lookup_ni(nsp1, &nrootg1, &ni1);
      lookup_rl(nsp1, &nrooti1, &nlevel1);

      lookup_ni(nsp2, &nrootg2, &ni2);
      lookup_rl(nsp2, &nrooti2, &nlevel2);

      nout1  = TRUE;
      nout2  = TRUE;
    }
  else if ( arg[0] == 'N' )
    {
      parg = &arg[1];
      if ( *parg == '=' ) parg++;
      if ( ! isdigit((int) *parg) ) cdoAbort("Wrong parameter: %s", arg);
      nlat1  = 2*atoi(parg);
      nlat2  = nlat1;
      nlon1  = nlat2nlon(nlat1);
      nlon2  = nlat2nlon(nlat2);
      nlat1  = nlon1 / 2;
      nlat2  = nlon2 / 2;
      ntr1   = (nlat1*2-1)/3;
      ntr2   = (nlat2*2-1)/2;
      ngp1   = nlon1*nlat1;
      ngp2   = nlon2*nlat2;

      nsp1   = NTR2NSP(ntr1);
      nsp2   = NTR2NSP(ntr2);

      lookup_ni(nsp1, &nrootg1, &ni1);
      lookup_rl(nsp1, &nrooti1, &nlevel1);

      lookup_ni(nsp2, &nrootg2, &ni2);
      lookup_rl(nsp2, &nrooti2, &nlevel2);

      nout1  = TRUE;
      nout2  = TRUE;
    }
  else if ( arg[0] == 'I' && arg[1] == 'C' && arg[2] == 'O' && arg[3] == 'N' )
    {
      parg = &arg[4];
      if ( *parg != 'R' ) cdoAbort("Wrong parameter: %s", arg);
      parg++;
      if ( ! isdigit((int) *parg) ) cdoAbort("Wrong parameter: %s", arg);
      nrooti1 = atoi(parg);
      nrooti2 = nrooti1;
      while ( isdigit((int) *parg) ) parg++;
      if ( *parg != 'L' ) cdoAbort("Wrong parameter: %s", arg);
      parg++;
      if ( ! isdigit((int) *parg) ) cdoAbort("Wrong parameter: %s", arg);
      nlevel1 = atoi(parg);
      nlevel2 = nlevel1;
      ngp_icon1 = NGP_ICON(nrooti1,nlevel1);
      ngp_icon2 = NGP_ICON(nrooti1,nlevel2);

      ntr1 = ngp2ntr(ngp_icon1);
      nsp1 = NTR2NSP(ntr1);
      ntr1 = NSP2NTR(nsp1);
      ntr2 = ntr1;

      nlat1 = ntr2nlat(ntr1);
      nlon1 = compnlon(nlat1);
      nlat1 = nlon1/2;

      nlat2 = ntr2nlat_linear(ntr2);
      nlon2 = compnlon(nlat2);
      nlat2 = nlon2/2;

      lookup_ni(nsp1, &nrootg1, &ni1);
      /* lookup_rl(nsp1, &nrooti1, &nlevel1);*/

      nrootg2 = nrootg1;
      ni2 = ni1;
      nrooti2 = nrooti1;
      nlevel2 = nlevel1;

      nout1  = TRUE;
      nout2  = TRUE;
    }
  else
    cdoAbort("Unsupported parameter: %s", arg);

  nsp1      = NTR2NSP(ntr1);
  nsp2      = NTR2NSP(ntr2);
  ngp1      = nlon1*nlat1;
  ngp2      = nlon2*nlat2;
  ngp_gme1  = NGP_GME(ni1);
  ngp_gme2  = NGP_GME(ni2);
  ngp_icon1 = NGP_ICON(nrooti1,nlevel1);
  ngp_icon2 = NGP_ICON(nrooti2,nlevel2);

  fprintf(stdout, "truncation     nsp  nlon  nlat      ngp  gme    ngp_gme  icon   ngp_icon\n");

  if ( nout1 ) fprintf(stdout, "   T%-4d  %8d %5d %5d %8d  ni%d %8d  R%dB%02d  %8d\n",
		       ntr1, nsp1, nlon1, nlat1, ngp1, ni1, ngp_gme1, nrooti1, nlevel1, ngp_icon1);

  if ( nout2 ) fprintf(stdout, "   TL%-4d %8d %5d %5d %8d  ni%d %8d  R%dB%02d  %8d\n",
		       ntr2, nsp2, nlon2, nlat2, ngp2, ni2, ngp_gme2, nrooti2, nlevel2, ngp_icon2);

  cdoFinish();

  return (0);
}
