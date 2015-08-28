#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef _DMEMORY_H
#  include "dmemory.h"
#endif

#define  OPENMP4  201307
#if defined(_OPENMP) && defined(OPENMP4) && _OPENMP >= OPENMP4
#define  HAVE_OPENMP4  1
#endif


#if defined(SX)
#  define  NFFT  1024
#else
#  define  NFFT  64
#endif

#ifndef  M_SQRT2
#define  M_SQRT2     1.41421356237309504880
#endif

#define  QUA  0.25
#define  QT5  0.559016994374947

#define  S36  0.587785252292473
#define  S60  0.866025403784437
#define  S72  0.951056516295154

#define  SQ2  0.707106781186547524401

#define  D60  (S60+S60)


long get_nfft(void)
{
  return ((long) NFFT);
}


void fft_set(double *trigs, long *ifax, long n)
{
  long j, k, nfax, len = n;
  long nhl;
  double del, angle;

  del = 4.0*asin(1.0) / n;
  nhl = n / 2;
  for ( k = 0; k < nhl; k++ )
    {
      angle = k * del;
      trigs[2*k  ] = cos(angle);
      trigs[2*k+1] = sin(angle);
    }  

  nfax = 0;
  for (k = 0; k < 9; ++k) ifax[k] = 0;

  ifax[9] = n;

  if    (n % 8 == 0)  { ifax[++nfax] = 8; n /= 8; }
  while (n % 6 == 0)  { ifax[++nfax] = 6; n /= 6; }
  while (n % 5 == 0)  { ifax[++nfax] = 5; n /= 5; }
  while (n % 4 == 0)  { ifax[++nfax] = 4; n /= 4; }
  while (n % 3 == 0)  { ifax[++nfax] = 3; n /= 3; }
  if    (n % 2 == 0)  { ifax[++nfax] = 2; n /= 2; }

  ifax[0] = nfax;

#if defined(CRAY)
#pragma _CRI novector
#endif
#if defined(SX)
#pragma vdir novector
#endif
#if defined(__uxp__)
#pragma loop scalar
#endif
  for ( k = 0; k < nfax / 2; k++ )
    {
      j = ifax[k + 1];
      ifax[k + 1] = ifax[nfax - k];
      ifax[nfax - k] = j;
    }

  if ( n > 8 )
    {
      fprintf(stderr, "fft does not work with len %ld\n", len);
      exit(1);
    }
}

static
int rpassc(double *a, double *b, double *c, double *d, double *trigs,
	   long inc1, long inc2, long inc3, long inc4,
	   long lot, long n, long ifac, long la)
{
  /*
     rpassc' - performs one pass through data as part;
     of multiple real fft (fourier synthesis) routine;

     a is first real input vector
     b is equivalent to a + la * inc1
     c is first real output vector
     d is equivalent to c + ifac * la * inc2
     trigs  is a precalculated list of sines & cosines
     inc1 is the addressing increment for a;
     inc2 is the addressing increment for c;
     inc3 is the increment between input vectors a;
     inc4 is the increment between output vectors c;
     lot is the number of vectors;
     n is the length of the vectors;
     ifac is the current factor of n;
     la is the product of previous factors;
     ierr is an error indicator:;
     0 - pass completed without error;
     2 - ifac not catered for;
     3 - ifac only catered for if la=n/ifac;
   */

  long i0, i1, i2, i3, i4, i5, i6, i7;
  long j0, j1, j2, j3, j4, j5, j6, j7;
  long ia, ib, ic, id, ie, iF;
  long ja, jb, jc, jd, je, jf;
  long i, j, k, l, m, ijk;
  long ibase, jbase;
  long iink, jink;
  long jump;
  long kstop;
  long kb, kc, kd, ke, kf;

  double c1, c2, c3, c4, c5;
  double s1, s2, s3, s4, s5;
  double qqrt5;
  double ssin36;
  double ssin72;

  double a10, a11, a20, a21;
  double b10, b11, b20, b21;

  m = n / ifac;
  iink = la * inc1;
  jink = la * inc2;
  jump = (ifac - 1) * jink;
  kstop = (n - ifac) / (2 * ifac);
  ibase = 0;
  jbase = 0;

  switch (ifac)
    {
    case 2:
      {
	double a0m1, b0p1;

	i0 = j0 = 0;
	i1 = i0 + inc1 * (m + m - la);
	j1 = j0 + jink;
	if (la != m)
	  {
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    c[j0 + j] = a[i0 + i] + a[i1 + i];
		    c[j1 + j] = a[i0 + i] - a[i1 + i];
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	    i0 += iink;
	    iink += iink;
	    i1 -= iink;
	    ibase = 0;
	    jbase += jump;
	    jump += jump + jink;

	    if (i0 != i1)
	      {
		for (k = la; k <= kstop; k += la)
		  {
		    kb = k + k;
		    c1 = trigs[kb  ];
		    s1 = trigs[kb+1];
		    ibase = 0;
		    for (l = 0; l < la; ++l)
		      {
			i = ibase;
			j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
			for (ijk = 0; ijk < lot; ++ijk)
			  {
			    a0m1 = a[i0 + i] - a[i1 + i];
			    b0p1 = b[i0 + i] + b[i1 + i];

			    c[j0 + j] = a[i0 + i] + a[i1 + i];
			    d[j0 + j] = b[i0 + i] - b[i1 + i];
			    c[j1 + j] = c1 * a0m1 - s1 * b0p1;
			    d[j1 + j] = s1 * a0m1 + c1 * b0p1;
			    i += inc3;
			    j += inc4;
			  }
			ibase += inc1;
			jbase += inc2;
		      }
		    i0 += iink;
		    i1 -= iink;
		    jbase += jump;
		  }		/* End FORK */
		if (i0 > i1)
		  return 0;
	      }			/* End (i0 != i1) */
	    ibase = 0;
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    c[j0 + j] = a[i0 + i];
		    c[j1 + j] = -b[i0 + i];
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	  }
	else			/* (la != m) */
	  {
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    c[j0 + j] = 2.0 * (a[i0 + i] + a[i1 + i]);
		    c[j1 + j] = 2.0 * (a[i0 + i] - a[i1 + i]);
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	  }
	return 0;
      }

    case 3:
      {
	double afa1, a1p2, a1m2, a0mm, a0mp;
	double bfa1, b1p2, b1m2, b0mm, b0mp;

	i0 = j0 = 0;
	i1 = i0 + inc1 * (m + m - la);
	i2 = i1;
	j1 = j0 + jink;
	j2 = j1 + jink;

	if (la != m)
	  {
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    afa1 = a[i0 + i] - 0.5 * a[i1 + i];
		    bfa1 = S60 * b[i1 + i];

		    c[j0 + j] = a[i0 + i] + a[i1 + i];
		    c[j1 + j] = afa1 - bfa1;
		    c[j2 + j] = afa1 + bfa1;
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	    i0 += iink;
	    iink += iink;
	    i1 += iink;
	    i2 -= iink;
	    jbase += jump;
	    jump += jump + jink;

	    if (i0 != i2)
	      {
		for (k = la; k <= kstop; k += la)
		  {
		    kb = k + k;
		    kc = kb + kb;
		    c1 = trigs[kb  ];
		    s1 = trigs[kb+1];
		    c2 = trigs[kc  ];
		    s2 = trigs[kc+1];
		    ibase = 0;
		    for (l = 0; l < la; ++l)
		      {
			i = ibase;
			j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
			for (ijk = 0; ijk < lot; ++ijk)
			  {
			    a1p2 = a[i0 + i] - 0.5 * (a[i1 + i] + a[i2 + i]);
			    b1m2 = b[i0 + i] - 0.5 * (b[i1 + i] - b[i2 + i]);
			    a1m2 = S60 * (a[i1 + i] - a[i2 + i]);
			    b1p2 = S60 * (b[i1 + i] + b[i2 + i]);

			    a0mm = a1p2 - b1p2;
			    a0mp = a1p2 + b1p2;
			    b0mm = b1m2 - a1m2;
			    b0mp = b1m2 + a1m2;

			    c[j0 + j] = a[i0 + i] + a[i1 + i] + a[i2 + i];
			    d[j0 + j] = b[i0 + i] + b[i1 + i] - b[i2 + i];
			    c[j1 + j] = c1 * a0mm - s1 * b0mp;
			    d[j1 + j] = s1 * a0mm + c1 * b0mp;
			    c[j2 + j] = c2 * a0mp - s2 * b0mm;
			    d[j2 + j] = s2 * a0mp + c2 * b0mm;
			    i += inc3;
			    j += inc4;
			  }
			ibase += inc1;
			jbase += inc2;
		      }
		    i0 += iink;
		    i1 += iink;
		    i2 -= iink;
		    jbase += jump;
		  }		/* End FORK */
		if (i0 > i2)
		  return 0;
	      }			/* End (i0 != i2) */
	    ibase = 0;
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    a0mp = 0.5 * a[i0 + i];
		    b0mp = S60 * b[i0 + i];

		    c[j0 + j] = a[i0 + i] + a[i1 + i];
		    c[j1 + j] = a0mp - a[i1 + i] - b0mp;
		    c[j2 + j] = a[i1 + i] - a0mp - b0mp;
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	  }
	else			/* (la != m) */
	  {
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    a0mp = 2.0 * a[i0 + i] - a[i1 + i];
		    b0mp = D60 * b[i1 + i];

		    c[j0 + j] = 2.0 * (a[i0 + i] + a[i1 + i]);
		    c[j1 + j] = a0mp - b0mp;
		    c[j2 + j] = a0mp + b0mp;
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	  }
	return 0;
      }

    case 4:
      {
	double a0m1, a0p2, a1p3, a0m2, a1m3, a0p2ma1p3, a0m2pb1p3, a0m2mb1p3;
	double b0p1, b0p2, b1p3, b0m2, b1m3, b0p2pa1m3, b0p2ma1m3, b0m2mb1m3;

	i0 = j0 = 0;
	i1 = i3 = i0 + inc1 * (m + m - la);
	i2 = i1 + inc1 * (m + m);
	j1 = j0 + jink;
	j2 = j1 + jink;
	j3 = j2 + jink;

	if (la != m)
	  {
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    a0p2 = a[i0 + i] + a[i2 + i];
		    a0m2 = a[i0 + i] - a[i2 + i];

		    c[j0 + j] = a0p2 + a[i1 + i];
		    c[j1 + j] = a0m2 - b[i1 + i];
		    c[j2 + j] = a0p2 - a[i1 + i];
		    c[j3 + j] = a0m2 + b[i1 + i];
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	    i0 += iink;
	    iink += iink;
	    i1 += iink;
	    i2 -= iink;
	    i3 -= iink;
	    jbase += jump;
	    jump += jump + jink;

	    if (i1 != i2)
	      {
		for (k = la; k <= kstop; k += la)
		  {
		    kb = k + k;
		    kc = kb + kb;
		    kd = kc + kb;
		    c1 = trigs[kb  ];
		    s1 = trigs[kb+1];
		    c2 = trigs[kc  ];
		    s2 = trigs[kc+1];
		    c3 = trigs[kd  ];
		    s3 = trigs[kd+1];
		    ibase = 0;
		    for (l = 0; l < la; ++l)
		      {
			i = ibase;
			j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
			for (ijk = 0; ijk < lot; ++ijk)
			  {
			    a0p2 = a[i0 + i] + a[i2 + i];
			    a0m2 = a[i0 + i] - a[i2 + i];
			    a1p3 = a[i1 + i] + a[i3 + i];
			    a1m3 = a[i1 + i] - a[i3 + i];
			    b0p2 = b[i0 + i] + b[i2 + i];
			    b0m2 = b[i0 + i] - b[i2 + i];
			    b1p3 = b[i1 + i] + b[i3 + i];
			    b1m3 = b[i1 + i] - b[i3 + i];

			    a0p2ma1p3 = a0p2 - a1p3;
			    a0m2pb1p3 = a0m2 + b1p3;
			    a0m2mb1p3 = a0m2 - b1p3;
			    b0p2pa1m3 = b0p2 + a1m3;
			    b0p2ma1m3 = b0p2 - a1m3;
			    b0m2mb1m3 = b0m2 - b1m3;

			    c[j0 + j] = a0p2 + a1p3;
			    d[j0 + j] = b0m2 + b1m3;
			    c[j2 + j] = c2 * a0p2ma1p3 - s2 * b0m2mb1m3;
			    d[j2 + j] = s2 * a0p2ma1p3 + c2 * b0m2mb1m3;
			    c[j1 + j] = c1 * a0m2mb1p3 - s1 * b0p2pa1m3;
			    d[j1 + j] = s1 * a0m2mb1p3 + c1 * b0p2pa1m3;
			    c[j3 + j] = c3 * a0m2pb1p3 - s3 * b0p2ma1m3;
			    d[j3 + j] = s3 * a0m2pb1p3 + c3 * b0p2ma1m3;
			    i += inc3;
			    j += inc4;
			  }
			ibase += inc1;
			jbase += inc2;
		      }
		    i0 += iink;
		    i1 += iink;
		    i2 -= iink;
		    i3 -= iink;
		    jbase += jump;
		  }		/* End FORK */
		if (i1 > i2)
		  return 0;
	      }			/* End (i1 != i2) */
	    ibase = 0;
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    a0m1 = a[i0 + i] - a[i1 + i];
		    b0p1 = b[i0 + i] + b[i1 + i];

		    c[j0 + j] = a[i0 + i] + a[i1 + i];
		    c[j2 + j] = b[i1 + i] - b[i0 + i];

		    c[j1 + j] = SQ2 * (a0m1 - b0p1);
		    c[j3 + j] = -SQ2 * (a0m1 + b0p1);
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	  }
	else			/* (la != m) */
	  {
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    a0p2 = a[i0 + i] + a[i2 + i];
		    a0m2 = a[i0 + i] - a[i2 + i];

		    c[j0 + j] = 2.0 * (a0p2 + a[i1 + i]);
		    c[j1 + j] = 2.0 * (a0m2 - b[i1 + i]);
		    c[j2 + j] = 2.0 * (a0p2 - a[i1 + i]);
		    c[j3 + j] = 2.0 * (a0m2 + b[i1 + i]);
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	  }
	return 0;
      }

    case 5:
      {
	double a1p2, a1m2, a0mm, a0mp, b136, b172, b236, b272;

	i0 = j0 = 0;
	i1 = i4 = i0 + inc1 * (m + m - la);
	i2 = i3 = i1 + inc1 * (m + m);
	j1 = j0 + jink;
	j2 = j1 + jink;
	j3 = j2 + jink;
	j4 = j3 + jink;

	if (la != m)
	  {
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    a1p2 = QUA * (a[i1 + i] + a[i2 + i]);
		    a1m2 = QT5 * (a[i1 + i] - a[i2 + i]);

		    a0mp = a[i0 + i] - a1p2 + a1m2;
		    a0mm = a[i0 + i] - a1p2 - a1m2;

		    b136 = b[i1 + i] * S36;
		    b172 = b[i1 + i] * S72;
		    b236 = b[i2 + i] * S36;
		    b272 = b[i2 + i] * S72;

		    c[j0 + j] = a[i0 + i] + a[i1 + i] + a[i2 + i];
		    c[j1 + j] = a0mp - b172 - b236;
		    c[j2 + j] = a0mm - b136 + b272;
		    c[j3 + j] = a0mm + b136 - b272;
		    c[j4 + j] = a0mp + b172 + b236;
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	    i0 += iink;
	    iink += iink;
	    i1 += iink;
	    i2 += iink;
	    i3 -= iink;
	    i4 -= iink;
	    jbase += jump;
	    jump += jump + jink;

	    if (i1 != i3)
	      {
		for (k = la; k <= kstop; k += la)
		  {
		    kb = k + k;
		    kc = kb + kb;
		    kd = kc + kb;
		    ke = kd + kb;
		    c1 = trigs[kb  ];
		    s1 = trigs[kb+1];
		    c2 = trigs[kc  ];
		    s2 = trigs[kc+1];
		    c3 = trigs[kd  ];
		    s3 = trigs[kd+1];
		    c4 = trigs[ke  ];
		    s4 = trigs[ke+1];
		    ibase = 0;
		    for (l = 0; l < la; ++l)
		      {
			i = ibase;
			j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
			for (ijk = 0; ijk < lot; ++ijk)
			  {
			    a10 = (a[i0 + i] - 0.25 * ((a[i1 + i] + a[i4 + i]) +
				                       (a[i2 + i] + a[i3 + i]))) +
			                        QT5 * ((a[i1 + i] + a[i4 + i]) -
				                       (a[i2 + i] + a[i3 + i]));
			    a20 = (a[i0 + i] - 0.25 * ((a[i1 + i] + a[i4 + i]) +
			                	       (a[i2 + i] + a[i3 + i]))) -
			                        QT5 * ((a[i1 + i] + a[i4 + i]) -
				                       (a[i2 + i] + a[i3 + i]));
			    b10 = (b[i0 + i] - 0.25 * ((b[i1 + i] - b[i4 + i]) +
						       (b[i2 + i] - b[i3 + i]))) +
			                        QT5 * ((b[i1 + i] - b[i4 + i]) -
						       (b[i2 + i] - b[i3 + i]));
			    b20 = (b[i0 + i] - 0.25 * ((b[i1 + i] - b[i4 + i]) +
						       (b[i2 + i] - b[i3 + i]))) -
			                        QT5 * ((b[i1 + i] - b[i4 + i]) -
						       (b[i2 + i] - b[i3 + i]));

			    a11 = S72 * (b[i1 + i] + b[i4 + i]) +
			          S36 * (b[i2 + i] + b[i3 + i]);
			    a21 = S36 * (b[i1 + i] + b[i4 + i]) -
			          S72 * (b[i2 + i] + b[i3 + i]);
			    b11 = S72 * (a[i1 + i] - a[i4 + i]) +
			          S36 * (a[i2 + i] - a[i3 + i]);
			    b21 = S36 * (a[i1 + i] - a[i4 + i]) -
			          S72 * (a[i2 + i] - a[i3 + i]);

			    c[j0 + j] = a[i0 + i] + ((a[i1 + i] + a[i4 + i]) +
						     (a[i2 + i] + a[i3 + i]));
			    d[j0 + j] = b[i0 + i] + ((b[i1 + i] - b[i4 + i]) +
						     (b[i2 + i] - b[i3 + i]));
			    c[j1 + j] = c1 * (a10 - a11) - s1 * (b10 + b11);
			    d[j1 + j] = s1 * (a10 - a11) + c1 * (b10 + b11);
			    c[j4 + j] = c4 * (a10 + a11) - s4 * (b10 - b11);
			    d[j4 + j] = s4 * (a10 + a11) + c4 * (b10 - b11);
			    c[j2 + j] = c2 * (a20 - a21) - s2 * (b20 + b21);
			    d[j2 + j] = s2 * (a20 - a21) + c2 * (b20 + b21);
			    c[j3 + j] = c3 * (a20 + a21) - s3 * (b20 - b21);
			    d[j3 + j] = s3 * (a20 + a21) + c3 * (b20 - b21);
			    i += inc3;
			    j += inc4;
			  }
			ibase += inc1;
			jbase += inc2;
		      }
		    i0 += iink;
		    i1 += iink;
		    i2 += iink;
		    i3 -= iink;
		    i4 -= iink;
		    jbase += jump;
		  }		/* End FORK */
		if (i1 > i3)
		  return 0;
	      }			/* End (i1 != i3) */
	    ibase = 0;
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    c[j0 + j] = a[i0 + i] + a[i1 + i] + a[i2 + i];
		    c[j1 + j] = (QT5 * (a[i0 + i] - a[i1 + i])
			     + (0.25 * (a[i0 + i] + a[i1 + i]) -
				        a[i2 + i])) - (S36 * b[i0 + i] +
						       S72 * b[i1 + i]);
		    c[j4 + j] = -(QT5 * (a[i0 + i] - a[i1 + i]) +
				(0.25 * (a[i0 + i] + a[i1 + i]) - a[i2 + i])) -
		                  (S36 * b[i0 + i] + S72 * b[i1 + i]);
		    c[j2 + j] =  (QT5 * (a[i0 + i] - a[i1 + i]) -
				(0.25 * (a[i0 + i] + a[i1 + i]) - a[i2 + i])) -
		                  (S72 * b[i0 + i] - S36 * b[i1 + i]);
		    c[j3 + j] = -(QT5 * (a[i0 + i] - a[i1 + i]) -
				(0.25 * (a[i0 + i] + a[i1 + i]) - a[i2 + i])) -
		                  (S72 * b[i0 + i] - S36 * b[i1 + i]);
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	  }
	else
	  {
	    qqrt5 = 2.0 * QT5;
	    ssin36 = 2.0 * S36;
	    ssin72 = 2.0 * S72;
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    c[j0 + j] = 2.0 * (a[i0 + i] + a[i1 + i] + a[i2 + i]);
		    c[j1 + j] =(2.0 * (a[i0 + i] - 0.25 * (a[i1 + i] + a[i2 + i])) +
		       qqrt5 * (a[i1 + i] - a[i2 + i])) - (ssin72 * b[i1 + i] +
							   ssin36 * b[i2 + i]);
		    c[j2 + j] =(2.0 * (a[i0 + i] - 0.25 * (a[i1 + i] + a[i2 + i])) -
		       qqrt5 * (a[i1 + i] - a[i2 + i])) - (ssin36 * b[i1 + i] -
							   ssin72 * b[i2 + i]);
		    c[j3 + j] =(2.0 * (a[i0 + i] - 0.25 * (a[i1 + i] + a[i2 + i])) -
		       qqrt5 * (a[i1 + i] - a[i2 + i])) + (ssin36 * b[i1 + i] -
							   ssin72 * b[i2 + i]);
		    c[j4 + j] =(2.0 * (a[i0 + i] - 0.25 * (a[i1 + i] + a[i2 + i])) +
		       qqrt5 * (a[i1 + i] - a[i2 + i])) + (ssin72 * b[i1 + i] +
							   ssin36 * b[i2 + i]);
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	  }
	return 0;
      }

    case 6:
      {
	ia = 0;
	ib = ia + (2 * m - la) * inc1;
	ic = ib + 2 * m * inc1;
	id = ic + 2 * m * inc1;
	ie = ic;
	iF = ib;
	ja = 0;
	jb = ja + jink;
	jc = jb + jink;
	jd = jc + jink;
	je = jd + jink;
	jf = je + jink;

	if (la != m)		/* go to 690 */
	  {
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    c[ja + j] = (a[ia + i] + a[id + i]) + (a[ib + i] + a[ic + i]);
		    c[jd + j] = (a[ia + i] - a[id + i]) - (a[ib + i] - a[ic + i]);
		    c[jb + j] =((a[ia + i] - a[id + i]) +
		          0.5 * (a[ib + i] - a[ic + i])) - S60 * (b[ib + i] +
								  b[ic + i]);
		    c[jf + j] =((a[ia + i] - a[id + i]) +
		          0.5 * (a[ib + i] - a[ic + i])) + S60 * (b[ib + i] +
								  b[ic + i]);
		    c[jc + j] =((a[ia + i] + a[id + i]) -
		          0.5 * (a[ib + i] + a[ic + i])) - S60 * (b[ib + i] -
								  b[ic + i]);
		    c[je + j] =((a[ia + i] + a[id + i]) -
		          0.5 * (a[ib + i] + a[ic + i])) + S60 * (b[ib + i] -
								  b[ic + i]);
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	    ia += iink;
	    iink += iink;
	    ib += iink;
	    ic += iink;
	    id -= iink;
	    ie -= iink;
	    iF -= iink;
	    jbase += jump;
	    jump += jump + jink;

	    if (ic != id)	/* go to 660 */
	      {
		for (k = la; k <= kstop; k += la)
		  {
		    kb = k + k;
		    kc = kb + kb;
		    kd = kc + kb;
		    ke = kd + kb;
		    kf = ke + kb;
		    c1 = trigs[kb  ];
		    s1 = trigs[kb+1];
		    c2 = trigs[kc  ];
		    s2 = trigs[kc+1];
		    c3 = trigs[kd  ];
		    s3 = trigs[kd+1];
		    c4 = trigs[ke  ];
		    s4 = trigs[ke+1];
		    c5 = trigs[kf  ];
		    s5 = trigs[kf+1];
		    ibase = 0;
		    for (l = 0; l < la; ++l)
		      {
			i = ibase;
			j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
			for (ijk = 0; ijk < lot; ++ijk)
			  {
			    a11 = a[ie + i] + a[ib + i] + a[ic + i] + a[iF + i];
			    a20 = a[ia + i] + a[id + i] - 0.5 * a11;
			    a21 = S60 * ((a[ie + i] + a[ib + i]) -
					 (a[ic + i] + a[iF + i]));
			    b11 = b[ib + i] - b[ie + i] + b[ic + i] - b[iF + i];
			    b20 = b[ia + i] - b[id + i] - 0.5 * b11;
			    b21 = S60 * ((b[ib + i] - b[ie + i]) -
					 (b[ic + i] - b[iF + i]));

			    c[ja + j] = a[ia + i] + a[id + i] + a11;
			    d[ja + j] = b[ia + i] - b[id + i] + b11;
			    c[jc + j] = c2 * (a20 - b21) - s2 * (b20 + a21);
			    d[jc + j] = s2 * (a20 - b21) + c2 * (b20 + a21);
			    c[je + j] = c4 * (a20 + b21) - s4 * (b20 - a21);
			    d[je + j] = s4 * (a20 + b21) + c4 * (b20 - a21);

			    a11 = (a[ie + i] - a[ib + i]) + (a[ic + i] - a[iF + i]);
			    b11 = (b[ie + i] + b[ib + i]) - (b[ic + i] + b[iF + i]);
			    a20 = (a[ia + i] - a[id + i]) - 0.5 * a11;
			    a21 = S60 * ((a[ie + i] - a[ib + i]) -
					 (a[ic + i] - a[iF + i]));
			    b20 = (b[ia + i] + b[id + i]) + 0.5 * b11;
			    b21 = S60 * ((b[ie + i] + b[ib + i]) +
					 (b[ic + i] + b[iF + i]));

			    c[jd + j] = c3 * (a[ia + i] - a[id + i] + a11)
			              - s3 * (b[ia + i] + b[id + i] - b11);
			    d[jd + j] = s3 * (a[ia + i] - a[id + i] + a11)
			              + c3 * (b[ia + i] + b[id + i] - b11);
			    c[jb + j] = c1 * (a20 - b21) - s1 * (b20 - a21);
			    d[jb + j] = s1 * (a20 - b21) + c1 * (b20 - a21);
			    c[jf + j] = c5 * (a20 + b21) - s5 * (b20 + a21);
			    d[jf + j] = s5 * (a20 + b21) + c5 * (b20 + a21);
			    i += inc3;
			    j += inc4;
			  }
			ibase += inc1;
			jbase += inc2;
		      }
		    ia += iink;
		    ib += iink;
		    ic += iink;
		    id -= iink;
		    ie -= iink;
		    iF -= iink;
		    jbase += jump;
		  }
		if (ic > id)
		  return 0;
	      }
	    ibase = 0;
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    c[ja + j] =  a[ib + i] + (a[ia + i] + a[ic + i]);
		    c[jd + j] =  b[ib + i] - (b[ia + i] + b[ic + i]);
		    c[jb + j] =  (S60 * (a[ia + i] - a[ic + i])) -
		                 (0.5 * (b[ia + i] + b[ic + i]) + b[ib + i]);
		    c[jf + j] = -(S60 * (a[ia + i] - a[ic + i])) -
		                 (0.5 * (b[ia + i] + b[ic + i]) + b[ib + i]);
		    c[jc + j] =   S60 * (b[ic + i] - b[ia + i]) +
		                 (0.5 * (a[ia + i] + a[ic + i]) - a[ib + i]);
		    c[je + j] =   S60 * (b[ic + i] - b[ia + i]) -
		                 (0.5 * (a[ia + i] + a[ic + i]) - a[ib + i]);
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	  }
	else
	  {
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    c[ja + j] = (2.0 * (a[ia + i] + a[id + i])) +
		                (2.0 * (a[ib + i] + a[ic + i]));
		    c[jd + j] = (2.0 * (a[ia + i] - a[id + i])) -
		                (2.0 * (a[ib + i] - a[ic + i]));
		    c[jb + j] = (2.0 * (a[ia + i] - a[id + i]) +
		                       (a[ib + i] - a[ic + i])) -
		                (D60 * (b[ib + i] + b[ic + i]));
		    c[jf + j] = (2.0 * (a[ia + i] - a[id + i]) +
		                       (a[ib + i] - a[ic + i])) +
		                (D60 * (b[ib + i] + b[ic + i]));
		    c[jc + j] = (2.0 * (a[ia + i] + a[id + i]) -
			               (a[ib + i] + a[ic + i])) -
		                (D60 * (b[ib + i] - b[ic + i]));
		    c[je + j] = (2.0 * (a[ia + i] + a[id + i]) -
		                       (a[ib + i] + a[ic + i])) +
                                (D60 * (b[ib + i] - b[ic + i]));
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	  }
	return 0;
      }

    case 8:
      {
	double a0p7, a1p5, a2p6, p073, p074, p152;
	double a0m7, a1m5, a2m6, m073, m074, m152;

	if (la != m)
	  return 3;
	i0 = 0;
	i1 = i0 + iink;
	i2 = i1 + iink;
	i3 = i2 + iink;
	i4 = i3 + iink;
	i5 = i4 + iink;
	i6 = i5 + iink;
	i7 = i6 + iink;
	j0 = 0;
	j1 = j0 + jink;
	j2 = j1 + jink;
	j3 = j2 + jink;
	j4 = j3 + jink;
	j5 = j4 + jink;
	j6 = j5 + jink;
	j7 = j6 + jink;

	for (l = 0; l < la; ++l)
	  {
	    i = ibase;
	    j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
	    for (ijk = 0; ijk < lot; ++ijk)
	      {
		a0p7 = a[i0 + i] + a[i7 + i];
		a0m7 = a[i0 + i] - a[i7 + i];
		a1p5 = a[i1 + i] + a[i5 + i];
		a1m5 = a[i1 + i] - a[i5 + i];
		a2p6 = a[i2 + i] + a[i6 + i];
		a2m6 = a[i2 + i] - a[i6 + i];

		p073 = a0p7 + a[i3 + i];
		m073 = a0p7 - a[i3 + i];

		p074 = 2.0 * (a0m7 + a[i4 + i]);
		m074 = 2.0 * (a0m7 - a[i4 + i]);

		p152 = M_SQRT2 * (a1m5 + a2p6);
		m152 = M_SQRT2 * (a1m5 - a2p6);

		c[j0 + j] = 2.0 * (p073 + a1p5);
		c[j4 + j] = 2.0 * (p073 - a1p5);
		c[j2 + j] = 2.0 * (m073 - a2m6);
		c[j6 + j] = 2.0 * (m073 + a2m6);

		c[j1 + j] = m074 + m152;
		c[j5 + j] = m074 - m152;
		c[j3 + j] = p074 - p152;
		c[j7 + j] = p074 + p152;
		i += inc3;
		j += inc4;
	      }
	    ibase += inc1;
	    jbase += inc2;
	  }
      }
    }
  return 0;
}

static
int qpassc(double *a, double *b, double *c, double *d, double *trigs,
	   long inc1, long inc2, long inc3, long inc4,
	   long lot, long n, long ifac, long la)
{
  /*
     qpassc - performs one pass through data as part
     of multiple real fft (fourier analysis) routine.

     a      is first real input vector
     b      is equivalent to a + ifac * la * inc1
     c      is first real output vector;
     d      is equivalent to c + la * inc2
     trigs  is a precalculated list of sines & cosines
     inc1   is the addressing increment for a
     inc2   is the addressing increment for c
     inc3   is the increment between input vectors a
     inc4   is the increment between output vectors c
     lot    is the number of vectors
     n      is the length of the vectors
     ifac   is the current factor of n
     la     is the product of previous factors
   */

  long i0, i1, i2, i3, i4, i5, i6, i7;
  long j0, j1, j2, j3, j4, j5, j6, j7;
  long ia, ib, ic;
  long ja, jb, jc;
  long i, j, k, l, m, ijk;
  long ibase, jbase;
  long iink, jink;
  long jump;
  long kstop;
  long kb, kc, kd, ke, kf;

  double a0, a1, a2, a3;
  double b0, b1, b2, b3;
  double c1, c2, c3, c4, c5;
  double s1, s2, s3, s4, s5;
  double w, x, y, z;

  m = n / ifac;
  iink = la * inc1;
  jink = la * inc2;
  jump = (ifac - 1) * iink;
  kstop = (n - ifac) / (2 * ifac);
  ibase = 0;
  jbase = 0;

  switch (ifac)
    {
    case 2:
      {
	i0 = j0 = 0;
	i1 = i0 + iink;
	j1 = j0 + inc2 * (m + m - la);
	if (la != m)
	  {
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    c[j0 + j] = a[i0 + i] + a[i1 + i];
		    c[j1 + j] = a[i0 + i] - a[i1 + i];
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	    j0    += jink;
	    jink  += jink;
	    j1    -= jink;
	    ibase += jump;
	    jump  += jump + iink;

	    if (j0 != j1)
	      {
		for (k = la; k <= kstop; k += la)
		  {
		    kb = k + k;
		    c1 = trigs[kb  ];
		    s1 = trigs[kb+1];
		    jbase = 0;
		    for (l = 0; l < la; ++l)
		      {
			i = ibase;
			j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
			for (ijk = 0; ijk < lot; ++ijk)
			  {
			    c[j0 + j] = a[i0 + i] + c1 * a[i1 + i] + s1 * b[i1 + i];
			    c[j1 + j] = a[i0 + i] - c1 * a[i1 + i] - s1 * b[i1 + i];
			    d[j0 + j] = c1 * b[i1 + i] - s1 * a[i1 + i] + b[i0 + i];
			    d[j1 + j] = c1 * b[i1 + i] - s1 * a[i1 + i] - b[i0 + i];
			    i += inc3;
			    j += inc4;
			  }
			ibase += inc1;
			jbase += inc2;
		      }
		    j0 += jink;
		    j1 -= jink;
		    ibase += jump;
		  }		/* End FORK */
		if (j0 > j1)
		  return 0;
	      }			/* End (i0 != i1) */
	    jbase = 0;
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    c[j0 + j] = a[i0 + i];
		    d[j1 + j] = -a[i1 + i];
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	  }
	else			/* (la != m) */
	  {
	    z = 1.0 / n;
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    c[j0 + j] = z * (a[i0 + i] + a[i1 + i]);
		    c[j1 + j] = z * (a[i0 + i] - a[i1 + i]);
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	  }
	return 0;
      }

    case 3:
      {
	ia = 0;
	ib = ia + iink;
	ic = ib + iink;

	ja = 0;
	jb = ja + inc2 * (m + m - la);
	jc = jb;

	if (la != m)		/* else 390 */
	  {
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    c[ja + j] = a[ia + i] + a[ib + i] + a[ic + i];
		    c[jb + j] = a[ia + i] - 0.5 * (a[ib + i] + a[ic + i]);
		    d[jb + j] = S60 * (a[ic + i] - a[ib + i]);
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	    ja += jink;
	    jink += jink;
	    jb += jink;
	    jc -= jink;
	    ibase += jump;
	    jump += jump + iink;

	    if (ja != jc)	/* else  360 */
	      {
		for (k = la; k <= kstop; k += la)
		  {
		    kb = k + k;
		    kc = kb + kb;
		    c1 = trigs[kb  ];
		    s1 = trigs[kb+1];
		    c2 = trigs[kc  ];
		    s2 = trigs[kc+1];
		    jbase = 0;
		    for (l = 0; l < la; ++l)
		      {
			i = ibase;
			j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
			for (ijk = 0; ijk < lot; ++ijk)
			  {
			    a1 = c1 * a[ib + i] + s1 * b[ib + i] +
                                 c2 * a[ic + i] + s2 * b[ic + i];
			    b1 = c1 * b[ib + i] - s1 * a[ib + i] +
                                 c2 * b[ic + i] - s2 * a[ic + i];
			    a2 = a[ia + i] - 0.5 * a1;
			    b2 = b[ia + i] - 0.5 * b1;
			    a3 = S60 * (c1 * a[ib + i] + s1 * b[ib + i] -
				        c2 * a[ic + i] - s2 * b[ic + i]);
			    b3 = S60 * (c1 * b[ib + i] - s1 * a[ib + i] -
				        c2 * b[ic + i] + s2 * a[ic + i]);

			    c[ja + j] = a[ia + i] + a1;
			    d[ja + j] = b[ia + i] + b1;
			    c[jb + j] = a2 + b3;
			    d[jb + j] = b2 - a3;
			    c[jc + j] = a2 - b3;
			    d[jc + j] = -b2 - a3;
			    i += inc3;
			    j += inc4;
			  }
			ibase += inc1;
			jbase += inc2;
		      }
		    ja += jink;
		    jb += jink;
		    jc -= jink;
		    ibase += jump;
		  }		/* End FORK */
		if (ja > jc)
		  return 0;
	      }			/* End (ia != ic) */
	    jbase = 0;
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    /* soweit */
		    c[ja + j] = a[ia + i] + 0.5 * (a[ib + i] - a[ic + i]);
		    d[ja + j] = -S60 * (a[ib + i] + a[ic + i]);
		    c[jb + j] = a[ia + i] - a[ib + i] + a[ic + i];
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	  }
	else			/* (la != m) */
	  {
	    z = 1.0 / n;
	    y = S60 / n;
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    c[ja + j] = z * (a[ia + i] + a[ib + i] + a[ic + i]);
		    c[jb + j] = z * (a[ia + i] - 0.5 * (a[ib + i] + a[ic + i]));
		    d[jb + j] = y * (a[ic + i] - a[ib + i]);
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	  }
	return 0;
      }

    case 4:
      {
	double a0p2, a1p3;

	i0 = 0;
	i1 = i0 + iink;
	i2 = i1 + iink;
	i3 = i2 + iink;
	j0 = 0;
	j1 = j0 + inc2 * (m + m - la);
	j2 = j1 + inc2 * (m + m);
	j3 = j1;

	if (la != m)		/*else go to 490 */
	  {
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    a0p2 = a[i0 + i] + a[i2 + i];
		    a1p3 = a[i1 + i] + a[i3 + i];

		    c[j0 + j] = a0p2 + a1p3;
		    c[j2 + j] = a0p2 - a1p3;

		    c[j1 + j] = a[i0 + i] - a[i2 + i];
		    d[j1 + j] = a[i3 + i] - a[i1 + i];
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	    j0 += jink;
	    jink += jink;
	    j1 += jink;
	    j2 -= jink;
	    j3 -= jink;
	    ibase += jump;
	    jump += jump + iink;

	    if (j1 != j2)	/* else go to 460; */
	      {
		for (k = la; k <= kstop; k += la)
		  {
		    kb = k + k;
		    kc = kb + kb;
		    kd = kc + kb;
		    c1 = trigs[kb  ];
		    s1 = trigs[kb+1];
		    c2 = trigs[kc  ];
		    s2 = trigs[kc+1];
		    c3 = trigs[kd  ];
		    s3 = trigs[kd+1];
		    jbase = 0;
		    for (l = 0; l < la; ++l)
		      {
			i = ibase;
			j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
			for (ijk = 0; ijk < lot; ++ijk)
			  {
			    a0 = a[i0 + i] + c2 * a[i2 + i] + s2 * b[i2 + i];
			    a2 = a[i0 + i] - c2 * a[i2 + i] - s2 * b[i2 + i];
			    b0 = b[i0 + i] + c2 * b[i2 + i] - s2 * a[i2 + i];
			    b2 = b[i0 + i] - c2 * b[i2 + i] + s2 * a[i2 + i];

			    a1 = c1 * a[i1 + i] + s1 * b[i1 + i] +
                                 c3 * a[i3 + i] + s3 * b[i3 + i];
			    a3 = c1 * a[i1 + i] + s1 * b[i1 + i] -
                                 c3 * a[i3 + i] - s3 * b[i3 + i];
			    b1 = c1 * b[i1 + i] - s1 * a[i1 + i] +
                                 c3 * b[i3 + i] - s3 * a[i3 + i];
			    b3 = c1 * b[i1 + i] - s1 * a[i1 + i] -
                                 c3 * b[i3 + i] + s3 * a[i3 + i];

			    c[j0 + j] = a0 + a1;
			    c[j2 + j] = a0 - a1;
			    d[j0 + j] = b0 + b1;
			    d[j2 + j] = b1 - b0;
			    c[j1 + j] = a2 + b3;
			    c[j3 + j] = a2 - b3;
			    d[j1 + j] = b2 - a3;
			    d[j3 + j] = -b2 - a3;
			    i += inc3;
			    j += inc4;
			  }
			ibase += inc1;
			jbase += inc2;
		      }
		    j0 += jink;
		    j1 += jink;
		    j2 -= jink;
		    j3 -= jink;
		    ibase += jump;
		  }		/* End FORK */
		if (j1 > j2)
		  return 0;
	      }			/* End (i1 != i2) */
	    jbase = 0;
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    c[j0 + j] =  a[i0 + i] + SQ2 * (a[i1 + i] - a[i3 + i]);
		    c[j1 + j] =  a[i0 + i] - SQ2 * (a[i1 + i] - a[i3 + i]);
		    d[j0 + j] = -a[i2 + i] - SQ2 * (a[i1 + i] + a[i3 + i]);
		    d[j1 + j] =  a[i2 + i] - SQ2 * (a[i1 + i] + a[i3 + i]);
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	  }
	else			/* (la != m) */
	  {
	    z = 1.0 / n;
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    a0p2 = a[i0 + i] + a[i2 + i];
		    a1p3 = a[i1 + i] + a[i3 + i];

		    c[j0 + j] = z * (a0p2 + a1p3);
		    c[j2 + j] = z * (a0p2 - a1p3);
		    c[j1 + j] = z * (a[i0 + i] - a[i2 + i]);
		    d[j1 + j] = z * (a[i3 + i] - a[i1 + i]);
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	  }
	return 0;
      }

    case 5:
      {
	double a1p4, a2p3, b1p4, b2p3, a025, b025, asps, bsps, a0pq, b0pq;
	double a1m4, a2m3, b1m4, b2m3, aqrt, bqrt, asms, bsms, a0mq, b0mq;

	i0 = 0;
	i1 = i0 + iink;
	i2 = i1 + iink;
	i3 = i2 + iink;
	i4 = i3 + iink;
	j0 = 0;
	j1 = j0 + inc2 * (m + m - la);
	j2 = j1 + inc2 * (m + m);
	j3 = j2;
	j4 = j1;

	if (la != m)		/* else go to 590; */
	  {
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    a1p4 = a[i1 + i] + a[i4 + i];
		    a1m4 = a[i1 + i] - a[i4 + i];
		    a2p3 = a[i2 + i] + a[i3 + i];
		    a2m3 = a[i2 + i] - a[i3 + i];

		    a025 = a[i0 + i] - 0.25 * (a1p4 + a2p3);
		    aqrt = QT5 * (a1p4 - a2p3);

		    c[j0 + j] = a[i0 + i] + a1p4 + a2p3;
		    c[j1 + j] = a025 + aqrt;
		    c[j2 + j] = a025 - aqrt;
		    d[j1 + j] = -S72 * a1m4 - S36 * a2m3;
		    d[j2 + j] = -S36 * a1m4 + S72 * a2m3;
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	    j0 += jink;
	    jink += jink;
	    j1 += jink;
	    j2 += jink;
	    j3 -= jink;
	    j4 -= jink;
	    ibase += jump;
	    jump += jump + iink;

	    if (j1 != j3)
	      {
		for (k = la; k <= kstop; k += la)
		  {
		    kb = k + k;
		    kc = kb + kb;
		    kd = kc + kb;
		    ke = kd + kb;
		    c1 = trigs[kb  ];
		    s1 = trigs[kb+1];
		    c2 = trigs[kc  ];
		    s2 = trigs[kc+1];
		    c3 = trigs[kd  ];
		    s3 = trigs[kd+1];
		    c4 = trigs[ke  ];
		    s4 = trigs[ke+1];
		    jbase = 0;
		    for (l = 0; l < la; ++l)
		      {
			i = ibase;
			j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
			for (ijk = 0; ijk < lot; ++ijk)
			  {
			    a1p4 = c1 * a[i1 + i] + s1 * b[i1 + i] + 
			           c4 * a[i4 + i] + s4 * b[i4 + i];
			    a1m4 = c1 * a[i1 + i] + s1 * b[i1 + i] - 
			           c4 * a[i4 + i] - s4 * b[i4 + i];
			    a2p3 = c2 * a[i2 + i] + s2 * b[i2 + i] + 
			           c3 * a[i3 + i] + s3 * b[i3 + i];
			    a2m3 = c2 * a[i2 + i] + s2 * b[i2 + i] -
			           c3 * a[i3 + i] - s3 * b[i3 + i];
			    b1p4 = c1 * b[i1 + i] - s1 * a[i1 + i] + 
			           c4 * b[i4 + i] - s4 * a[i4 + i];
			    b1m4 = c1 * b[i1 + i] - s1 * a[i1 + i] -
			           c4 * b[i4 + i] + s4 * a[i4 + i];
			    b2p3 = c2 * b[i2 + i] - s2 * a[i2 + i] +
			           c3 * b[i3 + i] - s3 * a[i3 + i];
			    b2m3 = c2 * b[i2 + i] - s2 * a[i2 + i] -
			           c3 * b[i3 + i] + s3 * a[i3 + i];

			    a025 = a[i0 + i] - 0.25 * (a1p4 + a2p3);
			    aqrt = QT5 * (a1p4 - a2p3);
			    b025 = b[i0 + i] - 0.25 * (b1p4 + b2p3);
			    bqrt = QT5 * (b1p4 - b2p3);

			    a0pq = a025 + aqrt;
			    a0mq = a025 - aqrt;
			    b0pq = b025 + bqrt;
			    b0mq = b025 - bqrt;

			    asps = S72 * a1m4 + S36 * a2m3;
			    asms = S36 * a1m4 - S72 * a2m3;
			    bsps = S72 * b1m4 + S36 * b2m3;
			    bsms = S36 * b1m4 - S72 * b2m3;

			    c[j0 + j] = a[i0 + i] + a1p4 + a2p3;
			    c[j1 + j] = a0pq + bsps;
			    c[j2 + j] = a0mq + bsms;
			    c[j3 + j] = a0mq - bsms;
			    c[j4 + j] = a0pq - bsps;
			    d[j0 + j] = b[i0 + i] + b1p4 + b2p3;
			    d[j1 + j] = b0pq - asps;
			    d[j2 + j] = b0mq - asms;
			    d[j3 + j] = -b0mq - asms;
			    d[j4 + j] = -b0pq - asps;
			    i += inc3;
			    j += inc4;
			  }
			ibase += inc1;
			jbase += inc2;
		      }
		    j0 += jink;
		    j1 += jink;
		    j2 += jink;
		    j3 -= jink;
		    j4 -= jink;
		    ibase += jump;
		  }		/* End FORK */
		if (j1 > j3)
		  return 0;
	      }			/* End (jb != jd) */
	    jbase = 0;
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    a1p4 = a[i1 + i] + a[i4 + i];
		    a1m4 = a[i1 + i] - a[i4 + i];
		    a2p3 = a[i2 + i] + a[i3 + i];
		    a2m3 = a[i2 + i] - a[i3 + i];

		    a025 = a[i0 + i] + 0.25 * (a1m4 - a2m3);
		    aqrt = QT5 * (a1m4 + a2m3);

		    c[j0 + j] = a025 + aqrt;
		    c[j1 + j] = a025 - aqrt;
		    c[j2 + j] = a[i0 + i] - a1m4 + a2m3;
		    d[j0 + j] = -S36 * a1p4 - S72 * a2p3;
		    d[j1 + j] = -S72 * a1p4 + S36 * a2p3;

		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	  }
	else
	  {
	    z = 1.0 / n;
	    y = QT5 / n;
	    x = S36 / n;
	    w = S72 / n;

	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    a1p4 = a[i1 + i] + a[i4 + i];
		    a1m4 = a[i1 + i] - a[i4 + i];
		    a2p3 = a[i2 + i] + a[i3 + i];
		    a2m3 = a[i2 + i] - a[i3 + i];

		    a025 = z * (a[i0 + i] - 0.25 * (a1p4 + a2p3));
		    aqrt = y * (a1p4 - a2p3);

		    c[j0 + j] = z * (a[i0 + i] + a1p4 + a2p3);
		    c[j1 + j] = a025 + aqrt;
		    c[j2 + j] = a025 - aqrt;
		    d[j1 + j] = -w * a1m4 - x * a2m3;
		    d[j2 + j] = w * a2m3 - x * a1m4;
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	  }
	return 0;
      }

    case 6:
      {
	double ab1a, ab2a, ab3a, ab4a, ab5a;
	double ab1b, ab2b, ab3b, ab4b, ab5b;
	double a0p3, a1p4, a1p5, a2p4, a2p5;
	double a0m3, a1m4, a1m5, a2m4, a2m5;
	double b1p4, b2p5;
	double b1m4, b2m5;
	double ap05, bp05, ap60, bp60;
	double am05, bm05, am60, bm60;

	i0 = 0;
	i1 = i0 + iink;
	i2 = i1 + iink;
	i3 = i2 + iink;
	i4 = i3 + iink;
	i5 = i4 + iink;
	j0 = 0;
	j1 = j0 + inc2 * (m + m - la);
	j2 = j1 + inc2 * (m + m);
	j3 = j2 + inc2 * (m + m);
	j4 = j2;
	j5 = j1;

	if (la != m)
	  {
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    a0p3 = a[i0 + i] + a[i3 + i];
		    a0m3 = a[i0 + i] - a[i3 + i];
		    a1p4 = a[i1 + i] + a[i4 + i];
		    a1m4 = a[i1 + i] - a[i4 + i];
		    a2p5 = a[i2 + i] + a[i5 + i];
		    a2m5 = a[i2 + i] - a[i5 + i];

		    c[j0 + j] = a0p3 + a1p4 + a2p5;
		    c[j3 + j] = a0m3 + a2m5 - a1m4;

		    c[j1 + j] = a0m3 - 0.5 * (a2m5 - a1m4);
		    c[j2 + j] = a0p3 - 0.5 * (a1p4 + a2p5);

		    d[j1 + j] = S60 * (-a2m5 - a1m4);
		    d[j2 + j] = S60 * (a2p5 - a1p4);
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	    j0 += jink;
	    jink += jink;
	    j1 += jink;
	    j2 += jink;
	    j3 -= jink;
	    j4 -= jink;
	    j5 -= jink;
	    ibase += jump;
	    jump += jump + iink;

	    if (j2 != j3)
	      {
		for (k = la; k <= kstop; k += la)
		  {
		    kb = k + k;
		    kc = kb + kb;
		    kd = kc + kb;
		    ke = kd + kb;
		    kf = ke + kb;
		    c1 = trigs[kb  ];
		    s1 = trigs[kb+1];
		    c2 = trigs[kc  ];
		    s2 = trigs[kc+1];
		    c3 = trigs[kd  ];
		    s3 = trigs[kd+1];
		    c4 = trigs[ke  ];
		    s4 = trigs[ke+1];
		    c5 = trigs[kf  ];
		    s5 = trigs[kf+1];
		    jbase = 0;
		    for (l = 0; l < la; ++l)
		      {
			i = ibase;
			j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
			for (ijk = 0; ijk < lot; ++ijk)
			  {
			    ab1a = c1 * a[i1 + i] + s1 * b[i1 + i];
			    ab1b = c1 * b[i1 + i] - s1 * a[i1 + i];
			    ab2a = c2 * a[i2 + i] + s2 * b[i2 + i];
			    ab2b = c2 * b[i2 + i] - s2 * a[i2 + i];
			    ab3a = c3 * a[i3 + i] + s3 * b[i3 + i];
			    ab3b = c3 * b[i3 + i] - s3 * a[i3 + i];
			    ab4a = c4 * a[i4 + i] + s4 * b[i4 + i];
			    ab4b = c4 * b[i4 + i] - s4 * a[i4 + i];
			    ab5a = c5 * a[i5 + i] + s5 * b[i5 + i];
			    ab5b = c5 * b[i5 + i] - s5 * a[i5 + i];

			    a1p4 = ab1a + ab4a;
			    a1m4 = ab1a - ab4a;
			    a2p5 = ab2a + ab5a;
			    a2m5 = ab2a - ab5a;

			    b1p4 = ab1b + ab4b;
			    b1m4 = ab1b - ab4b;
			    b2p5 = ab2b + ab5b;
			    b2m5 = ab2b - ab5b;

			    ap05 = a[i0 + i] + ab3a - 0.5 * (a1p4 + a2p5);
			    bp05 = b[i0 + i] + ab3b - 0.5 * (b1p4 + b2p5);
			    am05 = a[i0 + i] - ab3a - 0.5 * (a2m5 - a1m4);
			    bm05 = -b[i0 + i] + ab3b - 0.5 * (b1m4 - b2m5);

			    ap60 = S60 * (a2p5 - a1p4);
			    bp60 = S60 * (b2p5 - b1p4);
			    am60 = S60 * (-a2m5 - a1m4);
			    bm60 = S60 * (-b2m5 - b1m4);

			    c[j0 + j] = a[i0 + i] + ab3a + a1p4 + a2p5;
			    d[j0 + j] = b[i0 + i] + ab3b + b1p4 + b2p5;
			    c[j1 + j] = am05 - bm60;
			    d[j1 + j] = am60 - bm05;
			    c[j2 + j] = ap05 - bp60;
			    d[j2 + j] = ap60 + bp05;
			    c[j3 + j] = a[i0 + i] - ab3a - a1m4 + a2m5;
			    d[j3 + j] = -b[i0 + i] + ab3b + b1m4 - b2m5;
			    c[j4 + j] = ap05 + bp60;
			    d[j4 + j] = ap60 - bp05;
			    c[j5 + j] = am05 + bm60;
			    d[j5 + j] = am60 + bm05;
			    i += inc3;
			    j += inc4;
			  }
			ibase += inc1;
			jbase += inc2;
		      }
		    j0 += jink;
		    j1 += jink;
		    j2 += jink;
		    j3 -= jink;
		    j4 -= jink;
		    j5 -= jink;
		    ibase += jump;
		  }
		if (j2 > j3)
		  return 0;
	      }
	    jbase = 0;
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    a1p5 = a[i1 + i] + a[i5 + i];
		    a1m5 = a[i1 + i] - a[i5 + i];
		    a2p4 = a[i2 + i] + a[i4 + i];
		    a2m4 = a[i2 + i] - a[i4 + i];

		    c[j0 + j] =  a[i0 + i] + 0.5 * a2m4 + S60 * a1m5;
		    d[j0 + j] = -a[i3 + i] - 0.5 * a1p5 - S60 * a2p4;
		    c[j1 + j] =  a[i0 + i] - a2m4;
		    d[j1 + j] =  a[i3 + i] - a1p5;
		    c[j2 + j] =  a[i0 + i] + 0.5 * a2m4 - S60 * a1m5;
		    d[j2 + j] = -a[i3 + i] - 0.5 * a1p5 + S60 * a2p4;
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	  }
	else
	  {
	    z = 1.0 / n;
	    y = S60 / n;
	    for (l = 0; l < la; ++l)
	      {
		i = ibase;
		j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
		for (ijk = 0; ijk < lot; ++ijk)
		  {
		    a0p3 = a[i0 + i] + a[i3 + i];
		    a0m3 = a[i0 + i] - a[i3 + i];
		    a1p4 = a[i1 + i] + a[i4 + i];
		    a1m4 = a[i1 + i] - a[i4 + i];
		    a2p5 = a[i2 + i] + a[i5 + i];
		    a2m5 = a[i2 + i] - a[i5 + i];

		    c[j0 + j] = z * (a0p3 + a1p4 + a2p5);
		    c[j3 + j] = z * (a0m3 + a2m5 - a1m4);

		    c[j1 + j] = z * (a0m3 - 0.5 * (a2m5 - a1m4));
		    c[j2 + j] = z * (a0p3 - 0.5 * (a1p4 + a2p5));

		    d[j1 + j] = y * (-a2m5 - a1m4);
		    d[j2 + j] = y * (a2p5 - a1p4);
		    i += inc3;
		    j += inc4;
		  }
		ibase += inc1;
		jbase += inc2;
	      }
	  }
	return 0;
      }

    case 8:
      {
	double a0p4, a1p5, a2p6, a3p7;
	double a0m4, a1m5, a2m6, a3m7;

	if (la != m)
	  return 3;
	i0 = 0;
	i1 = i0 + iink;
	i2 = i1 + iink;
	i3 = i2 + iink;
	i4 = i3 + iink;
	i5 = i4 + iink;
	i6 = i5 + iink;
	i7 = i6 + iink;
	j0 = 0;
	j1 = j0 + jink;
	j2 = j1 + jink;
	j3 = j2 + jink;
	j4 = j3 + jink;
	j5 = j4 + jink;
	j6 = j5 + jink;
	j7 = j6 + jink;
	z = 1.0 / n;
	y = SQ2 / n;

	for (l = 0; l < la; ++l)
	  {
	    i = ibase;
	    j = jbase;
#if defined(SX)
#pragma vdir nodep
#endif
#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
	    for (ijk = 0; ijk < lot; ++ijk)
	      {
		a0p4 = a[i0 + i] + a[i4 + i];
		a0m4 = a[i0 + i] - a[i4 + i];
		a1p5 = a[i1 + i] + a[i5 + i];
		a1m5 = a[i1 + i] - a[i5 + i];
		a2p6 = a[i2 + i] + a[i6 + i];
		a2m6 = a[i2 + i] - a[i6 + i];
		a3p7 = a[i3 + i] + a[i7 + i];
		a3m7 = a[i3 + i] - a[i7 + i];

		c[j0 + j] =  z * (a0p4 + a1p5 + a2p6 + a3p7);
		c[j7 + j] =  z * (a0p4 - a1p5 + a2p6 - a3p7);

		c[j3 + j] =  z * (a0p4 - a2p6);
		c[j4 + j] =  z * (a3p7 - a1p5);

		c[j1 + j] =  z * a0m4 + y * (a1m5 - a3m7);
		c[j5 + j] =  z * a0m4 - y * (a1m5 - a3m7);
		c[j2 + j] = -z * a2m6 - y * (a1m5 + a3m7);
		c[j6 + j] =  z * a2m6 - y * (a1m5 + a3m7);
		i += inc3;
		j += inc4;
	      }
	    ibase += inc1;
	    jbase += inc2;
	  }
      }
    }
  return 0;
}

/* ====================== */
/* Fast Fourier Transform */
/* ====================== */
void fc2gp(double *restrict trig, long *restrict ifax, double *restrict fc, double *restrict gp, long nlat, long nlon, long nlev, long nfc)
{
  long fou, ia, ifac, k, la;
  long lat, lev, lon, rix, wix;
  double *wpt;
  long nx, nblox, nvex, nvex0, nb;
  long istart, i, j, ibase, jbase, jj, ii, ix;
  long *istartv;

  /* fc2gp performs fourier to gridpoint transforms using           */
  /* multiple fast fourier transform of length nlon                 */
  /*                                                                */
  /* fc   - real array of fourier coefficients fc[nlev][nfc][nlat]  */
  /* gp   - real array of gridpoints           gp[nlev][nlat][nlon] */
  /* nlat - Number of latitudes                                     */
  /* nlon - Number of longitudes                                    */
  /* nlev - Number of levels                                        */
  /* nfc  - Number of fourier coefficients on 1 latitude            */

  /* x(j) = sum(k=0,...,n-1)(c(k)*exp(2*i*j*k*pi/nlon))             */
  /*        where c(k) = a(k) + i*b(k) and c(n-k) = a(k)-i*b(k)     */

  if ( ifax[9] != nlon ) fprintf(stderr, "fc2gp: wrong initialization!\n");

  long nfax = ifax[0];

  long jump = (nlon + 2);
  long lot  = nlev * nlat;

  double *restrict wfc = (double*) malloc(lot*jump*sizeof(double));
#if ! defined(_OPENMP)
  double *restrict wgp = (double*) malloc(lot*jump*sizeof(double));
#endif

  for ( lev = 0; lev < nlev; ++lev )
    {
#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(lat, fou, wix, rix)
#endif
      for ( lat = 0; lat < nlat; ++lat )
	{
	  wix = jump * (lat + lev*nlat);
	  rix = lat + lev*nlat*nfc;
	  for ( fou = 0;   fou < nfc;  ++fou ) wfc[wix + fou] = fc[rix + fou*nlat];
	  for ( fou = nfc; fou < jump; ++fou ) wfc[wix + fou] = 0.0;
	  /*	  wfc[wix + 1] = 0.5 * wfc[wix]; */
	}
    }

  nx = nlon + 1;
  if ( nlon%2 == 1 ) nx = nlon;
  nblox = 1 + (lot-1)/NFFT;
  nvex  = lot - (nblox-1)*NFFT;
  nvex0 = nvex;

  istartv = (long*) malloc(nblox*sizeof(long));

  istart = 0;
  for ( nb = 0; nb < nblox; nb++ )
    {
      istartv[nb] = istart;
      istart = istart + nvex*jump;
      nvex = NFFT;
    }

#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(istart, nvex, ix, ii, jj, i, j, k, ia, la, ifac, ibase, jbase)
#endif
  for ( nb = 0; nb < nblox; nb++ )
    {
#if defined(_OPENMP)
      double *restrict wgp = (double*) malloc(lot*jump*sizeof(double));
#endif
      istart = istartv[nb];
      if ( nb == 0 ) nvex = nvex0;
      else           nvex = NFFT;

      i = istart;
#if defined(SX)
#pragma vdir nodep
#endif
      for ( j = 0; j < nvex; j++ )
	{
	  wfc[i+1] = 0.5*wfc[i];
	  i += jump;
	}
      if ( nlon%2 != 1 )
	{
	  i = istart + nlon;
	  for ( j = 0; j < nvex; j++ )
	    {
	      wfc[i] = 0.5*wfc[i];
	      i += jump;
	    }
	}

      ia = istart + 1;
      la = 1;
      for ( k = 0; k < nfax; ++k )
	{
	  ifac = ifax[k + 1];

	  if ( k & 1 )
	    rpassc(wgp, wgp+la, wfc+ia, wfc+ia+ifac*la, trig,
		   1, 1, nx, jump, nvex, nlon, ifac, la);
	  else
	    rpassc(wfc+ia, wfc+ia+la, wgp, wgp+ifac*la, trig,
		   1, 1, jump, nx, nvex, nlon, ifac, la);

	  la *= ifac;
	  ia = istart;
	}

      /* If necessary, copy results back to a */

      if ( nfax%2 != 0 )
	{
	  ibase = 0;
	  jbase = ia;
	  for ( jj = 0; jj < nvex; jj++ )
	    {
	      i = ibase;
	      j = jbase;
	      for ( ii = 0; ii < nlon; ii++ )
		{
		  wfc[j++] = wgp[i++];
		}
	      ibase = ibase + nx;
	      jbase = jbase + jump;
	    }
	 }

      /* Fill in zeros at end */

      ix = istart + nlon;
#if defined(SX)
#pragma vdir nodep
#endif
      for ( j = 0; j < nvex; j++ )
	{
          wfc[ix]   = 0.0;
          wfc[ix+1] = 0.0;
          ix = ix + jump;
	}

#if defined(_OPENMP)
      free(wgp);
#endif
    }

  wpt = wfc;

#if defined(_OPENMP)
#pragma omp parallel for default(shared) private(j, lon)
#endif
  for ( j = 0; j < lot; ++j )
    for ( lon = 0; lon < nlon; ++lon )
      gp[lon + j*nlon] = wpt[lon + j*jump];

  free(istartv);
#if ! defined(_OPENMP)
  free(wgp);
#endif
  free(wfc);
}


void gp2fc(double *trig, long *ifax, const double *restrict gp, double *restrict fc, long nlat, long nlon, long nlev, long nfc)
{
  long lot, fou, ia, ifac, jump, k, la;
  long lat, lev, lon, nfax, rix, wix;
  long nx, nblox, nvex, nb;
  long istart, i, j, ibase, jbase, jj, ii, ix, iz;

  /* gp2fc performs gridpoint to fourier transforms using           */
  /* multiple fast fourier transform of length nlon                 */
  /*                                                                */
  /* gp   - real array of gridpoints           gp[nlev][nlat][nlon] */
  /* fc   - real array of fourier coefficients fc[nlev][nfc][nlat]  */
  /* nlat - Number of latitudes                                     */
  /* nlon - Number of longitudes                                    */
  /* nlev - Number of levels                                        */
  /* nfc  - Number of fourier coefficients on 1 latitude            */

  /* a(k) =  (1/n) * sum(j=0,...,n-1)(x(j) * cos(2*j*k*pi/n))       */
  /* b(k) = -(1/n) * sum(j=0,...,n-1)(x(j) * sin(2*j*k*pi/n))       */

  if ( ifax[9] != nlon ) fprintf(stderr, "gp2fc: wrong initialization!\n");

  nfax = ifax[0];

  jump = (nlon + 2);
  lot  = nlev * nlat;

  double *restrict wfc = (double*) malloc(lot * jump * sizeof(double));
  double *restrict wgp = (double*) malloc(lot * jump * sizeof(double));

  rix = 0;
  wix = 0;
  for ( j = 0; j < lot; ++j )
    {
      for ( lon = 0; lon < nlon; ++lon )
	wgp[wix + lon] = gp[rix + lon];
      wgp[wix + nlon] = 0.0;
      wgp[wix + nlon + 1] = 0.0;
      rix += nlon;
      wix += jump;
    }

  nx = nlon + 1;
  if ( nlon%2 == 1 ) nx = nlon;
  nblox = 1 + (lot-1)/NFFT;
  nvex = lot - (nblox-1)*NFFT;

  istart = 0;
  for ( nb = 0; nb < nblox; nb++ )
    {
      ia = istart;
      la = nlon;

      for ( k = 0; k < nfax; ++k )
	{
	  ifac = ifax[nfax - k];
	  la /= ifac;
	  if (k & 1)
	    qpassc (wfc, wfc+ifac*la, wgp+ia, wgp+ia+la, trig,
		    1, 1, nx, jump, nvex, nlon, ifac, la);
	  else
	    qpassc (wgp+ia, wgp+ia+ifac*la, wfc, wfc+la, trig,
		    1, 1, jump, nx, nvex, nlon, ifac, la);
	  ia = istart + 1;
	}

      /* If necessary, copy results back to a */

      if ( nfax%2 != 0 )
	{
	  ibase = 0;
	  jbase = ia;
	  for ( jj = 0; jj < nvex; jj++ )
	    {
	      i = ibase;
	      j = jbase;
	      for ( ii = 0; ii < nlon; ii++ )
		{
		  wgp[j++] = wfc[i++];
		}
	      ibase = ibase + nx;
	      jbase = jbase + jump;
	    }
	 }

      /* Shift a(0) & fill in zero imag parts */

      ix = istart;
#if defined(SX)
#pragma vdir nodep
#endif
      for ( j = 0; j < nvex; j++ )
	{
          wgp[ix] = wgp[ix+1];
	  wgp[ix+1] = 0.0;
          ix = ix + jump;
	}

      if ( nlon%2 != 1 )
	{
	  iz = istart + (nlon+1);
	  for ( j = 0; j < nvex; j++ )
	    {
	      wgp[iz] = 0.0;
	      iz = iz + jump;
	    }
	}

      istart = istart + nvex*jump;
      nvex = NFFT;
    }

  const double *restrict wpt;
  double *restrict fct;

  for ( lev = 0; lev < nlev; ++lev )
    {
      for ( lat = 0; lat < nlat; ++lat )
	{
	  rix = jump * (lat + lev * nlat);
	  wix = lat + lev * nlat * nfc;
          wpt = wgp + rix;
          fct = fc + wix;
	  fct[0] = wpt[0];
	  fct[nlat] = 0.0;
	  for ( fou = 2; fou < nfc; ++fou )
	    fct[fou*nlat] = wpt[fou];
	}
    }

  free(wgp);
  free(wfc);
}
