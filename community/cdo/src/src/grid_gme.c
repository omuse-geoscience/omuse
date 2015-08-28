#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifndef M_PI
#define M_PI		3.14159265358979323846	/* pi */
#endif

struct dimension {
  int lower;
  int extent;
  int mult;
};

struct array {
  char *addr;
  int rank;
  int offset;
  struct dimension dim[3];  
};

struct cart {
  double x[3];
};

struct geo {
  double lon;
  double lat;
};

struct polygon {
  int type;
  struct geo center;
  struct geo boundary[6];
};

const double pid5   =   0.2*M_PI;
const double pid180 = 180.0/M_PI;

const int ispokes[12] = { 1,0,-1,-1,0,1,0,1,1,0,-1,-1,};

const int pentagon = 5;
const int hexagon  = 6;

/*****************************************************************************/

void factorni(int kni, int *kni2, int *kni3)
{
  /**********************************************************************/
  /* factorni computes the factors of the integer input kni, assuming   */ 
  /* that kni decomposes into kni3 factors (kni3 either 0 or 1) of "3"  */
  /* and kni2 (kni2 > 0) factors of "2". The subroutine returns the     */
  /* number of factors of "3", kni3, number of factors of "2", kni2.    */
  /* Bails out in case of error.                                        */
  /**********************************************************************/
  /*  Author: D. Majewski, DWD, January 2000                            */
  /**********************************************************************/ 
  /*  Input                                                             */
  /*  kni     INT   number of intervals on a main triangle side         */
  /**********************************************************************/
  /*  Output                                                            */
  /*  kni2    INT      exponent of "2", kni2 > 0                        */ 
  /*  kni3    INT      exponent of "3", either 0 or 1                   */
  /**********************************************************************/

  int mx;

  mx = kni;
  
  *kni2  = 0;
  *kni3  = 0;

  while (mx > 1) { 
      if (mx%2 == 0) {
	*kni2 = *kni2 + 1;
	mx    = mx/2;
      } else if (mx%3 == 0) {
	*kni3 = *kni3 + 1;
	mx    = mx/3;
      } else {
	/* error return */
      }
  }

  /* kni3 must not be greater than */
  
  if (*kni3 > 1) {
    /* error return */
  }
}

/*****************************************************************************/

static
int pow_ii(int x, int n)
{
  int pow;
  
  if (n <= 0) {
    if (n == 0 || x == 1)
      return 1;
    if (x != -1)
      return x == 0 ? 1/x : 0;
    n = -n;
  }
  for (pow = 1; ;) {
    if (n & 01)
      pow *= x;
    if (n >>= 1)
      x *= x;
    else
      break;
  }
  
  return (pow);
}

/*****************************************************************************/

static
struct cart circum_center(struct cart *v0, struct cart *v1, struct cart *v2)
{
  double cnorm;

  struct cart center;
  struct cart e1;
  struct cart e2;
  struct cart cu;

  double *ptmp1;
  double *ptmp2;
  double *ptmp3;

  int j;

  ptmp1 = ((double *)e1.x);
  ptmp2 = ((double *)v1->x);
  ptmp3 = ((double *)v0->x);

  for (j = 0; j < 3; j++) {
    {
      *ptmp1 = *ptmp2 - *ptmp3;
      ptmp1++;
    }
    ptmp3++;
    ptmp2++;
  }
  
  ptmp1 = ((double *)e2.x);
  ptmp2 = ((double *)v2->x);
  ptmp3 = ((double *)v0->x);
  for (j = 0; j < 3; j++) {
    { 
      *ptmp1 = *ptmp2 - *ptmp3;
      ptmp1++;
    }
    ptmp3++;
    ptmp2++;
  }

  cu.x[0] = e1.x[1]*e2.x[2] - e1.x[2]*e2.x[1];
  cu.x[1] = e1.x[2]*e2.x[0] - e1.x[0]*e2.x[2];
  cu.x[2] = e1.x[0]*e2.x[1] - e1.x[1]*e2.x[0];

  cnorm = sqrt(cu.x[0]*cu.x[0] 
	      +cu.x[1]*cu.x[1] 
	      +cu.x[2]*cu.x[2]);

  ptmp1 = ((double *)center.x);
  ptmp2 = ((double *)cu.x);
  for (j = 0; j < 3; j++) {
    {
      *ptmp1 = *ptmp2/cnorm;
      ptmp1++;
    }
    ptmp2++;
  }
  
  return (center);
}

/*****************************************************************************/

struct cart gc2cc(struct geo *position)
{
  double cln;
  double sln;
  double clt;
  double slt;
  
  struct cart x;

  sln = sin(position->lon);
  cln = cos(position->lon);
  slt = sin(position->lat);
  clt = cos(position->lat);

  x.x[0] = cln*clt;
  x.x[1] = sln*clt;
  x.x[2] = slt;

  return (x);
}

/*****************************************************************************/

static
struct geo cc2gc(struct cart *x)
{
  struct geo position;

  double tln;
  double tlt;
  double r;

  if ( !(fabs(x->x[0]) > 0.0) ) {
    if (x->x[1] >= 0.0) {
      position.lon = 0.5*M_PI;
    } else {
      position.lon = 1.5*M_PI;
    }
  } else {
    tln = x->x[1]/x->x[0];
    position.lon = atan(tln);
    if (x->x[0] < 0.0) {
      position.lon = position.lon + M_PI;
    }
    if (position.lon < 0.0) {
      position.lon = position.lon + 2*M_PI;
    }
  }

  r = sqrt(x->x[0]*x->x[0] + x->x[1]*x->x[1]);

  if ( !(fabs(r) > 0.0) ) {
    if (x->x[2] > 0.0) {
      position.lat =  0.5*M_PI;
    } else {
      position.lat = -0.5*M_PI;
    }
  } else {
    tlt = x->x[2]/r;
    position.lat = atan(tlt);
  }

  return (position);
}

/*****************************************************************************/

static
void boundary(struct polygon *poly, 
              int kip1s, int kip1e, int kip2s, int kip2e, int knd)
{
  struct polygon *ptmp1;

  struct cart c;
  struct cart x1;
  struct cart x2;
  struct cart v[6];

  int j1, j2, jd;
  int jm, jm1, jm2;

  struct array polyinfo;

  int tmp1, tmp2, tmp3, tmp4, tmp5;

  polyinfo.rank = 3;
  polyinfo.offset = 0;
  polyinfo.dim[0].lower = kip1s;
  tmp1 = kip1e - polyinfo.dim[0].lower + 1;
  if (tmp1 < 0) tmp1 = 0;
  polyinfo.dim[0].extent = tmp1;
  polyinfo.dim[0].mult = 1;
  polyinfo.offset -= polyinfo.dim[0].lower;
  tmp2 = tmp1;
  polyinfo.dim[1].lower = kip2s;
  tmp1 = kip2e - polyinfo.dim[1].lower + 1;
  if (tmp1 < 0) tmp1 = 0;
  polyinfo.dim[1].extent = tmp1;
  polyinfo.dim[1].mult = tmp2;
  polyinfo.offset -= polyinfo.dim[1].lower*polyinfo.dim[1].mult;
  tmp2 *= tmp1;
  polyinfo.dim[2].lower = 1;
  tmp1 = knd;
  if (tmp1 < 0) tmp1 = 0;
  polyinfo.dim[2].extent = tmp1;
  polyinfo.dim[2].mult = tmp2;
  polyinfo.offset -= polyinfo.dim[2].mult;
  tmp4 = polyinfo.dim[1].mult;
  tmp5 = polyinfo.dim[2].mult;
  tmp3 = polyinfo.offset;

  for (jd = 1; jd <= knd; jd++) {
    for (j2 = kip2s; j2 <= kip2e; j2++) {
      for (j1 = kip1s; j1 <= kip1e; j1++) {
	ptmp1 = &poly[j1+tmp4*j2+tmp5*jd+tmp3];
	c = gc2cc(&ptmp1->center);
	for (jm = 1; jm <= ptmp1->type; jm++) {
	  jm1 = jm;
	  jm2 = (jm % ptmp1->type) + 1;
	  x1 = gc2cc(&ptmp1->boundary[jm1-1]);
	  x2 = gc2cc(&ptmp1->boundary[jm2-1]);
	  if (jd < 6) {
	    v[jm-1] = circum_center(&c,&x1,&x2);
	  } else {
	    v[jm-1] = circum_center(&c,&x2,&x1);
	  }
	}
	if ( jd < 6 )
	  for(jm = 0; jm < ptmp1->type; jm++)
	    ptmp1->boundary[jm] = cc2gc(&v[jm]);
	else
	  for(jm = 0; jm < ptmp1->type; jm++)
	    ptmp1->boundary[ptmp1->type-jm-1] = cc2gc(&v[jm]);
      }
    }
  }

  return;
}

/*****************************************************************************/

static
void neighbours(double *px1, double *px2, 
                int kipx1s, int kipx1e, int kipx2s, int kipx2e, int kndx,
                struct polygon *poly, 
                int kip1s, int kip1e, int kip2s, int kip2e, int knd)
{
  struct polygon *ptmp1;

  int j1, j2, jd, jm, js1, js2;

  struct array px1info, px2info, polyinfo;

  int tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;
  int tmp10, tmp11, tmp12, tmp13, tmp14, tmp15;

  px1info.rank = 3;
  px1info.offset = 0;
  px1info.dim[0].lower = kipx1s;
  tmp1 =  kipx1e - px1info.dim[0].lower + 1;
  if (tmp1 < 0) tmp1 = 0;
  px1info.dim[0].extent = tmp1;
  px1info.dim[0].mult = 1;
  px1info.offset -= px1info.dim[0].lower;
  tmp2 = tmp1;
  px1info.dim[1].lower = kipx2s;
  tmp1 = kipx2e - px1info.dim[1].lower + 1;
  if (tmp1 < 0) tmp1 = 0;
  px1info.dim[1].extent = tmp1;
  px1info.dim[1].mult = tmp2;
  px1info.offset -= px1info.dim[1].lower*px1info.dim[1].mult;
  tmp2 *= tmp1;
  px1info.dim[2].lower = 1;
  tmp1 = kndx;
  if (tmp1 < 0) tmp1 = 0;
  px1info.dim[2].extent = tmp1;
  px1info.dim[2].mult = tmp2;
  px1info.offset -= px1info.dim[2].mult;
  tmp4 = px1info.dim[1].mult;
  tmp5 = px1info.dim[2].mult;
  tmp3 = px1info.offset;

  px2info.rank = 3;
  px2info.offset = 0;
  px2info.dim[0].lower = kipx1s;
  tmp6 = kipx1e - px2info.dim[0].lower + 1;
  if (tmp6 < 0) tmp6 = 0;
  px2info.dim[0].extent = tmp6;
  px2info.dim[0].mult = 1;
  px2info.offset -= px2info.dim[0].lower;
  tmp7 = tmp6;
  px2info.dim[1].lower = kipx2s;
  tmp6 = kipx2e - px2info.dim[1].lower + 1;
  if (tmp6 < 0) tmp6 = 0;
  px2info.dim[1].extent = tmp6;
  px2info.dim[1].mult = tmp7;
  px2info.offset -= px2info.dim[1].lower*px2info.dim[1].mult;
  tmp7 *= tmp6;
  px2info.dim[2].lower = 1;
  tmp6 =  kndx;
  if (tmp6 < 0) tmp6 = 0;
  px2info.dim[2].extent = tmp6;
  px2info.dim[2].mult = tmp7;
  px2info.offset -= px2info.dim[2].mult;
  tmp9 = px2info.dim[1].mult;
  tmp10 = px2info.dim[2].mult;
  tmp8 = px2info.offset;

  polyinfo.rank = 3;
  polyinfo.offset = 0;
  polyinfo.dim[0].lower = kip1s;
  tmp11 = kip1e - polyinfo.dim[0].lower + 1;
  if (tmp11 < 0) tmp11 = 0;
  polyinfo.dim[0].extent = tmp11;
  polyinfo.dim[0].mult = 1;
  polyinfo.offset -= polyinfo.dim[0].lower;
  tmp12 = tmp11;
  polyinfo.dim[1].lower = kip2s;
  tmp11 = kip2e - polyinfo.dim[1].lower + 1;
  if (tmp11 < 0) tmp11 = 0;
  polyinfo.dim[1].extent = tmp11;
  polyinfo.dim[1].mult = tmp12;
  polyinfo.offset -= polyinfo.dim[1].lower*polyinfo.dim[1].mult;
  tmp12 *= tmp11;
  polyinfo.dim[2].lower = 1;
  tmp11 =  knd;
  if (tmp11 < 0) tmp11 = 0;
  polyinfo.dim[2].extent = tmp11;
  polyinfo.dim[2].mult = tmp12;
  polyinfo.offset -= polyinfo.dim[2].mult;
  tmp14 = polyinfo.dim[1].mult;
  tmp15 = polyinfo.dim[2].mult;
  tmp13 = polyinfo.offset;

  for (jd = 1; jd <= kndx; jd++) {
    for (j2 = kipx2s+1; j2 <= kipx2e-1; j2++) {
      for (j1 = kipx1s+1; j1 <= kipx1e-1; j1++) {

	ptmp1 = &poly[j1+tmp14*j2+tmp15*jd+tmp13];

	ptmp1->center.lon = px1[j1+tmp4*j2+tmp5 *jd+tmp3];
	ptmp1->center.lat = px2[j1+tmp9*j2+tmp10*jd+tmp8];

	if (j1 == kipx1s+1 && j2 == kipx2s+1) {

	  ptmp1->type = pentagon;

	  for (jm = 1; jm <= 5; jm++) {

	    if (jd < 6) {
	      ptmp1->boundary[jm-1].lon 
		= px1[kipx1s+1+tmp4*(kipx2s+2)+tmp5 *(jm)+tmp3];
	      ptmp1->boundary[jm-1].lat 
		= px2[kipx1s+1+tmp9*(kipx2s+2)+tmp10*(jm)+tmp8];
	    } else {
	      ptmp1->boundary[jm-1].lon 
		= px1[ kipx1s+1+tmp4*( kipx2s+2)+tmp5*(jm+5)+tmp3];
	      ptmp1->boundary[jm-1].lat 
		= px2[ kipx1s+1+tmp9*( kipx2s+2)+tmp10*(jm+5)+tmp8];
	    }

	  }

	} else if (j1 == kipx1e-1 && j2 == kipx2s+1) {

	  ptmp1->type = pentagon;

	  ptmp1->boundary[0].lon = px1[kipx1e-1+tmp4*(kipx2s+2)+tmp5 *jd+tmp3];
	  ptmp1->boundary[0].lat = px2[kipx1e-1+tmp9*(kipx2s+2)+tmp10*jd+tmp8];
	  ptmp1->boundary[1].lon = px1[kipx1e-2+tmp4*(kipx2s+2)+tmp5 *jd+tmp3];
	  ptmp1->boundary[1].lat = px2[kipx1e-2+tmp9*(kipx2s+2)+tmp10*jd+tmp8];
	  ptmp1->boundary[2].lon = px1[kipx1e-2+tmp4*(kipx2s+1)+tmp5 *jd+tmp3];
	  ptmp1->boundary[2].lat = px2[kipx1e-2+tmp9*(kipx2s+1)+tmp10*jd+tmp8];
	  ptmp1->boundary[3].lon = px1[kipx1e-1+tmp4*(kipx2s)  +tmp5 *jd+tmp3];
	  ptmp1->boundary[3].lat = px2[kipx1e-1+tmp9*(kipx2s)  +tmp10*jd+tmp8];
	  ptmp1->boundary[4].lon = px1[kipx1e  +tmp4*(kipx2s+1)+tmp5 *jd+tmp3];
	  ptmp1->boundary[4].lat = px2[kipx1e  +tmp9*(kipx2s+1)+tmp10*jd+tmp8];

	} else if (j1 == kipx1s+1 && j2 == kipx2e-1) {

	  ptmp1->type = pentagon;

	  ptmp1->boundary[0].lon = px1[kipx1s+2+tmp4*(kipx2e-2)+tmp5 *jd+tmp3];
	  ptmp1->boundary[0].lat = px2[kipx1s+2+tmp9*(kipx2e-2)+tmp10*jd+tmp8];
	  ptmp1->boundary[1].lon = px1[kipx1s+2+tmp4*(kipx2e-1)+tmp5 *jd+tmp3];
	  ptmp1->boundary[1].lat = px2[kipx1s+2+tmp9*(kipx2e-1)+tmp10*jd+tmp8];
	  ptmp1->boundary[2].lon = px1[kipx1s+1+tmp4*(kipx2e)  +tmp5 *jd+tmp3];
	  ptmp1->boundary[2].lat = px2[kipx1s+1+tmp9*(kipx2e)  +tmp10*jd+tmp8];
	  ptmp1->boundary[3].lon = px1[kipx1s  +tmp4*(kipx2e-1)+tmp5 *jd+tmp3];
	  ptmp1->boundary[3].lat = px2[kipx1s  +tmp9*(kipx2e-1)+tmp10*jd+tmp8];
	  ptmp1->boundary[4].lon = px1[kipx1s+1+tmp4*(kipx2e-2)+tmp5 *jd+tmp3];
	  ptmp1->boundary[4].lat = px2[kipx1s+1+tmp9*(kipx2e-2)+tmp10*jd+tmp8];

	} else if (j1 == kipx1e-1 && j2 == kipx2e-1) {

	  ptmp1->type = pentagon;

	  ptmp1->boundary[0].lon = px1[kipx1e  +tmp4*(kipx2e)  +tmp5 *jd+tmp3];
	  ptmp1->boundary[0].lat = px2[kipx1e  +tmp9*(kipx2e)  +tmp10*jd+tmp8];
	  ptmp1->boundary[1].lon = px1[kipx1e-2+tmp4*(kipx2e)  +tmp5 *jd+tmp3];
	  ptmp1->boundary[1].lat = px2[kipx1e-2+tmp9*(kipx2e)  +tmp10*jd+tmp8];
	  ptmp1->boundary[2].lon = px1[kipx1e-2+tmp4*(kipx2e-1)+tmp5 *jd+tmp3];
	  ptmp1->boundary[2].lat = px2[kipx1e-2+tmp9*(kipx2e-1)+tmp10*jd+tmp8];
	  ptmp1->boundary[3].lon = px1[kipx1e-1+tmp4*(kipx2e-2)+tmp5 *jd+tmp3];
	  ptmp1->boundary[3].lat = px2[kipx1e-1+tmp9*(kipx2e-2)+tmp10*jd+tmp8];
	  ptmp1->boundary[4].lon = px1[kipx1e  +tmp4*(kipx2e-2)+tmp5 *jd+tmp3];
	  ptmp1->boundary[4].lat = px2[kipx1e  +tmp9*(kipx2e-2)+tmp10*jd+tmp8];

	} else {

	  for(jm = 1; jm <= 6; jm++) {

	    ptmp1->type = hexagon;

	    js1 = j1+ispokes[jm-1];
	    js2 = j2+ispokes[jm+5];

	    ptmp1->boundary[jm-1].lon = px1[js1+tmp4*js2+tmp5 *jd+tmp3];
	    ptmp1->boundary[jm-1].lat = px2[js1+tmp9*js2+tmp10*jd+tmp8];

	  }

	}

      }

    }

  }

  return;
}

/*****************************************************************************/

static
void xd(double *p, int kip1s, int kip1e, int kip2s, int kip2e, int knd,
	double *px, int kipx1s, int kipx1e, int kipx2s, int kipx2e, int kndx)
{
  int mi1sm1, mi1ep1, mi2sm1, mi2ep1;
  int mns, mpe, mpw, maw, mae, mpp;
  int j, j1, j2, jd;

  struct array pinfo, pxinfo;

  int tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9, tmp10;

  pinfo.rank = 3;
  pinfo.offset = 0;
  pinfo.dim[0].lower = kip1s;
  tmp1 =  kip1e - pinfo.dim[0].lower+1;
  if (tmp1 < 0) tmp1 = 0;
  pinfo.dim[0].extent = tmp1;
  pinfo.dim[0].mult = 1;
  pinfo.offset -= pinfo.dim[0].lower;
  tmp2 = tmp1;
  pinfo.dim[1].lower = kip2s;
  tmp1 =  kip2e - pinfo.dim[1].lower+1;
  if (tmp1 < 0) tmp1 = 0;
  pinfo.dim[1].extent = tmp1;
  pinfo.dim[1].mult = tmp2;
  pinfo.offset -= pinfo.dim[1].lower*pinfo.dim[1].mult;
  tmp2 *= tmp1;
  pinfo.dim[2].lower = 1;
  tmp1 =  knd;
  if (tmp1 < 0) tmp1 = 0;
  pinfo.dim[2].extent = tmp1;
  pinfo.dim[2].mult = tmp2;
  pinfo.offset -= pinfo.dim[2].mult;
  tmp4 = pinfo.dim[1].mult;
  tmp5 = pinfo.dim[2].mult;
  tmp3 = pinfo.offset;

  pxinfo.rank = 3;
  pxinfo.offset = 0;
  pxinfo.dim[0].lower = kipx1s;
  tmp6 =  kipx1e - pxinfo.dim[0].lower+1;
  if (tmp6 < 0) tmp6 = 0;
  pxinfo.dim[0].extent = tmp6;
  pxinfo.dim[0].mult = 1;
  pxinfo.offset -= pxinfo.dim[0].lower;
  tmp7 = tmp6;
  pxinfo.dim[1].lower = kipx2s;
  tmp6 =  kipx2e - pxinfo.dim[1].lower+1;
  if (tmp6 < 0) tmp6 = 0;
  pxinfo.dim[1].extent = tmp6;
  pxinfo.dim[1].mult = tmp7;
  pxinfo.offset -= pxinfo.dim[1].lower*pxinfo.dim[1].mult;
  tmp7 *= tmp6;
  pxinfo.dim[2].lower = 1;
  tmp6 =  kndx;
  if (tmp6 < 0) tmp6 = 0;
  pxinfo.dim[2].extent = tmp6;
  pxinfo.dim[2].mult = tmp7;
  pxinfo.offset -= pxinfo.dim[2].mult;
  tmp9 = pxinfo.dim[1].mult;
  tmp10 = pxinfo.dim[2].mult;
  tmp8 = pxinfo.offset;
  tmp1 = pxinfo.dim[0].extent;
  tmp2 = pxinfo.dim[1].extent;
  tmp6 = pxinfo.dim[2].extent;
  
  for(jd = 1;jd <= knd;jd++) {
    for(j2 = kip2s;j2 <= kip2e;j2++) {
      for(j1 = kip1s;j1 <= kip1e; j1++) {
	px[j1+tmp9*j2+tmp10*jd+tmp8] = p[j1+tmp4*j2+tmp5*jd+tmp3];
      }
    }
  }

  mi1sm1 =  kip1s-1;
  mi1ep1 =  kip1e+1;
  mi2sm1 =  kip2s-1;
  mi2ep1 =  kip2e+1;

  for(jd = 1; jd <= knd; jd++) {
    mns = (jd-1)/5;
    mpe = jd+1-(jd/(5*(1+mns)))*5;
    mpw = jd-1+((mns*10+6-jd)/(5*(1+mns)))*5;
    mae = jd+5-9*mns-5*(jd/10);
    maw = jd+4+((6-jd)/5)*5-9*mns;
    mpp = jd+3-((jd+2)/5)*5+5*mns;
    for(j = kip2s; j <= kip1e; j++) {
      px[j     +tmp9*mi2sm1+tmp10*jd+tmp8] 
	= p[kip1s+1  +tmp4*j          +tmp5*mpw+tmp3];
      px[mi1sm1+tmp9*(j+1) +tmp10*jd+tmp8] 
	= p[j-1      +tmp4*(kip2s+1)  +tmp5*mpe+tmp3];
      px[mi1ep1+tmp9*j     +tmp10*jd+tmp8] 
	= p[kip1e+1-j+tmp4*(kip2e-1)  +tmp5*maw+tmp3];
      px[j-1   +tmp9*mi2ep1+tmp10*jd+tmp8] 
	= p[kip1e-1  +tmp4*(kip2e+1-j)+tmp5*mae+tmp3];
    }
    px[mi1sm1+tmp9*kip2s +tmp10*jd+tmp8] = p[kip1s+1+tmp4*kip2s+tmp5*mpp+tmp3];
    px[kip1s +tmp9*mi2sm1+tmp10*jd+tmp8] = p[kip1s+1+tmp4*kip2s+tmp5*mpp+tmp3];

    px[mi1ep1+tmp9*mi2sm1+tmp10*jd+tmp8] = px[kip1e +tmp9*mi2sm1+tmp10*jd+tmp8];
    px[mi1sm1+tmp9*mi2ep1+tmp10*jd+tmp8] = px[mi1sm1+tmp9*kip2e +tmp10*jd+tmp8];

    px[mi1ep1+tmp9*kip2e +tmp10*jd+tmp8] = p[kip1e-1+tmp4*kip2s+tmp5*mae+tmp3];
    px[kip1e +tmp9*mi2ep1+tmp10*jd+tmp8] = p[kip1e-1+tmp4*kip2s+tmp5*mae+tmp3];

    px[mi1sm1+tmp9*mi2sm1+tmp10*jd+tmp8] = px[kip1s+tmp9*mi2sm1+tmp10*jd+tmp8];
    px[mi1ep1+tmp9*mi2ep1+tmp10*jd+tmp8] = px[kip1e+tmp9*mi2ep1+tmp10*jd+tmp8];
  }

  return;
}

/*****************************************************************************/

static
void tricntr(double *pxn, 
             int kig1s, int kig1e, int kig2s, int kig2e, int knd, 
             int kjd, int kni)
{
  int id1, id2, id3, ioffset;
  double r1, r2, r3;

  int j, mi1, mi2;
  double zxnorm;

  id1 = kig1e-kig1s+1;
  id2 = id1*(kig2e-kig2s+1);
  id3 = id2*3;
  ioffset = -(id1+id2+id3);
    
  for (j = 1; j <= 2; ++j) 
    {
      mi1 = j*kni/3;
      mi2 = (j-1)*kni+1;
      pxn[mi1+id1*(mi1+1)+id2*1+id3*kjd+ioffset] 
        = pxn[mi2-1+id1*(mi2)+id2*1+id3*kjd+ioffset] 
	+pxn[kni+id1*1+id2*1+id3*kjd+ioffset] 
	+pxn[0+id1*(kni+1)+id2*1+id3*kjd+ioffset];
      pxn[mi1+id1*(mi1+1)+id2*2+id3*kjd+ioffset] 
        = pxn[mi2-1+id1*(mi2)+id2*2+id3*kjd+ioffset] 
	+pxn[kni+id1*1+id2*2+id3*kjd+ioffset] 
	+pxn[0+id1*(kni+1)+id2*2+id3*kjd+ioffset];
      pxn[mi1+id1*(mi1+1)+id2*3+id3*kjd+ioffset] 
        = pxn[mi2-1+id1*(mi2)+id2*3+id3*kjd+ioffset] 
	+pxn[kni+id1*1+id2*3+id3*kjd+ioffset] 
	+pxn[0+id1*(kni+1)+id2*3+id3*kjd+ioffset];
      /* Normalize to unit-sphere */
      
      r1 = pxn[mi1+id1*(mi1+1)+id2*1+id3*kjd+ioffset];
      r2 = pxn[mi1+id1*(mi1+1)+id2*2+id3*kjd+ioffset];
      r3 = pxn[mi1+id1*(mi1+1)+id2*3+id3*kjd+ioffset];
      zxnorm = 1.0/sqrt(r1*r1+r2*r2+r3*r3);
      
      pxn[mi1+id1*(mi1+1)+id2*1+id3*kjd+ioffset] 
        = zxnorm*pxn[mi1+id1*(mi1+1)+id2*1+id3*kjd+ioffset];
      pxn[mi1+id1*(mi1+1)+id2*2+id3*kjd+ioffset] 
        = zxnorm*pxn[mi1+id1*(mi1+1)+id2*2+id3*kjd+ioffset];
      pxn[mi1+id1*(mi1+1)+id2*3+id3*kjd+ioffset] 
        = zxnorm*pxn[mi1+id1*(mi1+1)+id2*3+id3*kjd+ioffset];
    }
  
  return;
} /* tricntr */

/****************************************************************************/

static
void gcpt(double *pxn, 
          int kig1s, int kig1e, int kig2s, int kig2e, int knd, 
          int kjd, double pgamma, 
          int ki1, int kj1, int ki2, int kj2, int ki, int kj)
{
  int id1, id2, id3, ioffset;
  double r1, r2, r3;
    
  double zbeta, zalpha, zchord, ztheta;
    
  /* Calculate "zchord", the Cartesian distance between x1 and x2 */

  id1 = kig1e-kig1s+1;
  id2 = id1*(kig2e-kig2s+1);
  id3 = id2*3;
  ioffset = -(id1+id2+id3);

  r1 = (pxn[ki2+id1*kj2+id2*1+id3*kjd+ioffset] 
	-pxn[ki1+id1*kj1+id2*1+id3*kjd+ioffset]);
  r2 = (pxn[ki2+id1*kj2+id2*2+id3*kjd+ioffset] 
	-pxn[ki1+id1*kj1+id2*2+id3*kjd+ioffset]);
  r3 = (pxn[ki2+id1*kj2+id2*3+id3*kjd+ioffset] 
	-pxn[ki1+id1*kj1+id2*3+id3*kjd+ioffset]);
  
  zchord = sqrt((r1*r1)+(r2*r2)+(r3*r3));
    
  /* Calculate "ztheta", the great circle angle between x1 and x2 */
  ztheta = 2.0*asin(zchord*0.5);
    
  /* Calculate the weighting factors which follow from the condition */
  /* that x is a point on the unit-sphere, too. */
  zbeta = sin(pgamma*ztheta)/sin(ztheta);
  zalpha = sin((1.0-pgamma)*ztheta)/sin(ztheta);
    
  /* Store the (x,y,z) coordinates of the point x into the array pxn */
  pxn[ki+id1*kj+id2*1+id3*kjd+ioffset] 
    = zalpha*pxn[ki1+id1*kj1+id2*1+id3*kjd+ioffset] 
    +zbeta *pxn[ki2+id1*kj2+id2*1+id3*kjd+ioffset];
  pxn[ki+id1*kj+id2*2+id3*kjd+ioffset] 
    = zalpha*pxn[ki1+id1*kj1+id2*2+id3*kjd+ioffset] 
    +zbeta *pxn[ki2+id1*kj2+id2*2+id3*kjd+ioffset];
  pxn[ki+id1*kj+id2*3+id3*kjd+ioffset] 
    = zalpha*pxn[ki1+id1*kj1+id2*3+id3*kjd+ioffset] 
    +zbeta *pxn[ki2+id1*kj2+id2*3+id3*kjd+ioffset];

  return;
} /* gcpt */

/****************************************************************************/

static
void glo_coor(double *pxn, double *prlon, double *prlat, 
              int kig1s, int kig1e, int kig2s, int kig2e, int knd, 
              int kni2, int kni3)
{
  int id1, id2, id3, ioffset, joffset;

  double zsgn;
  int j1, j2;
  double zrlon, zcosw, zsinw;
  int jb, jd, ml, mm;
  double zgamma, zw;
  int mi1, mi2, ml2, ml3, mni;

  int *mcosv;

  /*
   * Calculate the Cartesian coordinates of the gridpoints of the
   * icosahedral grid on the unit sphere.  The grid resolution   
   * corresponds to a subdivision of the edges of the original   
   * icosahedral triangles into mni equal parts.                 
   */  

  if ( (mcosv = (int*) malloc(knd*sizeof(int))) == NULL ) 
    {
      perror("malloc mcosv");
      exit (-1);
    }

  id1 = kig1e-kig1s+1;
  id2 = id1*(kig2e-kig2s+1);
  id3 = id2*3;
  ioffset = -(id1+id2+id3);
  joffset = -(id1+id2);

  /* Compute angles associated with the icosahedron. */

  zw = acos(1.0/(sin(pid5)*2.0))*2.0;
  zcosw = cos(zw);
  zsinw = sin(zw);
  mni = pow_ii(2, kni2)*pow_ii(3, kni3);
  
  /*     Compute the local array mcosv, i.e. the meridian angle locations */

  for (jd = 1; jd <= knd; ++jd) 
    {
      if (jd % 2 == 1) 
        {
          mcosv[(jd+1)/2-1] = jd-2-knd*((jd-1)/7);
        } else {
          mcosv[jd/2+4]     = jd-2-knd*((jd-1)/7);
        }
    }

  /**************************************************************************/
  /*     Loop over the ten diamonds computing diamond vertices (x,y,z)      */
  /*     coordinates and then iteratively bisecting them kni2 times.        */
  /*     First a trisection is performed, if required (kni3=1).             */

  for (jd = 1; jd <= knd; ++jd) 
    {
      
      /*     Toggle the hemisphere */
      if (jd >= 6) 
        {
          zsgn = -1.0;
        } else {
          zsgn =  1.0;
        }
      
      /*     Compute the meridian angle for each diamond "home" vertex. */
      zrlon = mcosv[jd-1]*pid5;
      
      /*     Every diamond has one vertex at a pole (N or S). */
      /*     Label this point (0,1,,) in each diamond, and */
      /*     initialize it to have the (x,y,z) coordinates of */
      /*     the pole point on the unit sphere. */
      pxn[0+id1*1+id2*1+id3*jd+ioffset] = 0.0;
      pxn[0+id1*1+id2*2+id3*jd+ioffset] = 0.0;
      pxn[0+id1*1+id2*3+id3*jd+ioffset] = zsgn;
      
      /*     Now initialize the (x,y,z) coordinates of the "home" vertex, */
      /*     which defines which diamond we are talking about, and label  */
      /*     this point (mni,1,,). */
      pxn[mni+id1*1+id2*1+id3*jd+ioffset] = zsinw*cos(zrlon);
      pxn[mni+id1*1+id2*2+id3*jd+ioffset] = zsinw*sin(zrlon);
      pxn[mni+id1*1+id2*3+id3*jd+ioffset] = zcosw*zsgn;
      
      /*     Next initialize the (x,y,z) coordinates for the corner of the */
      /*     diamond on the same latitude as the (mni,1,,) vertex, which   */
      /*     is (0,mni+1,,) */
      pxn[0+id1*(mni+1)+id2*1+id3*jd+ioffset] = zsinw*cos(zrlon+2*pid5);
      pxn[0+id1*(mni+1)+id2*2+id3*jd+ioffset] = zsinw*sin(zrlon+2*pid5);
      pxn[0+id1*(mni+1)+id2*3+id3*jd+ioffset] = zcosw*zsgn;
      
      /*     Initialize the last diamond vertex, which is located */
      /*     in the opposite hemisphere as (mni,mni+1,,)          */
      pxn[mni+id1*(mni+1)+id2*1+id3*jd+ioffset] = zsinw*cos(zrlon+pid5);
      pxn[mni+id1*(mni+1)+id2*2+id3*jd+ioffset] = zsinw*sin(zrlon+pid5);
      pxn[mni+id1*(mni+1)+id2*3+id3*jd+ioffset] = -(zcosw*zsgn);
    
      /***********************************************************************/
      /*     First a trisection is performed, if required (kni3=1).          */
      
      if (kni3 == 1) 
        {
          ml3 = mni/3;
          
          /*     Trisect the rows of the diamond. */
          for (j1 = 1; j1 <= 2; ++j1) 
            {
              for (j2 = 1; j2 <= 2; ++j2) 
                {
                  mi1 = (j1-1)*mni;
                  mi2 = j2*ml3+1;
                  zgamma = (double) j2/3.0;
                  gcpt(pxn, kig1s, kig1e, kig2s, kig2e, knd, jd, 
                       zgamma, mi1, 1, mi1, mni+1, mi1, mi2);
                }
            }
          
          /*     Trisect the columns of the diamond. */
          for (j1 = 1; j1 <= 2; ++j1) 
            {
              for (j2 = 1; j2 <= 2; ++j2) 
                {
                  mi1 = j2*ml3;
                  mi2 = (j1-1)*mni+1;
                  zgamma = (double) j2/3.0;
                  gcpt(pxn, kig1s, kig1e, kig2s, kig2e, knd, jd,
                       zgamma, 0, mi2, mni, mi2, mi1, mi2);
                }
            }
      
          /*     Trisect the diagonal of the diamond. */
          for (j2 = 1; j2 <= 2; ++j2) 
            {
              mi1 = mni-j2*ml3;
              mi2 = j2*ml3+1;
              zgamma = (double) j2/(float)3.;
              gcpt(pxn, kig1s, kig1e, kig2s, kig2e, knd, jd, 
                   zgamma, mni, 1, 0, mni+1, mi1, mi2);
            }
          
          /*     Compute coordinates of icosahedral triangle centers. */
          tricntr(pxn, kig1s, kig1e, kig2s, kig2e, knd, jd, mni);
          
        }
      
      /***********************************************************************/
      /*     Find the coordinates of the triangle nodes by iteratively       */
      /*     bisecting the diamond intervals.                                */
      
      for (jb = 0; jb <= kni2-1; ++jb) 
        {
          mm = pow_ii(3, kni3)*pow_ii(2, jb);
          ml = mni/mm;
          ml2 = ml/2;
          
          /*     Compute the rows of the diamond. */
          
          for (j1 = 1; j1 <= mm+1; ++j1) 
            {
              for (j2 = 1; j2 <= mm; ++j2) 
                {
                  mi1 = (j1-1)*ml;
                  mi2 = (j2-1)*ml+ml2+1;
                  gcpt(pxn, kig1s, kig1e, kig2s, kig2e, knd, jd, 
                       0.5, mi1, mi2-ml2, mi1, mi2+ml2, mi1, mi2);
                }
            }
          
          /*     Compute the columns of diamond. */
          
          for (j1 = 1; j1 <= mm+1; ++j1) 
            {
              for (j2 = 1; j2 <= mm; ++j2) 
                {
                  mi1 = (j2-1)*ml+ml2;
                  mi2 = (j1-1)*ml+1;
                  gcpt(pxn, kig1s, kig1e, kig2s, kig2e, knd, jd,
                       0.5, mi1-ml2, mi2, mi1+ml2, mi2, mi1, mi2);
                }
            }
          
          /*     Compute the diagonals of the diamond. */
          
          for (j1 = 1; j1 <= mm; ++j1) 
            {
              for (j2 = 1; j2 <= mm; ++j2) 
                {
                  mi1 = (j1-1)*ml+ml2;
                  mi2 = (j2-1)*ml+ml2+1;
                  gcpt(pxn, kig1s, kig1e, kig2s, kig2e, knd, jd,
                       0.5, mi1-ml2, mi2+ml2, mi1+ml2, mi2-ml2, mi1, mi2);
                }
            }
          
        }
      
      /***********************************************************************/
      /* Set pxn to 0 if it is less than 2.5 e-14 to avoid round-off errors  */
      
      for (j2 = kig2s; j2 <= kig2e; ++j2) 
        {
          for (j1 = kig1s; j1 <= kig1e; ++j1) 
            {
              if (fabs(pxn[j1+id1*j2+id2*1+id3*jd+ioffset]) < 2.5e-14) 
                {
                  pxn[j1+id1*j2+id2*1+id3*jd+ioffset] = 0.0;
                }
              if (fabs(pxn[j1+id1*j2+id2*2+id3*jd+ioffset]) < 2.5e-14) 
                {
                  pxn[j1+id1*j2+id2*2+id3*jd+ioffset] = 0.0;
                }
              if (fabs(pxn[j1+id1*j2+id2*3+id3*jd+ioffset]) < 2.5e-14) 
                {
                  pxn[j1+id1*j2+id2*3+id3*jd+ioffset] = 0.0;
                }
            }
        }
    }
  
  /*************************************************************************/
  /*     Calculate the longitude "prlon" and the latitude "prlat";         */
  /*     only for the core of the diamonds, not the extended ones.         */
  
  for (jd = 1; jd <= knd; ++jd) 
    {
      for (j2 = kig2s; j2 <= kig2e; ++j2) 
        {
          for (j1 = kig1s; j1 <= kig1e; ++j1) 
            {
              prlon[j1+id1*j2+id2*jd+joffset] 
                = atan2(pxn[j1+id1*j2+id2*2+id3*jd+ioffset],
                        pxn[j1+id1*j2+id2*1+id3*jd+ioffset]+1.0e-20);
              prlat[j1+id1*j2+id2*jd+joffset] 
                = asin(pxn[j1+id1*j2+id2*3+id3*jd+ioffset]);
            }
        }
    }
  
  free(mcosv);
  
  return;
} /* glo_coor */

/*****************************************************************************/

static
void initmask(int *mask, int ni, int nd)
{
  int jd;

  int tmp1, tmp2, tmp3, tmp4, tmp5, tmp6, tmp7, tmp8, tmp9;

  int *ptmp1;
  char *ptmp2;

  struct array section;
  struct array maskinfo;

  maskinfo.rank = 3;
  maskinfo.offset = 0;
  maskinfo.dim[0].lower = 0;
  tmp1 =  ni+1;
  if (tmp1 < 0) tmp1 = 0;
  maskinfo.dim[0].extent = tmp1;
  maskinfo.dim[0].mult = 1;
  maskinfo.offset -= 0;
  tmp2 = tmp1;
  maskinfo.dim[1].lower = 1;
  tmp1 =  ni+1;
  if (tmp1 < 0) tmp1 = 0;
  maskinfo.dim[1].extent = tmp1;
  maskinfo.dim[1].mult = tmp2;
  maskinfo.offset -= maskinfo.dim[1].mult;
  tmp2 *= tmp1;
  maskinfo.dim[2].lower = 1;
  tmp1 =  nd;
  if (tmp1 < 0) tmp1 = 0;
  maskinfo.dim[2].extent = tmp1;
  maskinfo.dim[2].mult = tmp2;
  maskinfo.offset -= maskinfo.dim[2].mult;
  tmp4 = maskinfo.dim[1].mult;
  tmp5 = maskinfo.dim[2].mult;
  tmp3 = maskinfo.offset;
  tmp1 = maskinfo.dim[0].extent;
  tmp2 = maskinfo.dim[1].extent;
  tmp9 = maskinfo.dim[2].extent;

  ptmp1 = mask;
  for (tmp8 = 0; tmp8 < tmp9; tmp8++)  {
    for (tmp7 = 0; tmp7 < tmp2; tmp7++)  {
      for (tmp6 = 0; tmp6 < tmp1; tmp6++)  
	*ptmp1++ = 1;
    }
  }

  for (jd = 1; jd <= 10; jd++) {
    switch (jd) {
    case 1:
      break;
    case 3:
    case 4:
    case 2:
      tmp1 = 1;
      ptmp1 = mask;
      tmp2 = maskinfo.dim[0].extent;
      if (tmp2 < 0) tmp2 = 0;
      tmp1 *= maskinfo.dim[0].extent;
      tmp1 *= maskinfo.dim[1].extent;
      ptmp1 += tmp1*(jd-1);
      for (tmp6 = 0; tmp6 < tmp2; tmp6++)  
	*ptmp1++ = 0;
      break;
    case 5:
      tmp1 = 1;
      ptmp1 = mask;
      tmp2 = maskinfo.dim[0].extent;
      if (tmp2 < 0) tmp2 = 0;
      tmp1 *= maskinfo.dim[0].extent;
      tmp1 *= maskinfo.dim[1].extent;
      ptmp1 += tmp1*(jd-1);
      for (tmp6 = 0; tmp6 < tmp2; tmp6++)  
	*ptmp1++ = 0;

      tmp1 = 4;
      section.addr = (char *)mask;
      section.offset = 0;
      tmp1 *= maskinfo.dim[0].extent;
      section.dim[0].mult = tmp1;
      section.dim[0].extent = maskinfo.dim[1].extent;
      if (section.dim[0].extent < 0) section.dim[0].extent = 0;
      section.offset -= section.dim[0].mult;
      section.dim[0].lower = 1;
      tmp1 *= maskinfo.dim[1].extent;
      section.addr += tmp1*(jd-1);
      tmp2 = section.dim[0].extent;
      ptmp2 = section.addr;
      for (tmp6 = 0; tmp6 < tmp2; tmp6++)  {
	*((int *)ptmp2) = 0;
	ptmp2 += section.dim[0].mult;
      }
      break;
    case 6:

      tmp1 = 1;
      ptmp1 = mask;
      tmp2 = maskinfo.dim[0].extent;
      if (tmp2 < 0) tmp2 = 0;
      tmp1 *= maskinfo.dim[0].extent;
      ptmp1 += tmp1*ni;
      tmp1 *= maskinfo.dim[1].extent;
      ptmp1 += tmp1*(jd-1);
      for (tmp6 = 0; tmp6 < tmp2; tmp6++)  
	*ptmp1++ = 0;

      tmp1 = 4;
      section.addr = (char *)mask;
      section.offset = 0;
      section.addr += tmp1*ni;
      tmp1 *= maskinfo.dim[0].extent;
      section.dim[0].mult = tmp1;
      section.dim[0].extent = maskinfo.dim[1].extent;
      if (section.dim[0].extent < 0) section.dim[0].extent = 0;
      section.offset -= section.dim[0].mult;
      section.dim[0].lower = 1;
      tmp1 *= maskinfo.dim[1].extent;
      section.addr += tmp1*(jd-1);
      tmp2 = section.dim[0].extent;
      ptmp2 = section.addr;
      for (tmp6 = 0; tmp6 < tmp2; tmp6++)  {
	*((int *)ptmp2) = 0;
	ptmp2 += section.dim[0].mult;
      }
      break;
    case 8:
    case 9:
    case 7:
      tmp1 = 1;
      ptmp1 = mask;
      tmp2 = maskinfo.dim[0].extent;
      if (tmp2 < 0) tmp2 = 0;
      tmp1 *= maskinfo.dim[0].extent;
      tmp1 *= maskinfo.dim[1].extent;
      ptmp1 += tmp1*(jd-1);
      for (tmp6 = 0; tmp6 < tmp2; tmp6++)  
	*ptmp1++ = 0;

      tmp1 = 1;
      ptmp1 = mask;
      tmp2 = maskinfo.dim[0].extent;
      if (tmp2 < 0) tmp2 = 0;
      tmp1 *= maskinfo.dim[0].extent;
      ptmp1 += tmp1*ni;
      tmp1 *= maskinfo.dim[1].extent;
      ptmp1 += tmp1*(jd-1);
      for (tmp6 = 0; tmp6 < tmp2; tmp6++)  
	*ptmp1++ = 0;

      tmp1 = 4;
      section.addr = (char *)mask;
      section.offset = 0;
      section.addr += tmp1*ni;
      tmp1 *= maskinfo.dim[0].extent;
      section.dim[0].mult = tmp1;
      section.dim[0].extent = maskinfo.dim[1].extent;
      if (section.dim[0].extent < 0) section.dim[0].extent = 0;
      section.offset -= section.dim[0].mult;
      section.dim[0].lower = 1;
      tmp1 *= maskinfo.dim[1].extent;
      section.addr += tmp1*(jd-1);
      tmp2 = section.dim[0].extent;
      ptmp2 = section.addr;
      for (tmp6 = 0; tmp6 < tmp2; tmp6++)  {
	*((int *)ptmp2) = 0;
	ptmp2 += section.dim[0].mult;
      }
      break;
    case 10:

      tmp1 = 1;
      ptmp1 = mask;
      tmp2 = maskinfo.dim[0].extent;
      if (tmp2 < 0) tmp2 = 0;
      tmp1 *= maskinfo.dim[0].extent;
      tmp1 *= maskinfo.dim[1].extent;
      ptmp1 += tmp1*(jd-1);
      for (tmp6 = 0; tmp6 < tmp2; tmp6++)  
	*ptmp1++ = 0;

      tmp1 = 1;
      ptmp1 = mask;
      tmp2 = maskinfo.dim[0].extent;
      if (tmp2 < 0) tmp2 = 0;
      tmp1 *= maskinfo.dim[0].extent;
      ptmp1 += tmp1*ni;
      tmp1 *= maskinfo.dim[1].extent;
      ptmp1 += tmp1*(jd-1);
      for (tmp6 = 0; tmp6 < tmp2; tmp6++)  
	*ptmp1++ = 0;

      tmp1 = 4;
      section.addr = (char *)mask;
      section.offset = 0;
      tmp1 *= maskinfo.dim[0].extent;
      section.dim[0].mult = tmp1;
      section.dim[0].extent = maskinfo.dim[1].extent;
      if (section.dim[0].extent < 0) section.dim[0].extent = 0;
      section.offset -= section.dim[0].mult;
      section.dim[0].lower = 1;
      tmp1 *= maskinfo.dim[1].extent;
      section.addr += tmp1*(jd-1);
      tmp2 = section.dim[0].extent;
      ptmp2 = section.addr;
      for (tmp6 = 0; tmp6 < tmp2; tmp6++)  {
	*((int *)ptmp2) = 0;
	ptmp2 += section.dim[0].mult;
      }

      tmp1 = 4;
      section.addr = (char *)mask;
      section.offset = 0;
      section.addr += tmp1*ni;
      tmp1 *= maskinfo.dim[0].extent;
      section.dim[0].mult = tmp1;
      section.dim[0].extent = maskinfo.dim[1].extent;
      if (section.dim[0].extent < 0) section.dim[0].extent = 0;
      section.offset -= section.dim[0].mult;
      section.dim[0].lower = 1;
      tmp1 *= maskinfo.dim[1].extent;
      section.addr += tmp1*(jd-1);
      tmp2 = section.dim[0].extent;
      ptmp2 = section.addr;
      for (tmp6 = 0; tmp6 < tmp2; tmp6++)  {
	*((int *)ptmp2) = 0;
	ptmp2 += section.dim[0].mult;
      }
      break;
    }
  }

  return;
}

/*****************************************************************************/

void gme_grid_restore(double *p, int ni, int nd)
{
  int j, jd;
    
  int tmp1, tmp2, tmp3, tmp4, tmp5;

  struct array pinfo;

  pinfo.rank = 3;
  pinfo.offset = 0;
  pinfo.dim[0].lower = 0;
  tmp1 = ni+1;
  if (tmp1 < 0) tmp1 = 0;
  pinfo.dim[0].extent = tmp1;
  pinfo.dim[0].mult = 1;
  pinfo.offset -= 0;
  tmp2 = tmp1;
  pinfo.dim[1].lower = 1;
  tmp1 = ni+1;
  if (tmp1 < 0) tmp1 = 0;
  pinfo.dim[1].extent = tmp1;
  pinfo.dim[1].mult = tmp2;
  pinfo.offset -= pinfo.dim[1].mult;
  tmp2 *= tmp1;
  pinfo.dim[2].lower = 1;
  tmp1 =  nd;
  if (tmp1 < 0) tmp1 = 0;
  pinfo.dim[2].extent = tmp1;
  pinfo.dim[2].mult = tmp2;
  pinfo.offset -= pinfo.dim[2].mult;
  tmp4 = pinfo.dim[1].mult;
  tmp5 = pinfo.dim[2].mult;
  tmp3 = pinfo.offset;

  for (jd = 1; jd <= 10; jd++) {
    switch (jd) {
    case 1:
      break;
    case 3:
    case 4:
    case 2:
      for (j = 0; j <= ni; j++) {
	p[j+tmp4+tmp5*jd+tmp3] = p[tmp4*(j+1)+tmp5*(jd-1)+tmp3];
      }
      break;
    case 5:
      for (j = 0; j <= ni; j++) {
	p[j+tmp4+tmp5*jd+tmp3] = p[tmp4*(j+1)+tmp5*(jd-1)+tmp3];
      }
      for (j = 0; j <= ni; j++) {
	p[tmp4*(j+1)+tmp5*5+tmp3] = p[j+tmp4+tmp5+tmp3];
      }
      break;
    case 6:
      for (j = 0; j <= ni; j++) {
	p[j+tmp4*(ni+1)+tmp5*6+tmp3] = p[ni+tmp4*(ni+1-j)+tmp5*2+tmp3];
      }
      for (j = 0; j <= ni; j++) {
	p[ni+tmp4*(j+1)+tmp5*6+tmp3] = p[ni-j+tmp4*(ni+1)+tmp5+tmp3];
      }
      break;
    case 8:
    case 9:
    case 7:
      for (j = 0; j <= ni; j++) {
	p[j+tmp4*(ni+1)+tmp5*jd+tmp3] = p[ni+tmp4*(ni+1-j)+tmp5*(jd-4)+tmp3];
      }
      for (j = 0; j <= ni; j++) {
	p[ni+tmp4*(j+1)+tmp5*jd+tmp3] = p[ni-j+tmp4*(ni+1)+tmp5*(jd-5)+tmp3];
      }
      for (j = 0; j <= ni; j++) {
	p[j+tmp4+tmp5*jd+tmp3] = p[tmp4*(j+1)+tmp5*(jd-1)+tmp3];
      }
      break;
    case 10:
      for (j = 0; j <= ni; j++) {
	p[j+tmp4+tmp5*10+tmp3] = p[tmp4*(j+1)+tmp5*9+tmp3];
      }
      for (j = 0; j <= ni; j++) {
	p[tmp4*(j+1)+tmp5*10+tmp3] = p[j+tmp4+tmp5*6+tmp3];
      }
      for (j = 0; j <= ni; j++) {
	p[j+tmp4*(ni+1)+tmp5*10+tmp3] = p[ni+tmp4*(ni+1-j)+tmp5+tmp3];
      }
      for (j = 0; j <= ni; j++) {
	p[ni+tmp4*(j+1)+tmp5*10+tmp3] = p[ni-j+tmp4*(ni+1)+tmp5*5+tmp3];
      }
      break;
    }
  }

  return;
}

/*****************************************************************************/


void gme_grid(int lbounds, int gridsize, double *rlon, double *rlat,
	      double *blon, double *blat, int *imask,
              int ni, int nd, int ni2, int ni3)
{
  int im1s, im1e, im2s, im2e;
  int i, j;
  double *xn;
  double *rlonx, *rlatx;

  /* check gridsize */
  if ( (ni+1)*(ni+1)*nd != gridsize )
    {
      fprintf(stderr, "gme_grid: Calculation of the global GME grid failed (ni=%d)!\n", ni);
      if ( (ni+1)*(ni+1)*nd > gridsize )
	{
	  fprintf(stderr, "gme_grid: Resulting grid size is greater than the predetermined grid size of %d.\n", gridsize);
	  fprintf(stderr, "gme_grid: Maybe this is only a part of a global GME grid without further information.\n");
	}
      exit (-1);
    }

  xn    = (double*) malloc(gridsize*3*sizeof(double));
  rlonx = (double*) malloc((ni+3)*(ni+3)*nd*sizeof(double));
  rlatx = (double*) malloc((ni+3)*(ni+3)*nd*sizeof(double));

  im1s = 0;
  im1e = ni;
  im2s = 1;
  im2e = ni+1;

  glo_coor(xn, rlon, rlat, im1s, im1e, im2s, im2e, nd, ni2, ni3);

  xd(rlon,im1s,im1e,im2s,im2e,nd,rlonx,im1s-1,im1e+1,im2s-1,im2e+1,nd);
  xd(rlat,im1s,im1e,im2s,im2e,nd,rlatx,im1s-1,im1e+1,im2s-1,im2e+1,nd);

  initmask(imask,ni,nd);

  if ( lbounds )
    {
      struct polygon *poly;
      
      poly  = (struct polygon*) malloc((ni+1)*(ni+1)*nd*sizeof(struct polygon));

      neighbours(rlonx,rlatx,im1s-1,im1e+1,im2s-1,im2e+1,nd, poly,im1s,im1e,im2s,im2e,nd);

      boundary(poly,im1s,im1e,im2s,im2e,nd);

      for ( i = 0; i < gridsize; i++ )
	{
	  for ( j = 0; j < poly[i].type; j++ )
	    {
	      blon[i*6+j] = poly[i].boundary[j].lon;
	      blat[i*6+j] = poly[i].boundary[j].lat;
	    }
	  if ( poly[i].type == pentagon )
	    {
	      blon[i*6+5] = blon[i*6+4];
	      blat[i*6+5] = blat[i*6+4];
	    }
	}
      
      free(poly);
    }

  free(rlatx);
  free(rlonx);
  free(xn);
}

/*
int main(int argc, char *argv[])
{
  struct polygon *poly;

  double *xn, *rlat, *rlon;
  double *rlatx, *rlonx; 

  double *area, total_area;

  int *mask;

  int ni;
  int ni2, ni3;
  int nd;

  int im1s, im1e, im2s, im2e;

  int j1, j2, jd;

  nd = 10;

  ni = 2;

  factorni(ni, &ni2, &ni3);  

  if (( poly  = (struct polygon*) malloc((ni+1)*(ni+1)*nd*sizeof(struct polygon))) == NULL) {
    perror("malloc poly");
    exit (-1);
  } 
  
  if (( xn  = (double*) malloc((ni+1)*(ni+1)*3*nd*sizeof(double))) == NULL) {
    perror("malloc xn");
    exit (-1);
  } 
  
  if ((rlon = (double*) malloc((ni+1)*(ni+1)*nd*sizeof(double))) == NULL) {
    perror("malloc rlon");
    exit (-1);
  } 
  
  if ((rlat = (double*) malloc((ni+1)*(ni+1)*nd*sizeof(double))) == NULL) {
    perror("malloc rlat");
    exit (-1);
  } 

  if ((rlonx = (double*) malloc((ni+3)*(ni+3)*nd*sizeof(double))) == NULL) {
    perror("malloc rlonx");
    exit (-1);
  } 
  
  if ((rlatx = (double*) malloc((ni+3)*(ni+3)*nd*sizeof(double))) == NULL) {
    perror("malloc rlatx");
    exit (-1);
  } 

  if ((mask = (int*) malloc((ni+1)*(ni+1)*nd*sizeof(int))) == NULL) {
    perror("malloc mask");
    exit (-1);
  } 

  if (( area  = (double*) malloc((ni+1)*(ni+1)*nd*sizeof(double))) == NULL) {
    perror("malloc area");
    exit (-1);
  } 

  im1s = 0;
  im1e = ni;
  im2s = 1;
  im2e = ni+1;

  glo_coor(xn,rlon,rlat,im1s,im1e,im2s,im2e,nd,ni2,ni3);

  xd(rlon,im1s,im1e,im2s,im2e,nd,rlonx,im1s-1,im1e+1,im2s-1,im2e+1,nd);
  xd(rlat,im1s,im1e,im2s,im2e,nd,rlatx,im1s-1,im1e+1,im2s-1,im2e+1,nd);

  initmask(mask,ni,nd);

  {
    FILE *out;
    int id1, id2, ioffset;

    id1 = ni+1;
    id2 = id1*(ni+1);
    ioffset = -(id1+id2);
    
    if ((out  = fopen("mask.dat", "w")) == NULL) {
      perror("couldn't open mask.dat");
      exit (-1);
    }

    for (jd = 1; jd <= nd; jd++) {
      fprintf(out,"%d-------------------------------------------------\n", jd);
      for (j2 = 1; j2 <= ni+1; j2++) {
	for (j1 = 0; j1 <= ni; j1++) {
	      fprintf(out,"%8d", mask[j1+id1*j2+id2*jd+ioffset]);
	} 
        fprintf(out,"\n");
      }
    }

    fclose (out);
  }

  neighbours(rlonx,rlatx,im1s-1,im1e+1,im2s-1,im2e+1,nd,
             poly,im1s,im1e,im2s,im2e,nd);

  boundary(poly,im1s,im1e,im2s,im2e,nd);

  {
    FILE *out;
    int id1, id2, ioffset;
    int jm;

    id1 = ni+1;
    id2 = id1*(ni+1);
    ioffset = -(id1+id2);

    if ((out  = fopen("dual.dat", "w")) == NULL) {
      perror("couldn't open dual.dat");
      exit (-1);
    }

    for (jd = 1; jd <= nd; jd++) {
      for (j2 = 1; j2 <= ni+1; j2++) {
	for (j1 = 0; j1 <= ni; j1++) {
 	  if (mask[j1+id1*j2+id2*jd+ioffset]) {
	    fprintf(out,">\n");
	    for (jm = 1; jm <= poly[j1+id1*j2+id2*jd+ioffset].type; jm++) {
	      fprintf(out,"%8.2f%8.2f\n",
		 pid180*poly[j1+id1*j2+id2*jd+ioffset].boundary[jm-1].lon,
		 pid180*poly[j1+id1*j2+id2*jd+ioffset].boundary[jm-1].lat);
	    }
 	  } 
	}
      }
    }
    fclose (out);
  }

  {
    struct geo p1, p2, p3;

    struct cart c1, c2, c3;

    int id1, id2, ioffset;
    int jm, jl;
    
    id1 = ni+1;
    id2 = id1*(ni+1);
    ioffset = -(id1+id2);
   
    total_area = 0.0;

    for (jd = 1; jd <= nd; jd++) {
      for (j2 = 1; j2 <= ni+1; j2++) {
	for (j1 = 0; j1 <= ni; j1++) {
	  area[j1+id1*j2+id2*jd+ioffset] = 0.0;
 	  if (mask[j1+id1*j2+id2*jd+ioffset]) {
	    p3.lon = poly[j1+id1*j2+id2*jd+ioffset].center.lon; 
	    p3.lat = poly[j1+id1*j2+id2*jd+ioffset].center.lat;
	    c3 = gc2cc(&p3);
	    for (jm = 1; jm <= poly[j1+id1*j2+id2*jd+ioffset].type; jm++) {
	      jl = jm-1;
	      if (jm == poly[j1+id1*j2+id2*jd+ioffset].type) {
		p1.lon = poly[j1+id1*j2+id2*jd+ioffset].boundary[0].lon;
                p1.lat = poly[j1+id1*j2+id2*jd+ioffset].boundary[0].lat;
		p2.lon = poly[j1+id1*j2+id2*jd+ioffset].boundary[jl].lon;
                p2.lat = poly[j1+id1*j2+id2*jd+ioffset].boundary[jl].lat;
		c1 = gc2cc(&p1);
		c2 = gc2cc(&p2);
	      } else {
		p1.lon = poly[j1+id1*j2+id2*jd+ioffset].boundary[jl].lon;
                p1.lat = poly[j1+id1*j2+id2*jd+ioffset].boundary[jl].lat;
		p2.lon = poly[j1+id1*j2+id2*jd+ioffset].boundary[jl+1].lon;
                p2.lat = poly[j1+id1*j2+id2*jd+ioffset].boundary[jl+1].lat;
		c1 = gc2cc(&p1);
		c2 = gc2cc(&p2);
	      }
	      area[j1+id1*j2+id2*jd+ioffset] = area[j1+id1*j2+id2*jd+ioffset] + areas(&c1, &c2, &c3);
	    }
	    total_area = total_area+area[j1+id1*j2+id2*jd+ioffset];
	  }
	}
      } 
    }
  }

  free(poly);
  free(xn);
  free(rlon);
  free(rlat);
  free(rlonx);
  free(rlatx);
  free(mask);
  free(area);

  return(0);
}
*/
