#ifndef _SPECSPACE_H
#define _SPECSPACE_H

#include "afterburner.h"

typedef struct {
  long nlon;
  long nlat;
  long ntr;
  long poldim;
  long ifax[10];
  double *trig;
  double *poli;
  double *pold;
  double *pol2;    /* only for uv2dv  */
  double *pol3;    /* only for uv2dv  */
  double *coslat;  /* only for scaluv with uv2dv */
  double *rcoslat; /* only for scaluv with dv2uv */
}
SPTRANS;

typedef struct {
  int ntr;
  int fdim;
  double *f1;
  double *f2;
}
DVTRANS;

void dv2ps(const double * restrict div, double * restrict pot, long nlev, long ntr);

SPTRANS *sptrans_new(int nlon, int nlat, int ntr, int flag);
void sptrans_delete(SPTRANS *sptrans);

DVTRANS *dvtrans_new(int ntr);
void dvtrans_delete(DVTRANS *dvtrans);

void trans_uv2dv(SPTRANS *sptrans, int nlev,
		 int gridID1, double *gu, double *gv,
		 int gridID2, double *sd, double *svo);

void trans_dv2uv(SPTRANS *sptrans, DVTRANS *dvtrans, int nlev,
		 int gridID1, double *sd, double *svo,
		 int gridID2, double *gu, double *gv);

void grid2spec(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut);
void spec2grid(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut);
void four2spec(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut);
void spec2four(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut);
void four2grid(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut);
void grid2four(SPTRANS *sptrans, int gridIDin, double *arrayIn, int gridIDout, double *arrayOut);

void spec2spec(int gridIDin, double *arrayIn, int gridIDout, double *arrayOut);
void speccut(int gridIDin, double *arrayIn, double *arrayOut, int *waves);

void spcut(double *arrayIn, double *arrayOut, int ntr, int *waves);


#endif
