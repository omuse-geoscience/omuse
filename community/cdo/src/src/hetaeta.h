#ifndef _HETAETA_H
#define _HETAETA_H

void hetaeta(int ltq, int ngp, const int *imiss,
	     int nlev1, const double *restrict ah1, const double *restrict bh1,
             const double *restrict fis1, const double *restrict ps1, 
             const double *restrict t1, const double *restrict q1,
             int nlev2, const double *restrict ah2, const double *restrict bh2, 
             const double *restrict fis2, double *restrict ps2, 
             double *restrict t2, double *restrict q2,
	     int nvars, double *restrict *restrict vars1, double *restrict *restrict vars2,
	     double *restrict tscor, double *restrict pscor,
	     double *restrict secor);

#endif  /* _HETAETA_H */
