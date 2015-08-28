#ifndef _COLOR_H
#define _COLOR_H

typedef struct{
  double z_low, z_high, i_dz;
  int rgb_low[3], rgb_high[3], rgb_diff[3];
  int annot;
  int skip;
}
LUT;


typedef struct {	/* For back-, fore-, and nan-colors */
  int rgb[3];
  int skip;
}
BFN_COLOR;


typedef struct {
  int       ncolors;
  LUT      *lut;
  BFN_COLOR bfn[3];
}
CPT;


int cptRead(FILE *fp, CPT *cpt);
int cptWrite(FILE *fp, CPT cpt);
int cptWriteC(FILE *fp, CPT cpt, const char *name);

#endif  /* _COLOR_H */
