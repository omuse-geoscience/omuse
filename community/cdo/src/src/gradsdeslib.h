#ifndef  _GRADSDESLIB_H
#define  _GRADSDESLIB_H

#define  gaint     int
#define  gadouble  double
#define  gafloat   float
#define  galloc(x,y)  malloc(x)
#define  gree(x,y)    free(x)
#define  gaprnt(i,ch) printf("%s",ch)


/* Handling of missing data values. After the data I/O is done, 
   grid values are tested to see if they are within a small range 
   (+-value/EPSILON) of the missing value. If true, then the undef 
   mask is set to 0. If false, then the grid data values are good, 
   and the undef mask is set to 1. Everywhere else in the code, 
   undef tests are done on the mask values, not the data. */

#define EPSILON 1e5


/* Date/time structure */
struct dt {
  gaint yr;
  gaint mo;
  gaint dy;
  gaint hr;
  gaint mn;
};

/* Structure that describes a variable in a file.  These structures
   are built in arrays that are hung off of gafile structures.         */
struct gavar {
  char varnm[128];             /* Variable description.                */
  char abbrv[16];              /* Variable abbreviation.               */
  char longnm[257];            /* netcdf/hdf var name if different     */
  gadouble units[16];          /* Units indicator.                     
                                  Vals 0-7 are for variable codes:
                                  grib, non-float data, nc/hdf dims
                                  Vals  8-11 are for grib level codes  */
  gaint offset;                /* Offset in grid elements of the start
                                  of this variable within a time group
                                  within this file.                    */
  gaint recoff;                /* Record (XY grid) offset of the start
                                  of this variable within a time group */
  gaint ncvid;                 /* netcdf vid for this variable         */
  gaint sdvid;                 /* hdf vid for this variable            */
  gaint levels;                /* Number of levels for this variable.
                                  0 is special and indiates one grid is
                                  available for the surface only.      */
  gaint dfrm;                  /* format  type indicator
                                  1 - unsigned char
                                  4 - int                              */
  gaint var_t ;                /* variable t transform                 */
  gadouble scale;              /* scale factor for unpacking data      */
  gadouble add;                /* offset value for unpacking data      */
  gadouble undef;              /* undefined value                      */
  gaint vecpair;               /* Variable has a vector pair           */
  gaint isu;                   /* Variable is the u-component of a vector pair */
  gaint isdvar;                /* Variable is a valid data variable (for SDF files) */
  gaint nvardims;              /* Number of variable dimensions        */
  gaint vardimids[100];        /* Variable dimension IDs.                */
};

/* Sructure for string substitution in templating -- the %ch template.  
   This forms a linked list chained from pchsub1 in gafile */
struct gachsub {
  struct gachsub *forw;       /* Forward pointer */
  gaint t1;                   /* First time for this substitution */
  gaint t2;                   /* Last time.  -99 indicates open ended */
  char *ch;                   /* Substitution string */
};

/* Structure for ensemble metadata */
struct gaens {
  char name[16];             /* name of ensemble */
  gaint length;              /* length of time axis */
  struct dt tinit;           /* initial time */
  gaint gt;                  /* initial time in grid units */
  gaint grbcode[4];          /* grib2 codes */
};


#define  MAX_RECLEN   512
#define  MAX_NAMELEN  512
typedef struct {
  char name[MAX_NAMELEN];
  char dnam[MAX_NAMELEN];
  char title[MAX_NAMELEN];
  int bswap;
  long gsiz;                   /* Number of elements in a grid (x*y)    */
                               /* This is for actual grid on disk,
                                  not psuedo grid (when pp in force) */
  long tsiz;                   /* Number of elements in an entire time
                                  group (all variables at all levels
                                  for one time).                        */
  gaint trecs;                 /* Number of records (XY grids) per time
                                  group.                                */
  long fhdr;
  struct gachsub *pchsub1;     /* Pointer to first %ch substitution */
  gaint ncflg;                 /* 1==netcdf  2==hdfsds */
  long xyhdr;
  int seqflg;
  int yrflg;
  int zrflg;
  int pa2mb;
  int calendar;
  /* init !? */
  FILE *infile;                /* File pointer.                         */
  gaint type;                  /* Type of file:  1 = grid
                                                 2 = simple station
                                                 3 = mapped station
                                                 4 = defined grid       */
  gadouble undef;              /* Global undefined value for this file  */
  gadouble ulow,uhi;           /* Undefined limits for missing data test  */
  gaint dnum[5];               /* Dimension sizes for this file.        */
  gaint vnum;                  /* Number of variables.                  */
  struct gavar *pvar1;         /* Pointer to an array of structures.
                                  Each structure in the array has info
                                  about the specific variable.          */
  struct gaens *ens1;          /* pointer to array of ensemble structures */
  gaint wrap;                  /* The grid globally 'wraps' in X        */
  gadouble (*gr2ab[5]) (double *, double);
                               /* Addresses of routines to do conversion
                                  from grid coordinates to absolute
                                  coordinates for X, Y, Z.  All Date/time
                                  conversions handled by gr2t.          */
  gadouble (*ab2gr[5]) (double *, double);
                               /* Addresses of routines to do conversion
                                  from absolute coordinates to grid
                                  coordinates for X,Y,Z.  All date/time
                                  conversions handled by t2gr.          */
  gadouble *grvals[5];         /* Pointers to conversion information for
                                  grid-to-absolute conversion routines. */
  gadouble *abvals[5];         /* Pointers to conversion information for
                                  absolute-to-grid conversion routines. */
  gaint linear[5];             /* Indicates if a dimension has a linear
                                  grid/absolute coord transformation
                                  (Time coordinate always linear).      */
  gaint idxflg;                /* File records are indexed; 1==grib,station 2==grib2 */
  gaint flt64;                 /* 20120711 Uwe Schulzweida: added support for 64 bit floats */ 
  gaint tmplat;                /* File name templating:
                                   3==templating on E and T 
                                   2==templating only on E 
                                   1==templating only on T, or when 
                                      ddf has 'options template', but no % in dset 
                                   0==no templating  */
  gaint *fnums;                /* File number for each time */
  gaint fnumc;                 /* Current file number that is open */
  gaint fnume;                 /* Current ensemble file number that is open */
}
dsets_t;

void dsets_init(dsets_t *dsets);

int read_gradsdes(char *filename, dsets_t *pfi);

gadouble liconv (gadouble *, gadouble);
gadouble gr2lev (gadouble *, gadouble);
gadouble lev2gr (gadouble *, gadouble);
void gr2t (gadouble *, gadouble, struct dt *);
char *gafndt (char *, struct dt *, struct dt *, gadouble *, 
              struct gachsub *, struct gaens *, gaint, gaint, gaint *);

gaint cmpwrd (char *ch1, char *ch2);
char *intprs (char *ch, int *val);
void gabswp (void *r, gaint cnt);
void gabswp2 (void *r, gaint cnt);

#endif  /* _GRADSDESLIB_H */
