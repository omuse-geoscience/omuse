#ifndef _CGRIBEX_H
#define _CGRIBEX_H

#include <stdio.h>
#include <sys/types.h>

#define  GRIB_MISSVAL  -9.E33

/* GRIB1 Level Types */
#define  GRIB1_LTYPE_SURFACE               1
#define  GRIB1_LTYPE_CLOUD_BASE            2
#define  GRIB1_LTYPE_CLOUD_TOP             3
#define  GRIB1_LTYPE_ISOTHERM0             4
#define  GRIB1_LTYPE_TOA                   8
#define  GRIB1_LTYPE_SEA_BOTTOM            9
#define  GRIB1_LTYPE_ATMOSPHERE           10
#define  GRIB1_LTYPE_99                   99
#define  GRIB1_LTYPE_ISOBARIC            100
#define  GRIB1_LTYPE_MEANSEA             102
#define  GRIB1_LTYPE_ALTITUDE            103
#define  GRIB1_LTYPE_HEIGHT              105
#define  GRIB1_LTYPE_SIGMA               107
#define  GRIB1_LTYPE_SIGMA_LAYER         108
#define  GRIB1_LTYPE_HYBRID              109
#define  GRIB1_LTYPE_HYBRID_LAYER        110
#define  GRIB1_LTYPE_LANDDEPTH           111
#define  GRIB1_LTYPE_LANDDEPTH_LAYER     112
#define  GRIB1_LTYPE_ISENTROPIC          113
#define  GRIB1_LTYPE_SEADEPTH            160  /* Depth Below Sea Level                                 */
#define  GRIB1_LTYPE_LAKE_BOTTOM         162  /* Lake or River Bottom                                  */
#define  GRIB1_LTYPE_SEDIMENT_BOTTOM     163  /* Bottom Of Sediment Layer                              */
#define  GRIB1_LTYPE_SEDIMENT_BOTTOM_TA  164  /* Bottom Of Thermally Active Sediment Layer             */
#define  GRIB1_LTYPE_SEDIMENT_BOTTOM_TW  165  /* Bottom Of Sediment Layer Penetrated By Thermal Wave   */
#define  GRIB1_LTYPE_MIX_LAYER           166  /* Mixing Layer                                          */
#define  GRIB1_LTYPE_99_MARGIN          1000

/* GRIB1 Data representation type (Grid Type) [Table 6] */
#define  GRIB1_GTYPE_LATLON                0  /*  latitude/longitude                                   */
#define  GRIB1_GTYPE_LATLON_ROT           10  /*  rotated latitude/longitude                           */
#define  GRIB1_GTYPE_LATLON_STR           20  /*  stretched latitude/longitude                         */
#define  GRIB1_GTYPE_LATLON_ROTSTR        30  /*  rotated and stretched latitude/longitude             */
#define  GRIB1_GTYPE_GAUSSIAN              4  /*  gaussian grid                                        */
#define  GRIB1_GTYPE_GAUSSIAN_ROT         14  /*  rotated gaussian grid                                */
#define  GRIB1_GTYPE_GAUSSIAN_STR         24  /*  stretched gaussian grid                              */
#define  GRIB1_GTYPE_GAUSSIAN_ROTSTR      34  /*  rotated and stretched gaussian grid                  */
#define  GRIB1_GTYPE_LCC                   3  /*  Lambert conformal                                    */
#define  GRIB1_GTYPE_SPECTRAL             50  /*  spherical harmonics                                  */
#define  GRIB1_GTYPE_GME                 192  /*  hexagonal GME grid                                   */

/*
 *  Macros for the indicator section ( Section 0 )
 */
#define  ISEC0_GRIB_Len             (isec0[ 0])  /*  Number of octets in the GRIB message              */
#define  ISEC0_GRIB_Version         (isec0[ 1])  /*  GRIB edition number                               */


/*
 *  Macros for the product definition section ( Section 1 )
 */
#define  ISEC1_TABLE4_MINUTE      0
#define  ISEC1_TABLE4_HOUR        1
#define  ISEC1_TABLE4_DAY         2
#define  ISEC1_TABLE4_3HOURS     10
#define  ISEC1_TABLE4_6HOURS     11
#define  ISEC1_TABLE4_12HOURS    12
#define  ISEC1_TABLE4_QUARTER    13
#define  ISEC1_TABLE4_30MINUTES  14


#define  ISEC1_CodeTable            (isec1[ 0])  /*  Version number of code table                 */
#define  ISEC1_CenterID             (isec1[ 1])  /*  Identification of centre                     */
#define  ISEC1_ModelID              (isec1[ 2])  /*  Identification of model                      */
#define  ISEC1_GridDefinition       (isec1[ 3])  /*  Grid definition                              */
#define  ISEC1_Sec2Or3Flag          (isec1[ 4])  /*  Section 2 or 3 included                      */
#define  ISEC1_Parameter            (isec1[ 5])  /*  Parameter indicator                          */
#define  ISEC1_LevelType            (isec1[ 6])  /*  Type of level indicator                      */
#define  ISEC1_Level1               (isec1[ 7])  /*  Level 1                                      */
#define  ISEC1_Level2               (isec1[ 8])  /*  Level 2                                      */
#define  ISEC1_Year                 (isec1[ 9])  /*  Year of century (YY)                         */
#define  ISEC1_Month                (isec1[10])  /*  Month (MM)                                   */
#define  ISEC1_Day                  (isec1[11])  /*  Day (DD)                                     */
#define  ISEC1_Hour                 (isec1[12])  /*  Hour (HH)                                    */
#define  ISEC1_Minute               (isec1[13])  /*  Minute (MM)                                  */
#define  ISEC1_TimeUnit             (isec1[14])  /*  Time unit indicator                          */
#define  ISEC1_TimePeriod1          (isec1[15])  /*  P1 Time period                               */
#define  ISEC1_TimePeriod2          (isec1[16])  /*  P2 Time period                               */
#define  ISEC1_TimeRange            (isec1[17])  /*  Time range indicator                         */
#define  ISEC1_AvgNum               (isec1[18])  /*  Number of products included in an average    */
#define  ISEC1_AvgMiss              (isec1[19])  /*  Number of products missing from an average   */
#define  ISEC1_Century              (isec1[20])  /*  Century                                      */
#define  ISEC1_SubCenterID          (isec1[21])  /*  Subcenter identifier                         */
#define  ISEC1_DecScaleFactor       (isec1[22])  /*  Decimal scale factor                         */
#define  ISEC1_LocalFLag            (isec1[23])  /*  Flag field to indicate local use in isec1    */

#define  ISEC1_ECMWF_LocalExtension (isec1[36])
#define  ISEC1_ECMWF_Class          (isec1[37])


/*
 *  Macros for the grid definition section ( Section 2 )
 */
#define  ISEC2_GridType             (isec2[ 0])  /* Data representation type */

/* Triangular grids */

#define  ISEC2_GME_NI2              (isec2[ 1])  /*  Number of factor 2 in factorisation of Ni    */
#define  ISEC2_GME_NI3              (isec2[ 2])  /*  Number of factor 3 in factorisation of Ni    */
#define  ISEC2_GME_ND               (isec2[ 3])  /*  Nubmer of diamonds                           */
#define  ISEC2_GME_NI               (isec2[ 4])  /*  Number of tri. subdiv. of the icosahedron    */
#define  ISEC2_GME_AFlag            (isec2[ 5])  /*  Flag for orientation of diamonds (Table A)   */
#define  ISEC2_GME_LatPP            (isec2[ 6])  /*  Latitude of pole point                       */
#define  ISEC2_GME_LonPP            (isec2[ 7])  /*  Longitude of pole point                      */
#define  ISEC2_GME_LonMPL           (isec2[ 8])  /*  Longitude of the first diamond               */
#define  ISEC2_GME_BFlag            (isec2[ 9])  /*  Flag for storage sequence (Table B)          */

/* Spherical harmonic coeficients */

#define  ISEC2_PentaJ               (isec2[ 1])  /*  J pentagonal resolution parameter            */
#define  ISEC2_PentaK               (isec2[ 2])  /*  K pentagonal resolution parameter            */
#define  ISEC2_PentaM               (isec2[ 3])  /*  M pentagonal resolution parameter            */
#define  ISEC2_RepType              (isec2[ 4])  /*  Representation type                          */
#define  ISEC2_RepMode              (isec2[ 5])  /*  Representation mode                          */

/* Gaussian grids */

#define  ISEC2_NumLon               (isec2[ 1])  /*  Number of points along a parallel (Ni)       */
#define  ISEC2_NumLat               (isec2[ 2])  /*  Number of points along a meridian (Nj)       */
#define  ISEC2_FirstLat             (isec2[ 3])  /*  Latitude of the first grid point             */
#define  ISEC2_FirstLon             (isec2[ 4])  /*  Longitude of the first grid point            */
#define  ISEC2_ResFlag              (isec2[ 5])  /*  Resolution flag: 128 regular grid            */
#define  ISEC2_LastLat              (isec2[ 6])  /*  Latitude of the last grid point              */
#define  ISEC2_LastLon              (isec2[ 7])  /*  Longitude of the last grid point             */
#define  ISEC2_LonIncr              (isec2[ 8])  /*  i direction increment                        */
#define  ISEC2_LatIncr              (isec2[ 9])  /*  j direction increment                        */
#define  ISEC2_NumPar               (isec2[ 9])  /*  Number of parallels between a pole and the E.*/
#define  ISEC2_ScanFlag             (isec2[10])  /*  Scanning mode flags                          */
#define  ISEC2_NumVCP               (isec2[11])  /*  Number of vertical coordinate parameters     */

/* Lambert */
#define  ISEC2_Lambert_Lov          (isec2[ 6])  /*  Orientation of the grid                      */
#define  ISEC2_Lambert_dx           (isec2[ 8])  /*  X-direction grid length                      */
#define  ISEC2_Lambert_dy           (isec2[ 9])  /*  Y-direction grid length                      */
#define  ISEC2_Lambert_ProjFlag     (isec2[12])  /*  Projection centre flag                       */
#define  ISEC2_Lambert_LatS1        (isec2[13])  /*  First lat at which the secant cone cuts the sphere */
#define  ISEC2_Lambert_LatS2        (isec2[14])  /*  Second lat at which the secant cone cuts the sphere */
#define  ISEC2_Lambert_LatSP        (isec2[19])  /*  Latitude of the southern pole                */
#define  ISEC2_Lambert_LonSP        (isec2[20])  /*  Longitude of the southern pole               */


#define  ISEC2_Reduced              (isec2[16])  /* 0: regular, 1: reduced grid                   */

#define  ISEC2_RowLonPtr            (&isec2[22])
#define  ISEC2_RowLon(i)            (isec2[22+i]) /* Number of points along each parallel         */

/* */

#define  ISEC2_LatSP                (isec2[12])  /* Latitude of the southern pole of rotation     */
#define  ISEC2_LonSP                (isec2[13])  /* Longitude of the southern pole of rotation    */

#define  FSEC2_RotAngle             (fsec2[ 0])  /* Angle of rotation                             */
#define  FSEC2_StrFact              (fsec2[ 1])  /* Stretching factor                             */

/*
 *  Macros for the bit map section ( Section 3 )
 */
#define  ISEC3_PredefBitmap         (isec3[ 0])  /* Predefined bitmap                             */
#define  ISEC3_MissVal              (isec3[ 1])  /* Missing data value for integers               */
#define  FSEC3_MissVal              (fsec3[ 1])  /* Missing data value for floats                 */

/*
 *  Macros for the binary data section ( Section 4 )
 */
#define  ISEC4_NumValues            (isec4[ 0])  /* Number of data values for encode/decode       */
#define  ISEC4_NumBits              (isec4[ 1])  /* Number of bits used for each encoded value    */
#define  ISEC4_NumNonMissValues     (isec4[20])  /* Number of non-missing values                  */




void  gribFixZSE(int flag);     /* 1: Fix ZeroShiftError of simple packed spherical harmonics */
void  gribSetConst(int flag);   /* 1: Don't pack constant fields on regular grids */
void  gribSetDebug(int debug);  /* 1: Debugging */
void  gribSetRound(int round);
void  gribSetRefDP(double refval);
void  gribSetRefSP(float  refval);
void  gribSetValueCheck(int vcheck);


void  gribExSP(int *isec0, int *isec1, int *isec2, float *fsec2, int *isec3,
               float *fsec3, int *isec4, float *fsec4, int klenp, int *kgrib,
               int kleng, int *kword, char *hoper, int *kret);

void  gribExDP(int *isec0, int *isec1, int *isec2, double *fsec2, int *isec3,
               double *fsec3, int *isec4, double *fsec4, int klenp, int *kgrib,
               int kleng, int *kword, char *hoper, int *kret);


const char *cgribexLibraryVersion(void);

void  gribDebug(int debug);
void  gribSetCalendar(int calendar);

void  gribDateTime(int *isec1, int *date, int *time);
int   gribRefDate(int *isec1);
int   gribRefTime(int *isec1);
int   gribTimeIsFC(int *isec1);

void  gribPrintSec0(int *isec0);
void  gribPrintSec1(int *isec0, int *isec1);
void  gribPrintSec2DP(int *isec0, int *isec2, double *fsec2);
void  gribPrintSec2SP(int *isec0, int *isec2, float  *fsec2);
void  gribPrintSec3DP(int *isec0, int *isec3, double *fsec3);
void  gribPrintSec3SP(int *isec0, int *isec3, float  *fsec3);
void  gribPrintSec4DP(int *isec0, int *isec4, double *fsec4);
void  gribPrintSec4SP(int *isec0, int *isec4, float  *fsec4);
void  gribPrintSec4Wave(int *isec4);

void  gribPrintALL(int nrec, long offset, long recpos, long recsize, unsigned char *gribbuffer);
void  gribPrintPDS(int nrec, long recpos, long recsize, unsigned char *gribbuffer);
void  gribPrintGDS(int nrec, long recpos, long recsize, unsigned char *gribbuffer);
void  gribPrintBMS(int nrec, long recpos, long recsize, unsigned char *gribbuffer);
void  gribPrintBDS(int nrec, long recpos, long recsize, unsigned char *gribbuffer);
void  gribCheck1(int nrec, long recpos, long recsize, unsigned char *gribbuffer);
void  gribRepair1(int nrec, long recsize, unsigned char *gribbuffer);

int   gribGetZip(long recsize, unsigned char *gribbuffer, long *urecsize);

int   gribBzip(unsigned char *dbuf, long dbufsize, unsigned char *sbuf, long sbufsize);
int   gribZip(unsigned char *dbuf, long dbufsize, unsigned char *sbuf, long sbufsize);
int   gribUnzip(unsigned char *dbuf, long dbufsize, unsigned char *sbuf, long sbufsize);

int   gribOpen(const char *filename, const char *mode);
void  gribClose(int fileID);

int   gribRead(int fileID, unsigned char *buffer, size_t *buffersize);
int   gribWrite(int fileID, unsigned char *buffer, size_t buffersize);
off_t gribGetPos(int fileID);
int   gribGetSize(int fileID);
int   gribCheckSeek(int fileID, long *offset, int *version);
int   gribFileSeek(int fileID, long *offset);
int   gribReadSize(int fileID);
int   gribVersion(unsigned char *buffer, size_t buffersize);

int   grib_info_for_grads(off_t recpos, long recsize, unsigned char *gribbuffer, int *intnum, float *fltnum, off_t *bignum);

double calculate_pfactor(const double* spectralField, long fieldTruncation, long subsetTruncation);

#endif  /* _CGRIBEX_H */ 

