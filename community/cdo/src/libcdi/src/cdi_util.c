#include <stdio.h>

void cdiDecodeParam(int param, int *pnum, int *pcat, int *pdis)
{
  unsigned uparam = (unsigned)param;
  unsigned upnum;

  *pdis = 0xff   & uparam;
  *pcat = 0xff   & uparam >> 8;
  upnum = 0xffff & uparam >> 16;
  if ( upnum > 0x7fffU ) upnum = 0x8000U - upnum;
  *pnum = (int)upnum;
}


int cdiEncodeParam(int pnum, int pcat, int pdis)
{
  unsigned uparam, upnum;

  if ( pcat < 0 || pcat > 255 ) pcat = 255;
  if ( pdis < 0 || pdis > 255 ) pdis = 255;

  upnum = (unsigned)pnum;
  if ( pnum < 0 ) upnum = (unsigned)(0x8000 - pnum);

  uparam = upnum << 16 | (unsigned)(pcat << 8) | (unsigned)pdis;

  return ((int)uparam);
}


void cdiDecodeDate(int date, int *year, int *month, int *day)
{
  int idate;

  *year  =  date / 10000;
  idate  = date - *year*10000;
  if ( idate < 0 ) idate = -idate;
  *month = idate / 100;
  *day   = idate - *month*100;
}


int cdiEncodeDate(int year, int month, int day)
{
  int date;
  int iyear;

  iyear = year;
  if ( iyear < 0 ) iyear = -iyear;
  date = iyear*10000 + month*100 + day;
  if ( year < 0 ) date = -date;

  return (date);
}


void cdiDecodeTime(int time, int *hour, int *minute, int *second)
{
  int itime;

  *hour   = time / 10000;
  itime   = time - *hour*10000;
  *minute = itime / 100;
  *second = itime - *minute*100;
}


int cdiEncodeTime(int hour, int minute, int second)
{
  int time;

  time = hour*10000 + minute*100 + second;

  return (time);
}


void cdiParamToString(int param, char *paramstr, int maxlen)
{
  int dis, cat, num;
  int len;

  cdiDecodeParam(param, &num, &cat, &dis);

  if ( dis == 255 && (cat == 255 || cat == 0 ) )
    len = sprintf(paramstr, "%d", num);
  else  if ( dis == 255 )
    len = sprintf(paramstr, "%d.%d", num, cat);
  else
    len = sprintf(paramstr, "%d.%d.%d", num, cat, dis);

  if ( len > ( maxlen-1) )
    fprintf(stderr, "Internal problem (%s): size of input string is too small!\n", __func__);
}


char *cdiUnitNamePtr(int cdi_unit)
{
  char *cdiUnits[] = {
    /*  0 */  "undefined",
    /*  1 */  "Pa",
    /*  2 */  "hPa",
    /*  3 */  "mm",
    /*  4 */  "cm",
    /*  5 */  "dm",
    /*  6 */  "m",
  };
  char *name;
  int size = (int) (sizeof(cdiUnits)/sizeof(char *));

  if ( cdi_unit > 0 && cdi_unit < size )
    name = cdiUnits[cdi_unit];
  else
    name = NULL;

  return (name);
}

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
