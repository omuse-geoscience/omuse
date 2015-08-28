#ifndef  _TIMEBASE_H
#define  _TIMEBASE_H

#include <inttypes.h>

/* date format:  YYYYMMDD */
/* time format:  hhmmss   */

void decode_julday(int calendar, int julday, int *year, int *mon, int *day);
int  encode_julday(int calendar, int year, int month, int day);

int date_to_julday(int calendar, int date);
int julday_to_date(int calendar, int julday);

int time_to_sec(int time);
int sec_to_time(int secofday);

void   julday_add_seconds(int64_t seconds, int *julday, int *secofday);
void   julday_add(int days, int secs, int *julday, int *secofday);
double julday_sub(int julday1, int secofday1, int julday2, int secofday2, int *days, int *secs);

void encode_juldaysec(int calendar, int year, int month, int day, int hour, int minute, int second, int *julday, int *secofday);
void decode_juldaysec(int calendar, int julday, int secofday, int *year, int *month, int *day, int *hour, int *minute, int *second);

#endif  /* _TIMEBASE_H */
