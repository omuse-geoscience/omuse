#include <stdio.h>
#include "cdo.h"


void set_text_color(FILE *fp, int attr, int fg)
{
  int bg = -1;

  if ( fp == stdout && !COLOR_STDOUT ) return;
  if ( fp == stderr && !COLOR_STDERR ) return;

  fprintf(fp, "%c[%d", 0x1B, attr);
  if ( fg >= 0 )
    {
      fprintf(fp, ";%d", fg+30);
      if ( bg >= 0 ) fprintf(fp, ";%d", bg+40);
    }
  fprintf(fp, "m");
}


void reset_text_color(FILE *fp)
{
  int attr = RESET;

  if ( fp == stdout && !COLOR_STDOUT ) return;
  if ( fp == stderr && !COLOR_STDERR ) return;

  fprintf(fp, "%c[%dm", 0x1B, attr);
}
