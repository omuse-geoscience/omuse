#ifndef _TABLEPAR_H
#define _TABLEPAR_H

enum {
  TABLE_DUP_NAME = 1 << 0,
  TABLE_DUP_LONGNAME = 1 << 1,
  TABLE_DUP_UNITS = 1 << 2,
};

typedef struct
{
  int   id;	     /* Parameter number (GRIB) */
  int dupflags;      /* keep track of which attributes got strdup'ed */
  const char *name;	     /* Parameter name */
  const char *longname;    /* Parameter long name */
  const char *units;	     /* Parameter units */
}
PAR;


static void tableLink(int tableID, const PAR *pars, int npars);
int tableDef(int modelID, int tablegribID, const char *tablename);

#endif
/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
