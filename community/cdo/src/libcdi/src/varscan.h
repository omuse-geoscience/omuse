#ifndef _VARSCAN_H
#define _VARSCAN_H

#ifndef _GRID_H
#  include "grid.h"
#endif


void varAddRecord(int recID, int param, int gridID, int zaxistype, int lbounds,
		  int level1, int level2, int level_sf, int level_unit, int prec,
		  int *pvarID, int *plevelID, int tsteptype, int numavg, int ltype1, int ltype2,
		  const char *name, const char *stdname, const char *longname, const char *units,
                  const var_tile_t *tiles, int *tile_index);

void varDefVCT(size_t vctsize, double *vctptr);
void varDefZAxisReference(int nlev, int nvgrid, unsigned char uuid[CDI_UUID_SIZE]);

int  varDefZaxis(int vlistID, int zaxistype, int nlevels, double *levels, int lbounds,
		 double *levels1, double *levels2, int vctsize, double *vct, char *name,
		 char *longname, char *units, int prec, int mode, int ltype);

void varDefMissval(int varID, double missval);
void varDefCompType(int varID, int comptype);
void varDefInst(int varID, int instID);
int  varInqInst(int varID);
void varDefModel(int varID, int modelID);
int  varInqModel(int varID);
void varDefTable(int varID, int tableID);
int  varInqTable(int varID);
void varDefEnsembleInfo(int varID, int ens_idx, int ens_count, int forecast_type);

void varDefTypeOfGeneratingProcess(int varID, int typeOfGeneratingProcess);
void varDefProductDefinitionTemplate(int varID, int productDefinitionTemplate);


void varDefOptGribInt(int varID, int tile_index, long lval, const char *keyword);
void varDefOptGribDbl(int varID, int tile_index, double dval, const char *keyword);
int varOptGribNentries(int varID);

int  zaxisCompare(int zaxisID, int zaxistype, int nlevels, int lbounds, const double *levels, char *longname, char *units, int ltype);

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
