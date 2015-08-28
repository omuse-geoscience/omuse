#include <stdlib.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <map>

#if defined(__cplusplus)
extern "C" {
#endif
#include "cdi.h"
#if defined (__cplusplus)
}
#endif
#include "cdi.hpp"

/*
 * CdiGrid
 * {
 */
CdiGrid::CdiGrid() { gridID = -1; }
CdiGrid::CdiGrid(int gridid) {
  char _name[CHARSIZE];
  char _xname[CHARSIZE], _xlongname[CHARSIZE], _xstdname[CHARSIZE], _xunits[CHARSIZE];
  char _yname[CHARSIZE], _ylongname[CHARSIZE], _ystdname[CHARSIZE], _yunits[CHARSIZE];

  gridID = gridid;
  type   = gridInqType(gridID);
  /*  strcpy(_name,gridNamePtr(type)); */
  name = std::string(gridNamePtr(type));

  size  = gridInqSize(gridID);
  xsize = gridInqXsize(gridID);
  ysize = gridInqYsize(gridID);
  xvalues.resize(xsize);
  yvalues.resize(ysize);
//  std::cout << "(gridID:" << gridID<< ") xsize: " << xsize << " ysize:" << ysize << std::endl;

  getValues();
  getBounds();

  prec = gridInqPrec(gridID);

  gridInqXname(gridID    , _xname);
  gridInqXlongname(gridID, _xlongname);
  gridInqXstdname(gridID , _xstdname);
  gridInqXunits(gridID   , _xunits);
  gridInqYname(gridID    , _yname);
  gridInqYlongname(gridID, _ylongname);
  gridInqYstdname(gridID , _ystdname);
  gridInqYunits(gridID   , _yunits);
  xname     = _xname;
  xlongname = _xlongname;
  xstdname  = _xstdname;
  xunits    = _xunits;
  yname     = _yname;
  ylongname = _ylongname;
  ystdname  = _ystdname;
  yunits    = _yunits;

}
CdiGrid::~CdiGrid() { gridID = -1; }

void CdiGrid::getValues()
{
  double *xvals, *yvals;

  hasXValues = gridInqXvals(gridID, NULL);
  hasYValues = gridInqYvals(gridID, NULL);
  if ( hasXValues )
  {
    xvals = (double *) malloc(xsize*sizeof(double));

    gridInqXvals(gridID, xvals);
    std::copy(xvals, xvals + xsize, xvalues.begin());
    free(xvals);
  }
  else
  {
    std::cout << "grid " << gridID << " has no YValues" << std::endl;
  }

  if ( hasYValues )
  {
    yvals = (double *) malloc(ysize*sizeof(double));
    gridInqYvals(gridID, yvals);
    std::copy(yvals, yvals + ysize, yvalues.begin());
    free(yvals);
  }
  else
  {
    std::cout << "grid " << gridID << " has no YValues" << std::endl;
  }
}

void CdiGrid::getBounds()
{
  double *xbnds, *ybnds;

  hasBounds = (gridInqXbounds(gridID,NULL) && gridInqYbounds(gridID,NULL));
  if (hasBounds)
  {
    ncorner  = ( gridInqType(gridID) == GRID_UNSTRUCTURED ) ? gridInqNvertex(gridID) : 4;
    ncorner *= size;
    xbnds    = (double *) malloc((ncorner*size)*sizeof(double));
    ybnds    = (double *) malloc((ncorner*size)*sizeof(double));
    gridInqYbounds(gridID, xbnds);
    gridInqXbounds(gridID, ybnds);
    xbounds.resize(ncorner*size);
    ybounds.resize(ncorner*size);
    std::copy(xbnds, xbnds + ncorner*size, xbounds.begin());
    std::copy(ybnds, ybnds + ncorner*size, ybounds.begin());
    free(xbnds);
    free(ybnds);
  }
  else
  {
    ncorner = 0;
    std::cout << "grid " << gridID << " has no XYBounds" << std::endl;
    xbounds.resize(ncorner*size);
    ybounds.resize(ncorner*size);
  }
}
/*
 * pointer value functions {
 */
void CdiGrid::getValuesAsPointer(double *xvals, double *yvals)
{
  if (hasBounds)
  {
    std::copy(xvalues.begin(), xvalues.end(), xvals);
    std::copy(yvalues.begin(), yvalues.end(), yvals);
  }
  else
  {
    std::cout << "No Values available." << std::endl;
  }
}

void
CdiGrid::getBoundsAsPointer(double *xbnds, double *ybnds)
{
  if (hasBounds)
  {
    std::copy(xbounds.begin(), xbounds.end(), xbnds);
    std::copy(ybounds.begin(), ybounds.end(), ybnds);
  }
  else
  {
    std::cout << "No bounds available." << std::endl;
  }
}
  /* } */
/*
 * float pointer functions {
 */
void
CdiGrid::getFloatVals(float *xvalsF, float *yvalsF)
{
  for (int i = 0; i < xsize; i++)
    xvalsF[i] = (float) xvalues[i];
  for (int i = 0; i < ysize; i++)
    yvalsF[i] = (float) yvalues[i];
}

void
CdiGrid::getFloatBounds(float *xbndsF, float *ybndsF)
{
  if (hasBounds)
  {
    for (int i = 0; i < ncorner*size; i++)
    {
      xbndsF[i] = (float) xbounds[i];
      ybndsF[i] = (float) ybounds[i];
    }
  }
  else
  {
    std::cout << "No bounds available." << std::endl;
  }
}
  /* } */
/* } */

/*
 * CdiTaxis
 * {
 */
CdiTaxis::CdiTaxis() { taxisID = -1; }
CdiTaxis::CdiTaxis(int vlistID) {
  taxisID   = vlistInqTaxis(vlistID);
  type      = taxisInqType(taxisID);
  ntsteps   = vlistNtsteps(vlistID);
  hasBounds = taxisHasBounds(taxisID);
  vdate     = taxisInqVdate(taxisID);
  vtime     = taxisInqVtime(taxisID);
  rdate     = taxisInqRdate(taxisID);
  rtime     = taxisInqRtime(taxisID);
  calendar  = taxisInqCalendar(taxisID);
  unit      = taxisInqTunit(taxisID);
  unitname  = tunitNamePtr(taxisID);
}
CdiTaxis::~CdiTaxis() { if (taxisID >= 0) taxisID = -1; }
/* } */

/*
 * CdiZaxis
 * {
 */
CdiZaxis::CdiZaxis() { zaxisID = -1; }
CdiZaxis::CdiZaxis(int zaxisid) {
  char name[CHARSIZE], longname[CHARSIZE], units[CHARSIZE];

  zaxisID = zaxisid;
  size    = zaxisInqSize(zaxisID);
  prec    = zaxisInqPrec(zaxisID);
  type    = zaxisInqType(zaxisID);
  ltype   = zaxisInqLtype(zaxisID);

  plevels  = (double *) malloc(size*sizeof(double));
  plbounds = (double *) malloc(size*sizeof(double));
  pubounds = (double *) malloc(size*sizeof(double));
  pweights = (double *) malloc(size*sizeof(double));
  zaxisInqLevels(zaxisID  ,plevels);
  zaxisInqLbounds(zaxisID ,plbounds);
  zaxisInqUbounds(zaxisID ,pubounds);
  zaxisInqWeights(zaxisID ,pweights);
  levels.resize(size); lbounds.resize(size); ubounds.resize(size); weights.resize(size);
  std::copy(plevels,plevels+size,levels.begin());
  std::copy(plbounds,plbounds+size,lbounds.begin());
  std::copy(pubounds,pubounds+size,ubounds.begin());
  std::copy(pweights,pweights+size,weights.begin());

  zaxisInqName(zaxisid    ,name);
  zaxisInqLongname(zaxisid,longname);
  zaxisInqUnits(zaxisid   ,units);
}
CdiZaxis::~CdiZaxis() { if (zaxisID >= 0) zaxisID = -1; }
/* } */

/*
 * CdiVariable
 * {
 */
CdiVariable::CdiVariable() { size = -1; }
CdiVariable::CdiVariable(int streamid,int vlistid, int varid) {
  char _name[CHARSIZE],_longname[CHARSIZE], _units[CHARSIZE], _stdname[CHARSIZE];
  streamID = streamid;
  vlistID  = vlistid;
  varID    = varid;
  code     = vlistInqVarCode(vlistID, varID);
  vlistInqVar(vlistID        , varID, &gridID    , &zaxisID, &timeID);
  vlistInqVarName(vlistID    , varID, _name);
  vlistInqVarLongname(vlistID, varID, _longname);
  vlistInqVarStdname(vlistID , varID, _stdname);
  vlistInqVarUnits(vlistID   , varID, _units);
  taxisID  = vlistInqTaxis(vlistID);
  size     = vlistInqVarSize(vlistID, varID);
  datatype = vlistInqVarDatatype(vlistID,varID);
  missval  = vlistInqVarMissval(vlistID,varID);
  name     = _name;
  longname = _longname;
  stdname  = _stdname;
  units    = _units;
}
CdiVariable::~CdiVariable(){ size = -1; }

void
CdiVariable::sinfo() { std::cout << "code:"  << std::endl; }

double *
CdiVariable::getValuesAsPointer()
{
  int levelID = 0, tsID = 0, nmiss;
  int vdate, vtime, nrecs;
  double *field;

  nrecs = streamInqTimestep(streamID, tsID);
  vdate = taxisInqVdate(taxisID);
  vtime = taxisInqVtime(taxisID);
  field = (double *) malloc(grid.size*sizeof(double));
  streamReadVarSlice(streamID, varID, levelID, field, &nmiss);

  return field;
}

double **
CdiVariable::getValuesWithLevelAsPointer(int tsID)
{
  int    levelID, nmiss, nrecs;
  double **fieldWithLevel, *field;

  nrecs          = streamInqTimestep(streamID, tsID);
  fieldWithLevel = (double **) malloc(zaxis.size*sizeof(double *));
  field          = (double *)  malloc(grid.size*sizeof(double));
  for (levelID = 0; levelID < zaxis.size; levelID++)
  {
    streamReadVarSlice(streamID, varID, levelID, field, &nmiss);
    fieldWithLevel[levelID] = field;
  }
  free(field);
  return fieldWithLevel;
}

void
CdiVariable::getValues()
{
  double *field = getValuesAsPointer();
  values.resize(grid.size);
  std::copy(field, field + grid.size, values.begin());
  free(field);
}

void
CdiVariable::getValuesWithLevel(int tsID) {
  double **fieldWithLevel = getValuesWithLevelAsPointer(tsID);

  valuesWithLevel.resize(zaxis.size);
  for (int levelID = 0; levelID < zaxis.size; levelID++)
  {
    valuesWithLevel[levelID].resize(grid.size);
    std::copy(fieldWithLevel[levelID], fieldWithLevel[levelID] + grid.size, valuesWithLevel[levelID].begin());
  }
  free(fieldWithLevel);
}

std::vector<float>
CdiVariable::getFValues()
{
  if (values.empty()) getValues();
  std::vector<float> retval(values.begin(),values.end());
  return retval;
}


std::vector< std::vector<float> >
CdiVariable::getFValuesWithLevel(int tsID)
{
  if (valuesWithLevel.empty()) getValuesWithLevel();
  std::vector< std::vector<float> > retval;
  for (std::vector< std::vector<double> >::const_iterator it = valuesWithLevel.begin(); it != valuesWithLevel.end(); it++)
  {
    std::vector<float> fvalues((*it).begin(),(*it).end());
    retval.push_back(fvalues);
  }
  return retval;
}

/* } */

/*
 * Cdi
 * {
 */
Cdi::Cdi(const char *path)  {
  streamID = streamOpenRead(path);
  vlistID  = streamInqVlist(streamID);
  nvars    = vlistNvars(vlistID);

  getTaxes();
  getZaxes();
  getGrids();
  getVars();
}

Cdi::~Cdi() {}

void
Cdi::getTaxes() {
  int taxisID;
  ntaxes         = 1;
  taxisID        = vlistInqTaxis(vlistID);
  taxes[taxisID] = CdiTaxis(vlistID);
}

void
Cdi::getZaxes() {
  int zaxisID;
  nzaxes = vlistNzaxis(vlistID);
  for (int i = 0; i < nzaxes; i++)
  {
    zaxisID        = vlistZaxis(vlistID, i);
//    std::cout << "getZaxes : zaxisID=" << zaxisID << std::endl;
    zaxes[zaxisID] = CdiZaxis(zaxisID);
  }
}

void
Cdi::getGrids() {
  int gridID;
  ngrids = vlistNgrids(vlistID);
  for (int i = 0; i < ngrids; i++)
  {
    gridID        = vlistGrid(vlistID, i);
//    std::cout << "getGrids : gridID=" << gridID << std::endl;
    grids[gridID] = CdiGrid(gridID);
  }
}

void
Cdi::getVars() {
  char name[CHARSIZE];
  int varID, code;
//  std::cout << vlistID << std::endl;
  for (varID = 0; varID < nvars; varID++)
  {
    code = vlistInqVarCode(vlistID, varID);
    codes.push_back(code);

    vlistInqVarName(vlistID, varID, name);
    varnames.push_back(name);

    CdiVariable _var = CdiVariable(streamID,vlistID,varID);
    _var.grid        = grids[_var.gridID];
    _var.zaxis       = zaxes[_var.zaxisID];
    _var.taxis       = taxes[_var.taxisID];
  
    variables.push_back(_var);

    var[name] = _var;
    varByCode[code]        = _var;
  }
}

void Cdi::griddes() {
  //Cdo::griddes("|" + streamID);
}
/* } */
