#include <string>
#include <vector>
#include <map>
#define CHARSIZE      128

class CdiGrid {
  public:
    CdiGrid();
    CdiGrid(int gridid);
    ~CdiGrid();

    int gridID;
    int type, size, xsize, ysize, prec, ncorner;
    bool hasXValues, hasYValues, hasBounds;
    std::vector<double> xvalues, yvalues, xbounds, ybounds;


    std::string xname, xlongname, xstdname, xunits;
    std::string yname, ylongname, ystdname, yunits;      
    std::string name;

    void getValues();
    void getBounds();
    void getValuesAsPointer(double *xvals, double *yvals);
    void getBoundsAsPointer(double *xbnds, double *ybnds);
    void getFloatVals(float *xvals, float * yvals);
    void getFloatBounds(float *xbnds, float *ybnds);
      
  private:
    void determineGrid(int gridID);
};

class CdiTaxis {
  public:
    CdiTaxis();
    CdiTaxis(int vlistID);
    ~CdiTaxis();
    
    int taxisID;
    int ntsteps, unit;
    int rdate, rtime, vdate, vtime;
    int type, calendar, hasBounds;
    char name[CHARSIZE];
    const char *unitname;
};

class CdiZaxis {
  public:
    CdiZaxis();
    CdiZaxis(int zaxisid);
    ~CdiZaxis();

    int zaxisID;
    int type, ltype, size, prec;
    double *plevels, *plbounds, *pubounds, *pweights;
    std::vector<double> levels, lbounds, ubounds, weights;
    std::string name, longname, units;
};

class CdiVariable { 
  public:
    CdiVariable();
    CdiVariable(int streamid,int vlistid, int varid);
    ~CdiVariable();

    int    varID, zaxisID, gridID, taxisID, timeID, vlistID;
    int    size, code, datatype;
    int    streamID;
    std::string  name, longname, units, stdname;
    double missval;
            std::vector<double>   values;
    std::vector< std::vector<double> > valuesWithLevel;


    CdiGrid  grid;
    CdiZaxis zaxis;
    CdiTaxis taxis;

    void sinfo();
    void getValues(); 
    void getValuesWithLevel(int tsID = 0);
    std::vector<float>   getFValues();
    std::vector< std::vector<float> > getFValuesWithLevel(int tsID = 0);
    double  *getValuesAsPointer();
    double **getValuesWithLevelAsPointer(int tsID = 0);
};

class Cdi {
  public:
    Cdi(const char *path);
    ~Cdi();
    int streamID;
    int vlistID,nvars,nzaxes,ngrids,ntaxes,taxisID;

    std::vector <std::string> varnames;
    std::vector <int> codes;
    
    std::vector<CdiVariable>           variables;
    std::map<std::string, CdiVariable> var;
    std::map<int,         CdiVariable> varByCode;
    std::map<int        , CdiTaxis>    taxes;
    std::map<int        , CdiZaxis>    zaxes;
    std::map<int        , CdiGrid>     grids;

    void griddes();
  private:
    void getTaxes();
    void getZaxes();
    void getGrids();
    void getVars();
};
