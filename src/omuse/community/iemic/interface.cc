#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include "Ocean.H"
#include "Continuation.H"
#include "Utils.H"

#include "interface.hpp"

namespace
{
using namespace Utils;
using namespace Teuchos;

enum class Parameter : unsigned char
{ u = 0, v, w, p, t, s, count };

ParameterList&
insert_sublist(ParameterList& list, const std::string& name)
{
    list.set(name, ParameterList(name));
    return list.get<ParameterList>(name);
}

ParameterList
default_ocean_params()
{
    ParameterList params("Ocean");
    params.set("Load state", false);
    params.set("Save state", false);
    params.set("Input file", "ocean.h5");
    params.set("Output file", "ocean.h5");
    params.set("Store everything", false);

    auto& thcm_params = insert_sublist(params, "THCM");
    thcm_params.set("Problem Description", "North Atlantic");
    thcm_params.set("Global Bound xmin", 286.0);
    thcm_params.set("Global Bound xmax", 350.0);
    thcm_params.set("Global Bound ymin", 10.0);
    thcm_params.set("Global Bound ymax", 74.0);
    thcm_params.set("Periodic", false);
    thcm_params.set("Depth hdim", 4000.0);
    thcm_params.set("Grid Stretching qz", 1.0);
    thcm_params.set("Global Grid-Size n", 16);
    thcm_params.set("Global Grid-Size m", 16);
    thcm_params.set("Global Grid-Size l", 16);
    thcm_params.set("Topography", 1);
    thcm_params.set("Flat Bottom", false);
    thcm_params.set("Read Land Mask", false);
    thcm_params.set("Land Mask", "mask_natl16");
    thcm_params.set("Coupled Temperature", 0);
    thcm_params.set("Coupled Salinity", 0);

    auto& start_params = insert_sublist(thcm_params, "Starting Parameters");
    start_params.set("Combined Forcing", 0.0);
    start_params.set("Solar Forcing", 0.0);
    start_params.set("Salinity Forcing", 1.0);
    start_params.set("Wind Forcing", 1.0);
    start_params.set("Temperature Forcing", 10.0);
    start_params.set("SPL1", 2.0e3);
    start_params.set("SPL2", 0.01);

    thcm_params.set("Mixing", 1);
    thcm_params.set("Rho Mixing", false);
    thcm_params.set("Taper", 1);
    thcm_params.set("Restoring Temperature Profile", 1);
    thcm_params.set("Restoring Salinity Profile", 1);
    thcm_params.set("Freshwater Forcing", 1);
    thcm_params.set("Levitus T", 1);
    thcm_params.set("Levitus S", 1);
    thcm_params.set("Wind Forcing", 2);
    thcm_params.set("Scaling", "THCM");

    return params;
}

ParameterList
default_continuation_params()
{
    ParameterList params("Continuation parameters");
    params.set("continuation parameter", "Combined Forcing");
    params.set("initial step size", 7.0e-3);
    params.set("minimum step size", 1.0e-8);
    params.set("maximum step size", 1.0e-1);
    params.set("destination 0", 1.0);
    params.set("maximum number of steps", -1);
    params.set("Newton tolerance", 1.0e-4);
    params.set("destination tolerance", 1.0e-4);
    params.set("minimum desired Newton iterations", 4);
    params.set("maximum desired Newton iterations", 5);
    params.set("maximum Newton iterations", 7);
    params.set("enable backtracking", false);
    params.set("backtracking steps", 0);
    params.set("state tangent scaling", 1.0e-2);
    params.set("enable Newton Chord hybrid solve", false);

    return params;
}

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wglobal-constructors"
#pragma GCC diagnostic ignored "-Wexit-time-destructors"
RCP<ParameterList> oceanParams(new ParameterList());
RCP<ParameterList> continuationParams(new ParameterList());

RCP<Epetra_Comm> comm;
RCP<Ocean> ocean;
RCP<Continuation<RCP<Ocean>, RCP<ParameterList>>> continuation;
std::ofstream devNull("/dev/null");
#pragma GCC diagnostic pop

int32_t get_param(Parameter param, int *i, int *j, int *k, double *var, int n)
{
    auto paramCount = static_cast<size_t>(Parameter::count);
    auto N = ocean->getNdim();
    auto M = ocean->getMdim();
    auto L = ocean->getLdim();
    auto paramOffset = static_cast<unsigned char>(param);

    auto gvec = AllGather(*ocean->state_);
    auto vec = gvec->operator()(0);

    for (int x = 0; x < n; x++) {
        int idx = (paramCount * i[x])
                + (paramCount * N * j[x])
                + (paramCount * N * M * k[x]);
        var[x] = vec->operator[](idx);
    }
    return 0;
}

}

int32_t initialize()
{
    comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
    outFile = rcpFromRef(std::cout);
    cdataFile = rcpFromRef(devNull);
    return set_default_params();
}

int32_t set_default_params()
{
    oceanParams->operator=(default_ocean_params());
    continuationParams->operator=(default_continuation_params());
    return 0;
}

int32_t set_ocean_params(char *xmlFile)
{
    ParameterList params;
    updateParametersFromXmlFile(xmlFile, ptrFromRef(params));
    oceanParams->operator=(params);
    return 0;
}

int32_t get_ocean_params(char **)
{ return 1; }

int32_t set_continuation_params(char *xmlFile)
{
    ParameterList params;
    updateParametersFromXmlFile(xmlFile, ptrFromRef(params));
    continuationParams->operator=(params);
    return 0;
}

int32_t get_continuation_params(char **)
{ return 1; }

int32_t commit_parameters()
{
    using ContinuationType = Continuation<RCP<Ocean>, RCP<ParameterList>>;
    ocean = rcp(new Ocean(comm, oceanParams));
    ocean->setPar("Combined Forcing", 0.0);
    ocean->getState('V')->PutScalar(0.0);
    continuation = rcp(new ContinuationType(ocean, continuationParams));
    return 0;
}

int32_t recommit_parameters()
{ return 1; }

int32_t initialize_code()
{ return 0; }

int32_t step()
{
    int status = continuation->run();

    return status;
}

int32_t cleanup_code()
{
    continuation = Teuchos::null;
    ocean = Teuchos::null;
    comm = Teuchos::null;
    return 0;
}

int32_t get_u(int *i, int *j, int *k, double *var, int n)
{ return get_param(Parameter::u, i, j, k, var, n); }

int32_t get_v(int *i, int *j, int *k, double *var, int n)
{ return get_param(Parameter::v, i, j, k, var, n); }

int32_t get_w(int *i, int *j, int *k, double *var, int n)
{ return get_param(Parameter::w, i, j, k, var, n); }

int32_t get_p(int *i, int *j, int *k, double *var, int n)
{ return get_param(Parameter::p, i, j, k, var, n); }

int32_t get_t(int *i, int *j, int *k, double *var, int n)
{ return get_param(Parameter::t, i, j, k, var, n); }

int32_t get_s(int *i, int *j, int *k, double *var, int n)
{ return get_param(Parameter::s, i, j, k, var, n); }
