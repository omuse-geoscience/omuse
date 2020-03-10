#include <Epetra_Comm.h>
#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

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
ParameterList oceanParams;
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
        int idx = paramOffset
                + (paramCount * i[x])
                + (paramCount * N * j[x])
                + (paramCount * N * M * k[x]);
        var[x] = vec->operator[](idx);
    }
    return 0;
}
}

int32_t initialize()
{
    try {
        comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
        outFile = rcpFromRef(std::cout);
        cdataFile = rcpFromRef(devNull);
        return set_default_params();
    } catch (const std::exception& exc) {
        std::cout << exc.what() << std::endl;
    } catch (...) {
        std::cout << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t set_default_params()
{
    try {
        oceanParams = Ocean::getDefaultInitParameters();
        *continuationParams = default_continuation_params();
        return 0;
    } catch (const std::exception& exc) {
        std::cout << exc.what() << std::endl;
    } catch (...) {
        std::cout << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t commit_parameters()
{
    using ContinuationType = Continuation<RCP<Ocean>, RCP<ParameterList>>;

    try {
        ocean = rcp(new Ocean(comm, oceanParams));
        ocean->setPar("Combined Forcing", 0.0);
        ocean->getState('V')->PutScalar(0.0);
        continuation = rcp(new ContinuationType(ocean, continuationParams));

        return 0;
    } catch (const std::exception& exc) {
        std::cout << exc.what() << std::endl;
    } catch (...) {
        std::cout << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t recommit_parameters()
{ return 1; }

int32_t initialize_code()
{ return 0; }

int32_t test_grid(char *fileName)
{
    int32_t result = 0;
    auto N = ocean->getNdim();
    auto M = ocean->getMdim();
    auto L = ocean->getLdim();
    auto paramCount = static_cast<size_t>(Parameter::count);
    std::ofstream tmp(fileName ? fileName : "/dev/null");

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < M; j++) {
            for (int k = 0; k < L; k++) {
                for (int p = 0; p < paramCount; p++) {
                    int model;
                    int i_, j_, k_, p_;
                    int idx = p
                            + (paramCount * i)
                            + (paramCount * N * j)
                            + (paramCount * N * M * k);

                    ocean->gid2coord(idx, model, i_, j_, k_, p_);
                    if (i != i_ || j != j_ || k != k_ || p != p_) {
                        tmp << "Mismatch for gid: " << idx << std::endl
                            << "i: " << i << " " << i_ << std::endl
                            << "j: " << j << " " << j_ << std::endl
                            << "k: " << k << " " << k_ << std::endl
                            << "p: " << p << " " << p_ << std::endl;

                        result++;
                    }
                }
            }
        }
    }

    return result;
}

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

int32_t get_nrange(int *_min, int *_max)
{
  *_min = 0;
  *_max = ocean->getNdim() - 1;
  return 0;
}

int32_t get_mrange(int *_min, int *_max)
{
  *_min = 0;
  *_max = ocean->getMdim() - 1;
  return 0;
}

int32_t get_lrange(int *_min, int *_max)
{
  *_min = 0;
  *_max = ocean->getLdim() - 1;
  return 0;
}
