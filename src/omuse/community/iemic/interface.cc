#include <Epetra_Comm.h>
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "Ocean.H"
#include "Continuation.H"
#include "Utils.H"

#include "interface.hpp"
#include "paramset.hpp"

namespace
{
using namespace Utils;
using namespace Teuchos;

enum class Parameter : unsigned char
{ u = 0, v, w, p, t, s, count };

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wglobal-constructors"
#pragma GCC diagnostic ignored "-Wexit-time-destructors"
std::map<std::string, ParamSet> parameter_sets = {
    { "ocean", ParamSet("OMUSE Ocean Parameters", ParamSet::tag<Ocean>()) },
    { "continuation", ParamSet("OMUSE Continuation Parameters", ParamSet::tag<Continuation<RCP<Ocean>>>()) }
};

RCP<Epetra_Comm> comm;
RCP<Ocean> ocean;
RCP<Continuation<RCP<Ocean>>> continuation;
std::ofstream devNull("/dev/null");
#pragma GCC diagnostic pop

std::string resultString = "";

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
    using ContinuationType = Continuation<RCP<Ocean>>;

    try {
        ocean = rcp(new Ocean(comm, parameter_sets.at("ocean").get()));
        ocean->setPar("Combined Forcing", 0.0);
        ocean->getState('V')->PutScalar(0.0);
        continuation = rcp(new ContinuationType(ocean, parameter_sets.at("continuation").get()));

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

int32_t get_num_parameter_sets(int *_num)
{
    *_num = parameter_sets.size();
    return 0;
}

int32_t get_parameter_set_name(int i, char **name)
{
    try {
        auto it = parameter_sets.begin();

        if (i < 0 || i >= parameter_sets.size()) return -1;

        std::advance(it, i);

        resultString = it->first;
        *name = const_cast<char*>(resultString.c_str());

        return 0;
    } catch (const std::exception& exc) {
        std::cout << exc.what() << std::endl;
    } catch (...) {
        std::cout << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t get_num_parameters(char *param_set_name, char *param_name, int *_num)
{
    try {
        *_num = parameter_sets.at(param_set_name).get_num_params(param_name);
        return 0;
    } catch (const std::exception& exc) {
        std::cout << exc.what() << std::endl;
    } catch (...) {
        std::cout << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
get_parameter_name(char *param_set_name, char *param_name, int i, char **name)
{
    try {
        auto& paramSet = parameter_sets.at(param_set_name);
        resultString = paramSet.get_param_name(param_name, i);
        *name = const_cast<char*>(resultString.c_str());
        return 0;
    } catch (const std::exception& exc) {
        std::cout << exc.what() << std::endl;
    } catch (...) {
        std::cout << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t get_parameter_type(char *param_set_name, char *param_name, char **name)
{
    try {
        resultString = parameter_sets.at(param_set_name).get_param_type(param_name);
        *name = const_cast<char*>(resultString.c_str());
        return 0;
    } catch (const std::exception& exc) {
        std::cout << exc.what() << std::endl;
    } catch (...) {
        std::cout << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
get_bool_parameter(char *param_set_name, char *param_name, bool *result)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.get_param_value(param_name, *result);
        return 0;
    } catch (const std::exception& exc) {
        std::cout << exc.what() << std::endl;
    } catch (...) {
        std::cout << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
set_bool_parameter(char *param_set_name, char *param_name, bool val)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.set_param_value(param_name, val);
        return 0;
    } catch (const std::exception& exc) {
        std::cout << exc.what() << std::endl;
    } catch (...) {
        std::cout << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
get_default_bool_parameter(char *param_set_name, char *param_name, bool *result)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.get_default_param_value(param_name, *result);
        return 0;
    } catch (const std::exception& exc) {
        std::cout << exc.what() << std::endl;
    } catch (...) {
        std::cout << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
get_char_parameter(char *param_set_name, char *param_name, char *result)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.get_param_value(param_name, *result);
        return 0;
    } catch (const std::exception& exc) {
        std::cout << exc.what() << std::endl;
    } catch (...) {
        std::cout << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
set_char_parameter(char *param_set_name, char *param_name, char val)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.set_param_value(param_name, val);
        return 0;
    } catch (const std::exception& exc) {
        std::cout << exc.what() << std::endl;
    } catch (...) {
        std::cout << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
get_default_char_parameter(char *param_set_name, char *param_name, char *result)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.get_default_param_value(param_name, *result);
        return 0;
    } catch (const std::exception& exc) {
        std::cout << exc.what() << std::endl;
    } catch (...) {
        std::cout << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
get_double_parameter(char *param_set_name, char *param_name, double *result)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.get_param_value(param_name, *result);
        return 0;
    } catch (const std::exception& exc) {
        std::cout << exc.what() << std::endl;
    } catch (...) {
        std::cout << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
set_double_parameter(char *param_set_name, char *param_name, double val)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.set_param_value(param_name, val);
        return 0;
    } catch (const std::exception& exc) {
        std::cout << exc.what() << std::endl;
    } catch (...) {
        std::cout << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
get_default_double_parameter(char *param_set_name, char *param_name, double *result)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.get_default_param_value(param_name, *result);
        return 0;
    } catch (const std::exception& exc) {
        std::cout << exc.what() << std::endl;
    } catch (...) {
        std::cout << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
get_int_parameter(char *param_set_name, char *param_name, int *result)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.get_param_value(param_name, *result);
        return 0;
    } catch (const std::exception& exc) {
        std::cout << exc.what() << std::endl;
    } catch (...) {
        std::cout << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
set_int_parameter(char *param_set_name, char *param_name, int val)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.set_param_value(param_name, val);
        return 0;
    } catch (const std::exception& exc) {
        std::cout << exc.what() << std::endl;
    } catch (...) {
        std::cout << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
get_default_int_parameter(char *param_set_name, char *param_name, int *result)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.get_default_param_value(param_name, *result);
        return 0;
    } catch (const std::exception& exc) {
        std::cout << exc.what() << std::endl;
    } catch (...) {
        std::cout << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
get_string_parameter(char *param_set_name, char *param_name, char **result)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.get_param_value(param_name, resultString);
        *result = const_cast<char*>(resultString.c_str());
        return 0;
    } catch (const std::exception& exc) {
        std::cout << exc.what() << std::endl;
    } catch (...) {
        std::cout << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
set_string_parameter(char *param_set_name, char *param_name, char *val)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.set_param_value(param_name, std::string(val));
        return 0;
    } catch (const std::exception& exc) {
        std::cout << exc.what() << std::endl;
    } catch (...) {
        std::cout << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
get_default_string_parameter(char *param_set_name, char *param_name, char **result)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.get_default_param_value(param_name, resultString);
        *result = const_cast<char*>(resultString.c_str());
        return 0;
    } catch (const std::exception& exc) {
        std::cout << exc.what() << std::endl;
    } catch (...) {
        std::cout << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}
