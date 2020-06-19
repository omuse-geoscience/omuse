#include <Epetra_Comm.h>
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#endif

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>

#include "Ocean.H"
#include "Continuation.H"
#include "Utils.H"

#include <cstdint>
#include <map>

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
    { "Ocean", ParamSet("OMUSE Ocean Parameters", ParamSet::tag<Ocean>()) },
    { "Continuation", ParamSet("OMUSE Continuation Parameters", ParamSet::tag<Continuation<RCP<Ocean>>>()) }
};

RCP<Epetra_Comm> comm;
RCP<Ocean> ocean;
RCP<Continuation<RCP<Ocean>>> continuation;
std::ofstream devNull("/dev/null");
std::ofstream logStream("/dev/null");
std::ofstream outStream("/dev/stdout");
#pragma GCC diagnostic pop

std::string resultString = "";

int state_count=0;
map<int, RCP<Epetra_Vector>> states;
map<int, RCP<Epetra_Vector>> matrices;

int32_t get_param(RCP<Epetra_Vector> state, Parameter param, int *i, int *j, int *k, double *var, int n)
{
    auto paramCount = static_cast<size_t>(Parameter::count);
    auto N = ocean->getNdim();
    auto M = ocean->getMdim();
    auto L = ocean->getLdim();
    auto paramOffset = static_cast<unsigned char>(param);

    auto gvec = AllGather(*state);
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

int32_t set_log_file(char *filepath)
{
    try {
        std::ofstream logFile(filepath);
        logStream.swap(logFile);
        return 0;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t set_output_file(char *filepath)
{
    try {
        std::ofstream newOutFile(filepath);
        outStream.swap(newOutFile);
        return 0;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t _new_state(int *id)
{
    if( states.size()>INT_MAX/2 ) return -1;
    *id=state_count;
    while(states.count(*id)) *id=state_count++;
    RCP<Epetra_Vector> v = rcp(new Epetra_Vector(ocean->getState()->Map(), true));
    states[*id]=v;
    return 0;
}

int32_t _remove_state(int id)
{
  if(states.erase(id))
  {
      return 0;
  } else
  {
      return -1;
  }
}

int32_t _copy_state(int src, int target)
{
  if(! states.count(src) || ! states.count(target)) return -1;
  *states[target]=*states[src];
  return 0;
}

int32_t _to_str(int src, char **out)
{
  if(! states.count(src)) return -1;
  std::ostringstream ss;
  states[src]->Print(ss);
  resultString = ss.str();
  *out = const_cast<char*>(resultString.c_str());
  return 0;
}

int32_t _set_model_state(int src)
{
  if(! states.count(src)) return -1;
  *ocean->getState('V') = *states[src];
  return 0;
}

int32_t _get_rhs(int src, int target)
{
  if(! states.count(src) || ! states.count(target)) return -1;
  _set_model_state(src);
  ocean->computeRHS();
  *states[target] = *ocean->getRHS('V');
  return 0;
}

int32_t _update_state(int src1, int src2, double scal)
{
  if(! states.count(src1) || ! states.count(src2)) return -1;
  states[src1]->Update(1.0, *states[src2], scal);
  return 0;
}

int32_t _mul_state(int src, double x)
{
  if(! states.count(src)) return -1;
  states[src]->Scale(x);
  return 0;
}

int32_t _get_state_norm(int state, double *val)
{
  if(! states.count(state)) return -1;
  *val=norm(states[state]);
  return 0;
}

int32_t _dot(int src1, int src2, double *val)
{
  if(! states.count(src1) || ! states.count(src2)) return -1;
  *val = dot(states[src1], states[src2]);
  return 0;
}

int32_t _length(int state, int *val)
{
  if(! states.count(state)) return -1;
  *val = states[state]->GlobalLength();
  return 0;
}

int32_t _solve(int rhs, int target)
{
  if(! states.count(rhs) || ! states.count(target)) return -1;
  ocean->solve(states[rhs]);
  *states[target]=*ocean->getSolution('V');
  return 0;
}

int32_t _jacobian(int src)
{
  if(! states.count(src)) return -1;
  _set_model_state(src);
  ocean->computeJacobian();
  return 0;
}



int32_t initialize()
{
    try {
        comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
        outFile = rcpFromRef(outStream);
        cdataFile = rcpFromRef(devNull);
        return 0;
    } catch (const std::exception& exc) {
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t commit_parameters()
{
    try {
        ocean = rcp(new Ocean(comm, parameter_sets.at("Ocean").commit()));
        ocean->getState('V')->PutScalar(0.0);

        return 0;
    } catch (const std::exception& exc) {
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t commit_continuation_parameters()
{
    using ContinuationType = Continuation<RCP<Ocean>>;

    try {
        continuation = rcp(new ContinuationType(ocean, parameter_sets.at("Continuation").commit()));
        return 0;
    } catch (const std::exception& exc) {
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}


int32_t recommit_parameters()
{
    auto& ocean_params = parameter_sets.at("Ocean");
    try {
        ocean->setParameters(ocean_params.updates());
        ocean_params.update_committed_parameters(ocean->getParameters());

        return 0;
    } catch (const std::exception& exc) {
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t recommit_continuation_parameters()
{
    auto& continuation_params = parameter_sets.at("Continuation");
    try {
        continuation->setParameters(continuation_params.updates());
        continuation_params.update_committed_parameters(continuation->getParameters());

        return 0;
    } catch (const std::exception& exc) {
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

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

int32_t step_continuation()
{
    try {
        return continuation->run();
    } catch (const std::exception& exc) {
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t cleanup_code()
{
    continuation = Teuchos::null;
    ocean = Teuchos::null;
    comm = Teuchos::null;
    return 0;
}

int32_t
_load_xml_parameters(char *param_set_name, char *path)
{
    try {
        auto& paramSet = parameter_sets.at(param_set_name);
        paramSet.load_from_file(path);
        return 0;
    } catch (const std::exception& exc) {
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t save_xml_parameters(char *param_set_name, char *path)
{
    try {
        auto& paramSet = parameter_sets.at(param_set_name);
        paramSet.save_to_file(path);
        return 0;
    } catch (const std::exception& exc) {
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t get_u(int *i, int *j, int *k, double *var, int n)
{ return get_param(ocean->state_, Parameter::u, i, j, k, var, n); }

int32_t get_v(int *i, int *j, int *k, double *var, int n)
{ return get_param(ocean->state_, Parameter::v, i, j, k, var, n); }

int32_t get_w(int *i, int *j, int *k, double *var, int n)
{ return get_param(ocean->state_, Parameter::w, i, j, k, var, n); }

int32_t get_p(int *i, int *j, int *k, double *var, int n)
{ return get_param(ocean->state_, Parameter::p, i, j, k, var, n); }

int32_t get_t(int *i, int *j, int *k, double *var, int n)
{ return get_param(ocean->state_,Parameter::t, i, j, k, var, n); }

int32_t get_s(int *i, int *j, int *k, double *var, int n)
{ return get_param(ocean->state_, Parameter::s, i, j, k, var, n); }

// variants for any state
int32_t get_u_(int *i, int *j, int *k, int* sindex, double *var, int n)
{ return get_param(states[*sindex], Parameter::u, i, j, k, var, n); }

int32_t get_v_( int *i, int *j, int *k, int* sindex, double *var, int n)
{ return get_param(states[*sindex], Parameter::v, i, j, k, var, n); }

int32_t get_w_( int *i, int *j, int *k, int* sindex, double *var, int n)
{ return get_param(states[*sindex], Parameter::w, i, j, k, var, n); }

int32_t get_p_( int *i, int *j, int *k, int* sindex, double *var, int n)
{ return get_param(states[*sindex], Parameter::p, i, j, k, var, n); }

int32_t get_t_( int *i, int *j, int *k, int* sindex, double *var, int n)
{ return get_param(states[*sindex],Parameter::t, i, j, k, var, n); }

int32_t get_s_( int *i, int *j, int *k, int* sindex, double *var, int n)
{ return get_param(states[*sindex], Parameter::s, i, j, k, var, n); }


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
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t get_num_parameters(char *param_set_name, char *param_name, int *_num)
{
    try {
        *_num = parameter_sets.at(param_set_name).get_num_params(param_name);
        return 0;
    } catch (const std::exception& exc) {
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
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
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t _get_parameter_type(char *param_set_name, char *param_name, char **name)
{
    try {
        resultString = parameter_sets.at(param_set_name).get_param_type(param_name);
        *name = const_cast<char*>(resultString.c_str());
        return 0;
    } catch (const std::exception& exc) {
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
_get_bool_parameter(char *param_set_name, char *param_name, bool *result)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.get_param_value(param_name, *result);
        return 0;
    } catch (const std::exception& exc) {
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
_set_bool_parameter(char *param_set_name, char *param_name, bool val)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.set_param_value(param_name, val);
        return 0;
    } catch (const std::exception& exc) {
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
_get_default_bool_parameter(char *param_set_name, char *param_name, bool *result)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.get_default_param_value(param_name, *result);
        return 0;
    } catch (const std::exception& exc) {
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
_get_char_parameter(char *param_set_name, char *param_name, char **result)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        char paramValue;
        paramset.get_param_value(param_name, paramValue);
        resultString.clear();
        resultString.push_back(paramValue);
        *result = const_cast<char*>(resultString.c_str());
        return 0;
    } catch (const std::exception& exc) {
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
_set_char_parameter(char *param_set_name, char *param_name, char *val)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.set_param_value(param_name, *val);
        return 0;
    } catch (const std::exception& exc) {
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
_get_default_char_parameter(char *param_set_name, char *param_name, char **result)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        char paramValue;
        paramset.get_default_param_value(param_name, paramValue);
        resultString.clear();
        resultString.push_back(paramValue);
        *result = const_cast<char*>(resultString.c_str());
        return 0;
    } catch (const std::exception& exc) {
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
_get_double_parameter(char *param_set_name, char *param_name, double *result)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.get_param_value(param_name, *result);
        return 0;
    } catch (const std::exception& exc) {
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
_set_double_parameter(char *param_set_name, char *param_name, double val)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.set_param_value(param_name, val);
        return 0;
    } catch (const std::exception& exc) {
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
_get_default_double_parameter(char *param_set_name, char *param_name, double *result)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.get_default_param_value(param_name, *result);
        return 0;
    } catch (const std::exception& exc) {
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
_get_int_parameter(char *param_set_name, char *param_name, int *result)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.get_param_value(param_name, *result);
        return 0;
    } catch (const std::exception& exc) {
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
_set_int_parameter(char *param_set_name, char *param_name, int val)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.set_param_value(param_name, val);
        return 0;
    } catch (const std::exception& exc) {
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
_get_default_int_parameter(char *param_set_name, char *param_name, int *result)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.get_default_param_value(param_name, *result);
        return 0;
    } catch (const std::exception& exc) {
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
_get_string_parameter(char *param_set_name, char *param_name, char **result)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.get_param_value(param_name, resultString);
        *result = const_cast<char*>(resultString.c_str());
        return 0;
    } catch (const std::exception& exc) {
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
_set_string_parameter(char *param_set_name, char *param_name, char *val)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.set_param_value(param_name, std::string(val));
        return 0;
    } catch (const std::exception& exc) {
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t
_get_default_string_parameter(char *param_set_name, char *param_name, char **result)
{
    try {
        auto& paramset = parameter_sets.at(param_set_name);
        paramset.get_default_param_value(param_name, resultString);
        *result = const_cast<char*>(resultString.c_str());
        return 0;
    } catch (const std::exception& exc) {
        logStream << exc.what() << std::endl;
    } catch (...) {
        logStream << "Encountered unexpected C++ exception!" << std::endl;
    }

    return -1;
}

int32_t get_state_norm(double *val)
{
  *val=norm(ocean->getState('V'));
  return 0;
}



// norm rhs (cont)
// parameter (cont)
// psimin, psimax (ocean writeData)
// # newton iterations (continuation writeData)
