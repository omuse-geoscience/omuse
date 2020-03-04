#ifndef INTERFACE_HPP
#define INTERFACE_HPP

#include <cstdint>

int32_t initialize();

int32_t set_default_params();
int32_t get_ocean_params(char **xmlFile);
int32_t set_ocean_params(char *xmlFile);
int32_t get_continuation_params(char **xmlFile);
int32_t set_continuation_params(char *xmlFile);

int32_t commit_parameters();
int32_t recommit_parameters();
int32_t initialize_code();
int32_t test_grid(const char* logFile);
int32_t step();
int32_t cleanup_code();

int32_t get_u(int *i, int *j, int *k, double *var, int n);
int32_t get_v(int *i, int *j, int *k, double *var, int n);
int32_t get_w(int *i, int *j, int *k, double *var, int n);
int32_t get_p(int *i, int *j, int *k, double *var, int n);
int32_t get_t(int *i, int *j, int *k, double *var, int n);
int32_t get_s(int *i, int *j, int *k, double *var, int n);
#endif
