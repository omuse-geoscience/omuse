#ifndef INTERFACE_HPP
#define INTERFACE_HPP

#include <cstdint>

int32_t initialize();

int32_t commit_parameters();
int32_t recommit_parameters();
int32_t initialize_code();
int32_t test_grid(char* logFile);
int32_t step();
int32_t run_continuation();
int32_t cleanup_code();

int32_t get_u(int *i, int *j, int *k, double *var, int n);
int32_t get_v(int *i, int *j, int *k, double *var, int n);
int32_t get_w(int *i, int *j, int *k, double *var, int n);
int32_t get_p(int *i, int *j, int *k, double *var, int n);
int32_t get_t(int *i, int *j, int *k, double *var, int n);
int32_t get_s(int *i, int *j, int *k, double *var, int n);

int32_t get_num_parameter_sets(int *_num);
int32_t get_parameter_set_name(int i, char **name);
int32_t get_num_parameters(char *param_set_name, char *param_name, int *_num);
int32_t get_parameter_name(char *param_set_name, char *param_name, int i, char **name);
int32_t get_parameter_type(char *param_set_name, char *param_name, char **name);

int32_t get_bool_parameter(char *param_set_name, char *param_name, bool *result);
int32_t set_bool_parameter(char *param_set_name, char *param_name, bool val);
int32_t get_default_bool_parameter(char *param_set_name, char *param_name, bool *result);

int32_t get_char_parameter(char *param_set_name, char *param_name, char *result);
int32_t set_char_parameter(char *param_set_name, char *param_name, char val);
int32_t get_default_char_parameter(char *param_set_name, char *param_name, char *result);

int32_t get_double_parameter(char *param_set_name, char *param_name, double *result);
int32_t set_double_parameter(char *param_set_name, char *param_name, double val);
int32_t get_default_double_parameter(char *param_set_name, char *param_name, double *result);

int32_t get_int_parameter(char *param_set_name, char *param_name, int *result);
int32_t set_int_parameter(char *param_set_name, char *param_name, int val);
int32_t get_default_int_parameter(char *param_set_name, char *param_name, int *result);

int32_t get_string_parameter(char *param_set_name, char *param_name, char **result);
int32_t set_string_parameter(char *param_set_name, char *param_name, char *val);
int32_t get_default_string_parameter(char *param_set_name, char *param_name, char **result);
#endif
