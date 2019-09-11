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
int32_t step();
int32_t cleanup_code();
#endif
