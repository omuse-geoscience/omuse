#ifndef STREAM_CKSUM_H
#define STREAM_CKSUM_H

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "var_cksum.h"

struct cksum_table *
cksum_stream(const char *fname, size_t *table_len);

#endif
