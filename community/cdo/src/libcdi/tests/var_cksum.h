#ifndef VAR_CKSUM_H
#define VAR_CKSUM_H

#include <inttypes.h>
#include <stdlib.h>

struct cksum_table
{
  int code;
  uint32_t cksum;
};

/* returns EXIT_SUCCESS if a contains the same entries with the same
 * checksums as b, EXIT_FAILURE otherwise */
int
compare_checksums(struct cksum_table a[], size_t a_size, const char *src_a,
                  struct cksum_table b[], size_t b_size, const char *src_b);

#endif
