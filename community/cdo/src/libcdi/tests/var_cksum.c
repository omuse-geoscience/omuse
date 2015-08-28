#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include "var_cksum.h"

#include <stdio.h>
#include <stdlib.h>

int
compare_checksums(struct cksum_table a[], size_t a_size, const char *src_a,
                  struct cksum_table b[], size_t b_size, const char *src_b)
{
  int checked_a[a_size], checked_b[b_size];
  size_t i, j;
  int retcode = EXIT_SUCCESS;

  for (i = 0; i < a_size; ++i)
    checked_a[i] = 0;
  for (j = 0; j < b_size; ++j)
    checked_b[j] = 0;

  for (j = 0; j < b_size; ++j)
    for (i = 0; i < a_size; ++i)
      if (a[i].code == b[j].code)
      {
        if (a[i].cksum != b[j].cksum)
          {
            fprintf(stderr, "checksum error for varID %d, code %d!\n"
                    "%08lx != %08lx\n", (int)i, a[i].code,
                    (unsigned long)a[i].cksum, (unsigned long)b[j].cksum);
            retcode = EXIT_FAILURE;
          }
        checked_a[i] = 1;
        checked_b[j] = 1;
        break;
      }

  for (i = 0; i < a_size; ++i)
    if (!checked_a[i])
      {
        fprintf(stderr, "variable %d, code %d from %s not checked!\n",
                (int)i, a[i].code, src_a);
        retcode = EXIT_FAILURE;
      }
  for (j = 0; j < b_size; ++j)
    if (!checked_b[j])
      {
        fprintf(stderr, "variable %d, code %d from %s not checked!\n",
                (int)j, b[j].code, src_b);
        retcode = EXIT_FAILURE;
      }
  return retcode;
}

