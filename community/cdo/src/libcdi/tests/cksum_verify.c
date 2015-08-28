#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#if defined __cplusplus && !defined __STDC_FORMAT_MACROS
#define __STDC_FORMAT_MACROS
#endif

#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "cksum.h"

enum {
  block_size = 16,
};

int
main()
{
  unsigned char *test_data, *init_block;
  size_t num_blocks = 1000;
  if (!(init_block = calloc(block_size * block_size, 1U))
      || !(test_data = calloc((size_t)block_size * num_blocks, 1U)))
    return EXIT_FAILURE;
  /* this is supposed to be non-random */
  srand48(5L);
  init_block[7] = 15U;
  /* repeat block and rotate */
  for (size_t i = 1; i < block_size; ++i)
    {
      memcpy(init_block + block_size * i,
             init_block + block_size - i, i);
      memcpy(init_block + block_size * i + i, init_block,
             block_size - i);
    }
  for (size_t i = 0; i < num_blocks; ++i)
    {
      size_t block_idx = ((size_t)lrand48()) % block_size;
      memcpy(test_data + i, init_block + block_idx * block_size,
             block_size);
    }
  uint32_t cksum_result = memcrc(test_data, num_blocks * block_size);
  if (cksum_result != UINT32_C(0xc47779cd))
    {
      printf("unexpected crc result: 0x%8"PRIx32"\n", cksum_result);
      return EXIT_FAILURE;
    }
  return EXIT_SUCCESS;
}

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
