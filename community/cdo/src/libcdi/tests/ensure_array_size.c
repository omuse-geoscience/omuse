#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <stdlib.h>

#include "ensure_array_size.h"
#include "error.h"

void
realloc_array(void **array, size_t elem_size, size_t *curr_array_size,
              size_t requested_size)
{
  const size_t array_inc_size = (1024 + elem_size - 1)/ elem_size;
  *curr_array_size = array_inc_size
    * ((requested_size + array_inc_size) / array_inc_size);
  *array = realloc(*array, *curr_array_size * elem_size);
  if (!*array)
    xabort("reallocation failed");
}
