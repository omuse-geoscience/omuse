#ifdef HAVE_CONFIG_H
# include "config.h"
#endif

#include <assert.h>
#include <inttypes.h>
#include <stdlib.h>
#include <stdio.h>

#include "ensure_array_size.h"
#include "var_cksum.h"
#include "stream_cksum.h"

static struct cksum_table *
read_table(const char *table_fname, size_t *table_len)
{
  struct cksum_table *table = NULL;
  FILE *tablefp;
  size_t table_size = 0, table_used = 0;
  uint32_t cksum_temp;
  int code;

  if (!(tablefp = fopen(table_fname, "r")))
    {
      perror("failed to open table file");
      *table_len = (size_t)-1;
      return NULL;
    }
  while (fscanf(tablefp, "%08"PRIx32" %d\n", &cksum_temp, &code) == 2)
    {
      ENSURE_ARRAY_SIZE(table, table_size, table_used + 1);
      table[table_used].code = code;
      table[table_used].cksum = cksum_temp;
      ++table_used;
    }
  fclose(tablefp);
  *table_len = table_used;
  return table;
}


int main(int argc, char *argv[])
{
  char *fname = "example.grb", *table_fname = "example.cksum";

  if (argc > 1)
    fname = argv[1];
  if (argc > 2)
    table_fname = argv[2];

  // compute checksums from data file
  size_t nvars;
  struct cksum_table *file_vars = cksum_stream(fname, &nvars);
  if (!file_vars)
    exit(EXIT_FAILURE);
  // check checksums from table file
  int retcode;
  {
    size_t num_ref_entries;
    struct cksum_table *ref_var_table
      = read_table(table_fname, &num_ref_entries);
    if (num_ref_entries == (size_t)-1)
      exit(EXIT_FAILURE);
    retcode
      = compare_checksums(file_vars, nvars, "file",
                          ref_var_table, num_ref_entries, "reference table");
    free(ref_var_table);
  }

  return retcode;
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
