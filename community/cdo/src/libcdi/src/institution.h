#ifndef INSTITUTION_H
#define INSTITUTION_H

int
instituteUnpack(void *buf, int size, int *position, int originNamespace,
                void *context, int force_id);

void instituteDefaultEntries(void);

#endif

/*
 * Local Variables:
 * c-file-style: "Java"
 * c-basic-offset: 2
 * indent-tabs-mode: nil
 * show-trailing-whitespace: t
 * require-trailing-newline: t
 * End:
 */
