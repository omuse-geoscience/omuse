#ifndef _COUNTER_H
#define _COUNTER_H

typedef struct {
  double cputime;
  char   mark[32];
}
counter_t;

void counter_start(counter_t *counter);
void counter_stop(counter_t *counter);
double counter_cputime(counter_t counter);

#endif  /* _COUNTER_H */
