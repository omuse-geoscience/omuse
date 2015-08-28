#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>

#include "cdi.h"
#include "dmemory.h"

#include "pio_impl.h"
#include "pio_util.h"

struct cons
{
  void * val;
  struct cons * next;
};

struct listSet {
  struct cons *head, *tail;
  valDestroyFunction valDestroy;
  eqPredicate keyCompare;
  int count;
};

listSet *listSetNew( valDestroyFunction vD, eqPredicate kC )
{
  listSet *myq;

  myq = (listSet*) xmalloc( sizeof ( listSet ));

  myq->head = NULL;
  myq->tail = NULL;
  myq->valDestroy = vD;
  myq->keyCompare = kC;
  myq->count = 0;

  return myq;
}

void
listSetDelete(listSet *q)
{
  struct cons *curr, *succ;

  if ( q->head )
    {
      curr = q->head;

      while ( curr )
        {
          succ = curr->next;
          ( *( q->valDestroy )) ( curr->val );
          free ( curr );
          curr = succ;
        }
    }

  free ( q );

  return;
}

int
listSetAdd(listSet *q, void *v)
{
  struct cons *newCons;

  {
    struct cons *p;
    for (p = q->head; p; p = p->next)
      // ensure unique keys
      if (q->keyCompare(v, p->val))
        return -1;
  }

  if ((newCons = (struct cons*) malloc(sizeof(struct cons))) == NULL)
    {
      perror ( "pio_listSet: listSetAdd (): Not enough memory" );
      /* FIXME: why not abort? */
      return 1;
    }

  newCons->val = v;
  newCons->next = NULL;


  if ( q->tail != NULL)
    q->tail->next = newCons;
  else
    q->head = newCons;

  q->tail = newCons;

  return q->count++;
}

int
listSetRemove(listSet *q, int (*predicate)(void *, void *),
              void *data)
{
  struct cons **p;

  for (p = &q->head; *p; p = &(*p)->next)
    if (predicate((*p)->val, data))
      {
        struct cons *rem = *p;
        if (rem == q->tail) q->tail = NULL;
        int iret = q->valDestroy(rem->val);
        *p = rem->next;
        free(rem);
        return iret;
      }
  return -1;
}

void *
listSetGet(listSet *q, int (*predicate)(void *, void *), void *data)
{
  struct cons *p;
  xassert(q && predicate);
  for (p = q->head; p; p = p->next)
    if (predicate(p->val, data))
      return p->val;
  return NULL;
}

bool
listSetIsEmpty(listSet *q)
{
  return q->head == NULL;
}


void
listSetForeach(listSet *q, void (*func)(void *elem, void *data), void *data)
{
  struct cons *p;
  for (p = q->head; p; p = p->next)
    func(p->val, data);
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
