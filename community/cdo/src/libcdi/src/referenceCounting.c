#include "referenceCounting.h"

#include "error.h"

void cdiRefObject_construct(CdiReferencedObject* me)
{
  me->destructor = cdiRefObject_destruct;
  me->refCount = 1;
}

void cdiRefObject_retain(CdiReferencedObject* me)
{
  size_t oldCount = me->refCount++;
  xassert(oldCount && "A reference counted object was used after it was destructed.");
}

void cdiRefObject_release(CdiReferencedObject* me)
{
  size_t oldCount = me->refCount--;
  xassert(oldCount && "A reference counted object was released too often.");
  if(oldCount == 1)
    {
      me->destructor(me);
      free(me);
    }
}

void cdiRefObject_destruct(CdiReferencedObject* me)
{
  (void)me;
  /* Empty for now, but that's no reason not to call it! */
}
