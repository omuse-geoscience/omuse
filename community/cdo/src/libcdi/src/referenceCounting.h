#ifndef INCLUDE_GUARD_CDI_REFERENCE_COUNTING
#define INCLUDE_GUARD_CDI_REFERENCE_COUNTING

#include "error.h"

#include <sys/types.h>
#include <stdlib.h>

/*
This is a base class for all objects that need reference counting.
A CdiReferencedObject has a reference count of one when it is constructed, refObjectRetain() increments the reference count, refObject Release() decrements it.
When the reference count reaches zero, the destructor function is called before the memory of the object is deallocated with free().

>>> Warning <<<
This code is currently not thread-safe.

We are currently using the C99 standard, which does not have atomic types.
Also, there are still tons of systems out there that have a gcc without wrong C11 atomics support
(__STDC_NO_ATOMICS__ not defined even though stdatomics.h is not even present).
Consequently, it is impossible to write preprocessor code to even check for the presence of atomic types.
So, we have two options: provide multithreading support by means of locks, or wait a year or two before doing this right.
I, for one, prefer doing things right.
*/
typedef struct CdiReferencedObject CdiReferencedObject;
struct CdiReferencedObject {
  //protected:
    void (*destructor)(CdiReferencedObject* me);  //Subclass constructors should set this to their own destructor.

  //private:    //Subclasses may read it to determine whether there is only one reference, though.
    size_t refCount;
};

void cdiRefObject_construct(CdiReferencedObject* me);
void cdiRefObject_retain(CdiReferencedObject* me);
void cdiRefObject_release(CdiReferencedObject* me);
void cdiRefObject_destruct(CdiReferencedObject* me);

#endif
