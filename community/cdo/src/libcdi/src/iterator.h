/*
 * This file is for the use of iterator.c and the CdiIterator subclasses only.
 */

#ifndef INCLUDE_GUARD_CDI_ITERATOR_INT_H
#define INCLUDE_GUARD_CDI_ITERATOR_INT_H

#include "cdi.h"

#include <stdbool.h>

/*
class CdiIterator

An iterator is an object that identifies the position of one record in a file, where a record is defined as the data belonging to one level, timestep, and variable.
Using iterators to read a file can be significantly faster than using streams, because they can avoid building an index of the file.
For file formats like grib that do not provide an index within the file, this makes the difference between reading the file once or reading the file twice.

CdiIterator is an abstract base class. Which derived class is used depends on the type of the file. The class hierarchy currently looks like this:

    CdiIterator <|--+-- CdiFallbackIterator
                    |
                    +-- CdiGribIterator

The fallback implementation currently uses the stream interface of CDI under the hood to provide full functionality for all filetypes for which no iterator implementation exists yet.
*/
//TODO[NH]: Debug messages, print function.

struct CdiIterator {
  int filetype;      //This is used to dispatch calls to the correct subclass.
  bool isAdvanced;    //Used to catch inquiries before the first call to CdiIteratorNextField(). //XXX: Advanced is probably not a good word (initialized?)

  //The metadata that can be accessed by the inquiry calls.
  //While theoretically redundant, these fields allow the handling of most inquiry calls within the base class.
  //Only the name is excempted because it needs an allocation.
  //These fields are set by the subclasses in the xxxIterNextField() method.
  int datatype, timesteptype;
  int gridId;
  CdiParam param;

  //The status information for reading/advancing is added in the subclasses.
};

void baseIterConstruct(CdiIterator* me, int filetype);
const char* baseIter_constructFromString(CdiIterator* me, const char* description);     //Returns a pointer past the end of the parsed portion of the description string.
void baseIterDestruct(CdiIterator* me);

#endif
