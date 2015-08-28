/*
  This file is part of CDO. CDO is a collection of Operators to
  manipulate and analyse Climate model Data.

  Copyright (C) 2003-2015 Uwe Schulzweida, <uwe.schulzweida AT mpimet.mpg.de>
  See COPYING file for copying and redistribution conditions.

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#ifndef _KVLIST_H
#define _KVLIST_H

void *kvlParseFile(const char *filename);

void kvlDelete(void *kvlist);

int kvlGetNumLists(void *kvlist);
const char *kvlGetListName(void *kvlist, int listID);

int kvlGetListNumElements(void *kvlist, int listID);
const char *kvlGetListElementName(void *kvlist, int listID, int elemID);
const char *kvlGetListElementValue(void *kvlist, int listID, int elemID);

#endif  /* _KVLIST_H */
