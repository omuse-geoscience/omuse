#ifdef HAVE_CONFIG_H
#  include "config.h"
#endif

#include "cdi.h"
#include "dmemory.h"
#include "grid.h"
#include "institution.h"
#include "model.h"
#include "cdi_int.h"
#include "vlist.h"
#include "namespace.h"
#include "serialize.h"
#include "resource_handle.h"
#include "resource_unpack.h"
#include "taxis.h"
#include "zaxis.h"

/*****************************************************************************/

void reshUnpackResources(char * unpackBuffer, int unpackBufferSize,
                         void *context)
{
  int updateType, resH, originNamespace;
  int unpackBufferPos = 0;
  int numAssociations = 0, sizeAssociations = 16;
  struct streamAssoc *associations
    = (struct streamAssoc *)xmalloc(sizeof (associations[0])
                                    * (size_t)sizeAssociations);

  {
    int msgHdr[2];
    serializeUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos,
                    &msgHdr, 2, DATATYPE_INT, context);
    if (msgHdr[0] != START)
      xabort("error parsing resource serialization buffer");
    originNamespace = msgHdr[1];
  }
  while ( unpackBufferPos < unpackBufferSize )
    {
      serializeUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos,
                      &updateType, 1, DATATYPE_INT, context);
      if (updateType == END)
        break;
      switch (updateType)
	{
	case GRID:
	  gridUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos,
                     originNamespace, context, 1);
	  break;
	case ZAXIS:
	  zaxisUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos,
                      originNamespace, context, 1);
	  break;
	case TAXIS:
	  taxisUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos,
                      originNamespace, context, 1);
	  break;
	case INSTITUTE:
          instituteUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos,
                          originNamespace, context, 1);
	  break;
	case MODEL:
          modelUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos,
                      originNamespace, context, 1);
	  break;
	case STREAM:
          if (sizeAssociations == numAssociations)
            associations
              = xrealloc(associations,
                         sizeof (associations[0]) * (size_t)(sizeAssociations *= 2));
	  associations[numAssociations]
            = streamUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos,
                           originNamespace, context);
          ++numAssociations;
	  break;
	case VLIST:
          vlistUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos,
                      originNamespace, context, 1);
	  break;
        case RESH_DELETE:
          serializeUnpack(unpackBuffer, unpackBufferSize, &unpackBufferPos,
                          &resH, 1, DATATYPE_INT, context);
          reshDestroy(namespaceAdaptKey(resH, originNamespace));
          break;
	default:
	  xabort("Invalid/unexpected serialization type %d or transfer error!",
                 updateType);
	}
    }
  for (int i = 0; i < numAssociations; ++i)
    {
      cdiStreamSetupVlist(stream_to_pointer(associations[i].streamID),
                          namespaceAdaptKey(associations[i].vlistID,
                                            originNamespace),
                          namespaceAdaptKey(associations[i].vlistIDorig,
                                            originNamespace));
    }
  free(associations);
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
