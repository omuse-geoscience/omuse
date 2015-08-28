#ifndef _SUBTYPE_H
#define _SUBTYPE_H


enum {
  /* subtype attributes wrt. TILES */
  SUBTYPE_ATT_TILEINDEX                 = 0,
  SUBTYPE_ATT_TOTALNO_OF_TILEATTR_PAIRS = 1,
  SUBTYPE_ATT_TILE_CLASSIFICATION       = 2,
  SUBTYPE_ATT_NUMBER_OF_TILES           = 3,
  SUBTYPE_ATT_NUMBER_OF_ATTR            = 4,
  SUBTYPE_ATT_TILEATTRIBUTE             = 5,
/* No. of different constants in the enumeration
   "subtype_attributes" */
  nSubtypeAttributes
};


/* Literal constants corresponding to the different constants of the
   enumeration "subtype_attributes". */
extern const char * const cdiSubtypeAttributeName[];

/* Data type specifying an attribute of a subtype (for example an
   attribute of a set of TILES) or an attribute of a subtype entry
   (for example an attribute of a single TILE). This data type is part
   of a linked list. */
struct subtype_attr_t {
  int   key, val;                                /* key/value pair */
  struct subtype_attr_t* next;                   /* next element in linked list */
};


/* Data type specifying a single entry of a subtype, for example a
   single TILE in a set of TILES. */
struct subtype_entry_t {
  int                     self;                  /* list entry index (0,...,nentries-1) */
  struct subtype_entry_t *next;                  /* next node in linked list */

  /* linked list with attributes for this subtype entry, ordered by its key values*/
  struct subtype_attr_t  *atts;
};


/* Data type specifying a variable subtype, for example a list of
   TILES. This can be interpreted as an additional axis like the
   vertical axis. */
typedef struct  {
  int                     self;                  /* resource handler ID */
  int                     subtype;               /* subtype kind: TILES, ... */
  int                     nentries;              /* counter: total no. of entries in list */

  struct subtype_entry_t  globals;               /* global attributes */

  /* list of subtype entries, e.g. the list of tiles, ordered by entry->self. */
  struct subtype_entry_t *entries;
  /* currently active subtype, e.g. GRIB2 tile index (for example for
     stream/vlist accesses): */
  int                     active_subtype_index;
} subtype_t;




/* prototypes: allocation and destruction */
void  subtypeAllocate(subtype_t **subtype_ptr2, int subtype);
int   subtypePush(subtype_t *subtype_ptr);
void  subtypeDestroyPtr(void *ptr);
void  subtypeDuplicate(subtype_t *subtype_ptr, subtype_t **dst);
struct subtype_entry_t* subtypeEntryInsert(subtype_t* head);

/* prototypes: accessing global attributes */
void  subtypePrint(int subtypeID);
void  subtypePrintPtr(subtype_t* subtype_ptr);
void  subtypeDefGlobalDataP(subtype_t *subtype_ptr, int key, int val);
void  subtypeDefGlobalData(int subtypeID, int key, int val);
int   subtypeGetGlobalData(int subtypeID, int key);
int   subtypeGetGlobalDataP(subtype_t *subtype_ptr, int key);
int   subtypeComparePtr(int s1_ID, subtype_t *s2);

/* prototypes: accessing subtype entries */
void  subtypeDefEntryDataP(struct subtype_entry_t *subtype_entry_ptr, int key, int val);


/* prototypes: tile implementations */
void  tilesetInsertP(subtype_t *s1, subtype_t *s2);

/* Construct a new subtype for a tile set. If a corresponding subtype
 * already exists, then we return this subtype ID instead. */
int vlistDefTileSubtype(int vlistID, subtype_t *tiles);

/* Insert a trivial one-tile-subtype */
int vlistInsertTrivialTileSubtype(int vlistID);


#endif
