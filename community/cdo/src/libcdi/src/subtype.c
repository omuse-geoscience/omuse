/* Subroutines and data structures for storing "subtypes".             */
/*                                                                     */
/* A subtype is, for example, a list of TILES. This can be interpreted */
/* as an additional axis like the vertical axis.                       */
/*                                                                     */
/* @author 02/2015 F. Prill, DWD                                       */
/*                                                                     */
/*  DATA LAYOUT:                                                       */
/*                                                                     */
/*  A subtype contains several "subtype entries", each of which        */
/*  contains a linked list of subtype attributes.                      */
/*                                                                     */
/*  The number of subtype entries is not specified in advance, but the */
/*  list of entries is itself dynamically growing. There is no         */
/*  guaranteed ordering of the entries, therefore each entry must be   */
/*  identifiable by its attributes.                                    */
/*                                                                     */
/*  [subtype_t]                                                        */
/*      |                                                              */
/*      |------- globals                  [subtype_entry_t]            */
/*      |          |--- atts              [subtype_attr_t]             */
/*      |                                                              */
/*      |------- entries                                               */
/*                 |- entry #0                                         */
/*                 |  |--- atts              [subtype_attr_t]          */
/*                 |- entry #1                                         */
/*                 |  |--- atts              [subtype_attr_t]          */
/*                 |- entry #2                                         */
/*                 .  |--- atts              [subtype_attr_t]          */
/*                 .                                                   */

#if defined (HAVE_CONFIG_H)
#  include "config.h"
#endif

#include "cdi.h"
#include "cdi_int.h"
#include "subtype.h"
#include "vlist.h"

/* Literal constants corresponding to the different subtypes of the
   enumeration "subtype_kind". */
const char* subtypeName[] = {
  "tileset"
};

const char * const cdiSubtypeAttributeName[] = {
  "tileIndex",
  "totalNumberOfTileAttributePairs",
  "tileClassification",
  "numberOfTiles",
  "numberOfTileAttributes",
  "tileAttribute"
};


/* prototypes: */
static int    subtypeCompareP    (subtype_t *z1, subtype_t *z2);
static void   subtypeDestroyP    ( void * subtype_ptr );
static void   subtypePrintP      ( void * subtype_ptr, FILE * fp );
static int    subtypeGetPackSize ( void * subtype_ptr, void *context);
static void   subtypePack        ( void * subtype_ptr, void * buffer, int size, int *pos, void *context);
static int    subtypeTxCode      ( void );

static const resOps subtypeOps = {
  (int (*) (void *, void *)) subtypeCompareP,
  (void (*)(void *))         subtypeDestroyP,
  (void (*)(void *, FILE *)) subtypePrintP,
  (int (*) (void *, void *)) subtypeGetPackSize,
                             subtypePack,
                             subtypeTxCode
};

enum {
  differ = 1,
};



/* ------------------------------------------------------------------- */
/* SUBROUTINES FOR ATTRIBUTE LISTS				       */
/* ------------------------------------------------------------------- */


int attribute_to_index(const char *key)
{
  if (key == NULL)  Error("Internal error!");
  for (int i=0; i<nSubtypeAttributes; i++)
    if ( strcmp(key, cdiSubtypeAttributeName[i]) == 0 ) return i;
  return -1;
}



/*
  @Function  subtypeAttrNewList
  @Title     Create new linked list of subtype attributes.
  @EndFunction
*/
struct subtype_attr_t* subtypeAttrNewList(struct subtype_entry_t* head, int key, int val)
{
  if (head == NULL)  Error("Internal error!");
  struct subtype_attr_t *ptr = (struct subtype_attr_t*) malloc(sizeof(struct subtype_attr_t));
  if(NULL == ptr)  Error("Node creation failed");
  ptr->key   = key;
  ptr->val   = val;
  ptr->next  = NULL;

  head->atts = ptr;
  return ptr;
}


/*
  @Function  subtypeAttrInsert

  @Title Add subtype attribute to linked list, s.t. the result is a
         smallest-to-largest ordered list.
  @EndFunction
*/
struct subtype_attr_t* subtypeAttrInsert(struct subtype_entry_t* head, int key, int val)
{
  if (head == NULL)  Error("Internal error!");
  if (head->atts == NULL)  return (subtypeAttrNewList(head, key, val));

  /* create new attribute */
  struct subtype_attr_t* ptr = (struct subtype_attr_t*) malloc(sizeof(struct subtype_attr_t));
  if(NULL == ptr)    Error("Node creation failed");

  ptr->key   = key;
  ptr->val   = val;
  ptr->next  = NULL;

  /* find the right place for insertion: */
  if (head->atts->key >= key) {
    /* insert at position 0 */
    ptr->next = head->atts;
    head->atts = ptr;
  } else {
    struct subtype_attr_t** predec = &head->atts;
    while (((*predec)->next != NULL) && ((*predec)->next->key < key)) {
      predec = &((*predec)->next);
    }
    ptr->next = (*predec)->next;
    (*predec)->next = ptr;  
  }
  return ptr;
}


/* Recursively free a linked list with attributes. */
void subtypeAttrDestroy(struct subtype_attr_t* head)
{
  if (head == NULL) return;
  subtypeAttrDestroy(head->next);
  free(head);
  head = NULL; 
}


/* Find an attribute in linked list by its key or return NULL
   otherwise. */
struct subtype_attr_t* subtypeAttrFind(struct subtype_attr_t* head, int key)
{
  if (head == NULL) 
    return NULL;
  else if (head->key == key)
    return head;
  else
    return subtypeAttrFind(head->next, key);
}


/* Recursively compares two subtype attribute lists under the implicit
   assumptions that both lists are ordered by their keys and that keys
   are unique. */
static int subtypeAttsCompare(struct subtype_attr_t *a1, struct subtype_attr_t *a2)
{
  if ((a1 == NULL) && (a2 == NULL)) 
    return 0;
  else if ((a1 == NULL) && (a2 != NULL)) 
    {
      return differ;
    }
  else if ((a1 != NULL) && (a2 == NULL)) 
    {
      return differ;
    }

  if (a1->key != a2->key) 
    {
      return differ;
    }
  if (a1->val != a2->val)
    return differ;
    
  return subtypeAttsCompare(a1->next, a2->next);
}


/* (Recursively) duplicate linked list of attributes. */
void subtypeAttsDuplicate(struct subtype_attr_t *a1, struct subtype_entry_t* dst)
{
  if (a1 == NULL)  return;
  /* duplicate "a1->key", "a1->val" */
  subtypeAttsDuplicate(a1->next, dst);
  (void) subtypeAttrInsert(dst, a1->key, a1->val);
}



/* ------------------------------------------------------------------- */
/* SUBROUTINES FOR LIST OF ENTRIES				       */
/* ------------------------------------------------------------------- */


/*
  @Function  subtypeEntryNewList
  @Title     Create new linked list of subtype entries.
  @EndFunction
*/
struct subtype_entry_t* subtypeEntryNewList(subtype_t* head)
{
  struct subtype_entry_t *ptr = (struct subtype_entry_t*) malloc(sizeof(struct subtype_entry_t));
  if(NULL == ptr)  Error("Node creation failed");
  ptr->atts      = NULL;
  ptr->next      = NULL;
  head->entries  = ptr;
  head->nentries = 0;
  ptr->self      = head->nentries++;
  return ptr;
}


/*
  @Function  subtypeEntryInsert

  @Title Add subtype entry to the head of a linked list.
  @EndFunction
*/
struct subtype_entry_t* subtypeEntryInsert(subtype_t* head)
{
  if (head == NULL)  Error("Internal error!");
  if (head->entries == NULL)  return (subtypeEntryNewList(head));

  /* create new entry */
  struct subtype_entry_t* ptr = (struct subtype_entry_t*) malloc(sizeof(struct subtype_entry_t));
  if(NULL == ptr)    Error("Node creation failed");

  ptr->atts     = NULL;
  ptr->self     = head->nentries++;

  /* find the right place for insertion: */
  if (head->entries->self >= ptr->self) {
    /* insert at position 0 */
    ptr->next     = head->entries;
    head->entries = ptr;
  } else {
    struct subtype_entry_t** predec = &head->entries;
    while (((*predec)->next != NULL) && ((*predec)->next->self < ptr->self)) {
      predec = &((*predec)->next);
    }
    ptr->next = (*predec)->next;
    (*predec)->next = ptr;  
  }
  return ptr;
}


/*
  @Function  subtypeEntryAppend

  @Title Append subtype entry to the end of a linked list.
  @EndFunction
*/
struct subtype_entry_t* subtypeEntryAppend(subtype_t* head)
{
  if (head == NULL)  Error("Internal error!");
  if (head->entries == NULL)  return (subtypeEntryNewList(head));

  /* create new entry */
  struct subtype_entry_t* ptr = (struct subtype_entry_t*) malloc(sizeof(struct subtype_entry_t));
  if(NULL == ptr)    Error("Node creation failed");

  ptr->atts     = NULL;
  ptr->next     = NULL;
  ptr->self     = head->nentries++;

  /* find last position of linked list */
  struct subtype_entry_t* prec_ptr = head->entries;
  while (prec_ptr->next != NULL)
    prec_ptr = prec_ptr->next;

  prec_ptr->next  = ptr;
  return ptr;
}


/* Recursively free a list of subtype entries. */
void subtypeEntryDestroy(struct subtype_entry_t *entry)
{
  if (entry == NULL) return;
  subtypeEntryDestroy(entry->next);
  subtypeAttrDestroy(entry->atts);
  free(entry);
  entry = NULL;
}


/* Compares two subtype entries. */
static int subtypeEntryCompare(struct subtype_entry_t *e1, struct subtype_entry_t *e2)
{
  if (e1 == NULL)  Error("Internal error!");
  if (e2 == NULL)  Error("Internal error!");
  return 
    (e1->self == e2->self) && 
    subtypeAttsCompare(e1->atts, e2->atts);
}


/* (Recursively) duplicate list of entries. */
void subtypeEntryDuplicate(struct subtype_entry_t *a1, subtype_t* dst)
{
  if (a1 == NULL) return;
  /* append entry to dst pointer */
  struct subtype_entry_t *ptr = subtypeEntryAppend(dst);
  /* duplicate attributes */
  subtypeAttsDuplicate(a1->atts, ptr);
  ptr->self = a1->self;
  /* call next link in linked list */
  subtypeEntryDuplicate(a1->next, dst);
}



/* ------------------------------------------------------------------- */
/* SUBROUTINES FOR THE SUBTYPE ITSELF				       */
/* ------------------------------------------------------------------- */

/* Print-out subtype data structure together with its attributes. */
void subtypePrintKernel(subtype_t *subtype_ptr, FILE *fp)
{
  if (subtype_ptr == NULL)  Error("Internal error!");
  fprintf(fp, "# %s (subtype ID %d)\n", subtypeName[subtype_ptr->subtype], subtype_ptr->self);
  /* print global attributes of this subtype */
  struct subtype_attr_t* ptr = subtype_ptr->globals.atts;
  if (ptr != NULL)  fprintf(fp, "#\n# global attributes:\n");
  while (ptr != NULL) {
    fprintf(fp, "#   %-40s   (%2d) : %d\n", cdiSubtypeAttributeName[ptr->key], ptr->key, ptr->val);    
    ptr = ptr->next;
  } 
  /* print attributes for each subtype */
  fprintf(fp, "# %d local entries:\n", subtype_ptr->nentries);
  struct subtype_entry_t *entry = subtype_ptr->entries;
  while (entry != NULL) {
    fprintf(fp, "# subtype entry %d\n", entry->self);
    ptr = entry->atts;
    if (ptr != NULL)  fprintf(fp, "#   attributes:\n");
    while (ptr != NULL) {
      fprintf(fp, "#     %-40s (%2d) : %d\n", cdiSubtypeAttributeName[ptr->key], ptr->key, ptr->val);    
      ptr = ptr->next;
    } 
    entry = entry->next;
  }
  fprintf(fp, "\n");
}


/* Compares two subtype data structures. Pointer version of this
   method. */
static int subtypeCompareP(subtype_t *s1, subtype_t *s2)
{
  xassert(s1 && s2);
  if (s1->subtype != s2->subtype) return differ;
  if (subtypeEntryCompare(&s1->globals, &s2->globals) != 0) return differ;

  struct subtype_entry_t *entry1 = s1->entries;
  struct subtype_entry_t *entry2 = s2->entries;
  while ((entry1 != NULL) && (entry2 != NULL)) {
    if (subtypeEntryCompare(entry1, entry2) != 0)  return differ;
    entry1 = entry1->next;
    entry2 = entry2->next;
  }
  /* compare list lengths: */
  if ((entry1 != NULL) || (entry2 != NULL))  return differ;
  return 0;
}


/* Clean up data structure. */
static void subtypeDestroyP(void *ptr)
{
  subtype_t *subtype_ptr = (subtype_t*) ptr;
  /* destroy global attributes */
  subtypeAttrDestroy(subtype_ptr->globals.atts);
  /* destroy list of subtype entries */
  subtypeEntryDestroy(subtype_ptr->entries);
  subtype_ptr->entries = NULL;
  free(subtype_ptr);
  subtype_ptr = NULL;
}


/* Non-static wrapper function for "subtypeDestroyP". */
void subtypeDestroyPtr(void *ptr)
{
  subtypeDestroyP(ptr);
}


/* Non-static wrapper function for "subtypeCompareP". */
int subtypeComparePtr(int s1_ID, subtype_t *s2)
{
  subtype_t *subtype_ptr = reshGetVal(s1_ID, &subtypeOps);
  if (subtype_ptr == NULL)  Error("Internal error");
  return subtypeCompareP(subtype_ptr,s2);
}


/* Print-out subtype data structure together with its attributes.
   Pointer version of this method. */
static void subtypePrintP(void * subtype_ptr, FILE * fp)
{  subtypePrintKernel(subtype_ptr, fp); }



/* Print-out subtype data structure together with its attributes. */
void subtypePrintPtr(subtype_t* subtype_ptr)
{
  subtypePrintKernel(subtype_ptr, stdout);
}


/* Fill subtype data structure with default values. */
static void subtypeDefaultValue(subtype_t *subtype_ptr)
{
  if (subtype_ptr == NULL)  Error("Internal error!");
  subtype_ptr->self                 = CDI_UNDEFID;
  subtype_ptr->nentries             = 0;
  subtype_ptr->entries              = NULL;
  subtype_ptr->globals.atts         = NULL;
  subtype_ptr->globals.next         = NULL;
  subtype_ptr->globals.self         = -1;
  subtype_ptr->active_subtype_index = 0;
}


void subtypeAllocate(subtype_t **subtype_ptr2, int subtype)
{
  /* allocate new subtype */
  (*subtype_ptr2) = (subtype_t *) malloc(sizeof(subtype_t));
  subtype_t* subtype_ptr = *subtype_ptr2;
  subtypeDefaultValue(subtype_ptr);
  subtype_ptr->subtype = subtype;
  subtype_ptr->self    = CDI_UNDEFID;
}


/* Create a copy of an existing subtype data structure. */
void subtypeDuplicate(subtype_t *subtype_ptr, subtype_t **dst_ptr)
{
  if (subtype_ptr == NULL)  Error("Internal error!");
  subtypeAllocate(dst_ptr, subtype_ptr->subtype);
  subtype_t *dst = (*dst_ptr);
  /* create duplicate of subtype globals */
  subtypeAttsDuplicate(subtype_ptr->globals.atts, &dst->globals);
  dst->globals.self = subtype_ptr->globals.self;
  /* create duplicate of subtype entries */
  subtypeEntryDuplicate( subtype_ptr->entries, dst);
}


/* Register subtype object at resource handler. */
int subtypePush(subtype_t *subtype_ptr)
{
  if (subtype_ptr == NULL)  Error("Internal error!");
  subtype_ptr->self = reshPut(subtype_ptr, &subtypeOps);
  return subtype_ptr->self; /* subtypeID */
}



/* Sets an attribute for a subtype (for example a set of TILES). If
   the attribute has already been defined, then its value is
   overwritten. */
void subtypeDefGlobalDataP(subtype_t *subtype_ptr, int key, int val)
{
  if (subtype_ptr == NULL)  Error("Internal error!");
  /* find entry in linked list or append otherwise */
  struct subtype_attr_t* att_ptr = subtypeAttrFind(subtype_ptr->globals.atts, key);
  if (att_ptr == NULL) 
    subtypeAttrInsert(&subtype_ptr->globals, key, val);
  else
    att_ptr->val = val;
}


/* Sets an attribute for a subtype (for example a set of TILES). If
   the attribute has already been defined, then its value is
   overwritten. */
void subtypeDefGlobalData(int subtypeID, int key, int val)
{
  subtype_t *subtype_ptr = reshGetVal(subtypeID, &subtypeOps);
  subtypeDefGlobalDataP(subtype_ptr, key, val);
}


/* Retrieves an attribute for a subtype (for example a set of TILES).
   If the attribute has not been defined, then return -1. */
int subtypeGetGlobalDataP(subtype_t *subtype_ptr, int key)
{
  if (subtype_ptr == NULL)  Error("Internal error!");
  /* find entry in linked list */
  struct subtype_attr_t* att_ptr = subtypeAttrFind(subtype_ptr->globals.atts, key);
  if (att_ptr == NULL) 
    return -1;
  else
    return att_ptr->val;
}


/* Retrieves an attribute for a subtype (for example a set of TILES) .
   If the attribute has not been defined, then return -1. */
int subtypeGetGlobalData(int subtypeID, int key)
{
  subtype_t *subtype_ptr = reshGetVal(subtypeID, &subtypeOps);
  return subtypeGetGlobalDataP(subtype_ptr, key);
}


/* Sets an attribute for a single subtype entry (e.g. a single TILE).
   If the attribute has already been defined, then its value is
   overwritten. */
void subtypeDefEntryDataP(struct subtype_entry_t *subtype_entry_ptr, int key, int val)
{
  if (subtype_entry_ptr == NULL)  Error("Internal error!");
  /* find entry in linked list or append otherwise */
  struct subtype_attr_t* att_ptr = subtypeAttrFind(subtype_entry_ptr->atts, key);
  if (att_ptr == NULL) 
    subtypeAttrInsert(subtype_entry_ptr, key, val);
  else
    att_ptr->val = val;
}



/* ------------------------------------------------------------------- */
/* IMPLEMENTATIONS FOR KEY-VALUE-PAIR QUERIES			       */
/* ------------------------------------------------------------------- */


/* Generate a "query object" out of a key-value pair. */
subtype_query_t keyValuePair(const char* key, int value)
{
  subtype_query_t result;
  result.nAND = 1;
  result.key_value_pairs[0][0] = attribute_to_index(key);
  result.key_value_pairs[1][0] = value;
  if (CDI_Debug) {
    Message("key  %s matches %d", key, result.key_value_pairs[0][0]);
    Message("%d --?-- %d", result.key_value_pairs[0][0], result.key_value_pairs[1][0]);
  }
  return result;
}


/* Generate an AND-combined "query object" out of two previous query
   objects. */
subtype_query_t matchAND(subtype_query_t q1, subtype_query_t q2)
{
  if ((q1.nAND + q2.nAND) > MAX_KV_PAIRS_MATCH)  Error("Internal error");
  subtype_query_t result;
  result.nAND = q1.nAND;
  for (int i=0; i<q1.nAND; i++)
    {
      result.key_value_pairs[0][i] = q1.key_value_pairs[0][i];
      result.key_value_pairs[1][i] = q1.key_value_pairs[1][i];
    }
  for (int i=0; i<q2.nAND; i++)
    {
      result.key_value_pairs[0][result.nAND] = q2.key_value_pairs[0][i];
      result.key_value_pairs[1][result.nAND] = q2.key_value_pairs[1][i];
      result.nAND++;
    }

  if (CDI_Debug) {
    Message("combined criterion:");
    for (int i=0; i<result.nAND; i++)
      Message("%d --?-- %d", result.key_value_pairs[0][i], result.key_value_pairs[1][i]);
  }
  return result;
}



/* ------------------------------------------------------------------- */
/* SPECIFIC IMPLEMENTATIONS FOR TILE SETS			       */
/* ------------------------------------------------------------------- */


/* Integrate tile set "s2" into the tile set "subtype1_ID":

   Insert all entries set 2 to set 1 together with its attributes.
*/
void tilesetInsertP(subtype_t *s1, subtype_t *s2)
{
  if (s1 == NULL)  Error("Internal error!");
  if (s2 == NULL)  Error("Internal error!");
  struct subtype_entry_t 
    *entry1 = s1->entries,
    *entry2 = s2->entries;
  struct subtype_attr_t *att_ptr2;

  /* test all entries of set 2 against set 1, to check if entry
     already exists: */
  if (subtypeAttsCompare(s1->globals.atts, s2->globals.atts) != differ) 
    {
      while (entry1 != NULL) {
	int found = 1;
	entry2 = s2->entries;
	while (entry2 != NULL) {
	  found &= (subtypeAttsCompare(entry1->atts, entry2->atts) != differ);
	  entry2 = entry2->next;
	}
	if (found) 
	  {
	    return;
	  }
	entry1 = entry1->next;
      }
    
      entry2 = s2->entries;
      while (entry2 != NULL) {
	entry1 = subtypeEntryInsert(s1);
	
	att_ptr2 = entry2->atts;
	while (att_ptr2 != NULL) {
	  (void) subtypeAttrInsert(entry1, att_ptr2->key, att_ptr2->val);
	  att_ptr2 = att_ptr2->next;
	}
	entry2 = entry2->next;
      }
    }
  else
    {
      fprintf(stderr, "\n# SUBTYPE A:\n");
      subtypePrintKernel(s1, stderr);
      fprintf(stderr, "\n# SUBTYPE B:\n");
      subtypePrintKernel(s2, stderr);
      Error("Attempting to insert subtype entry into subtype with different global attributes!");
    }
}



/* ------------------------------------------------------------------- */
/* IMPLEMENTATIONS FOR ROUTINES VISIBLE THROUGH CDI.H		       */
/* ------------------------------------------------------------------- */


/*
  @Function  subtypeCreate
  @Title     Create a variable subtype
  
  @Prototype int subtypeCreate(int subtype)
  @Parameter
  @Item  subtype  The type of the variable subtype, one of the set of predefined CDI variable subtypes.
  The valid CDI variable subtypes are @func{SUBTYPE_TILES}
  
  @Description
  The function @func{subtypeCreate} creates a variable subtype.
  
  @Result
  @func{subtypeCreate} returns an identifier to the variable subtype.
  
  @EndFunction
*/
int subtypeCreate(int subtype)
{
  if ( CDI_Debug )  Message("subtype: %d ", subtype);
  Message("subtype: %d ", subtype);

  /* allocate new subtype */
  subtype_t *subtype_ptr;
  subtypeAllocate(&subtype_ptr, subtype);
  /* register object at resource handler */
  return subtypePush(subtype_ptr);
}


/* Print-out subtype data structure together with its attributes. */
void subtypePrint(int subtypeID)
{
  subtype_t *subtype_ptr = reshGetVal(subtypeID, &subtypeOps);
  subtypePrintKernel(subtype_ptr, stdout);
}


/* Compares two subtype data structures. */
int subtypeCompare(int subtypeID1, int subtypeID2)
{
  subtype_t *subtype_ptr1 = reshGetVal(subtypeID1, &subtypeOps);
  subtype_t *subtype_ptr2 = reshGetVal(subtypeID2, &subtypeOps);
  return subtypeCompareP(subtype_ptr1,subtype_ptr2);
}


/*  Get the size of a subtype (e.g. no. of tiles). */
int subtypeInqSize(int subtypeID)
{
  if ( subtypeID == CDI_UNDEFID )
    {
      return 0;
    }
  else
    {
      subtype_t *subtype_ptr = reshGetVal(subtypeID, &subtypeOps);
      return subtype_ptr->nentries;
    }
}


/* Get the currently active index of a subtype (e.g. current tile index). */
int subtypeInqActiveIndex(int subtypeID)
{
  if (subtypeID == CDI_UNDEFID)  return 0;
  subtype_t *subtype_ptr = reshGetVal(subtypeID, &subtypeOps);
  return subtype_ptr->active_subtype_index;
}


/* Set the currently active index of a subtype (e.g. current tile index). */
void subtypeDefActiveIndex(int subtypeID, int index)
{
  subtype_t *subtype_ptr = reshGetVal(subtypeID, &subtypeOps);
  subtype_ptr->active_subtype_index = index;
}


/* subtypeInqSubEntry: Returns subtype entry ID for a given
   criterion. */
int subtypeInqSubEntry(int subtypeID, subtype_query_t criterion)
{
  subtype_t *subtype_ptr = reshGetVal(subtypeID, &subtypeOps);
  struct subtype_entry_t *entry = subtype_ptr->entries;
  /* loop over all entries of this subtype */
  while (entry != NULL) {
    {
      int match = 1;
      /* test if this entry matches ALL criteria. */
      for (int j=0; (j<criterion.nAND) && (match); j++)
	{
	  if (CDI_Debug)  Message("check criterion %d :  %d --?-- %d", j,
				  criterion.key_value_pairs[0][j], criterion.key_value_pairs[1][j]);
	  struct subtype_attr_t* att_ptr = 
	    subtypeAttrFind(entry->atts, criterion.key_value_pairs[0][j]);
	  if (att_ptr == NULL)
	    {
	      match = 0;
	      if (CDI_Debug)  Message("did not find %d", criterion.key_value_pairs[0][j]);
	    }
	  else
	    {
	      if (CDI_Debug)  Message("found %d", criterion.key_value_pairs[0][j]);
	      match &= (att_ptr->val == criterion.key_value_pairs[1][j]);
	    }
	}
      if (match) return entry->self;
    }
    entry = entry->next;
  }
  return CDI_UNDEFID;
}


int subtypeInqTile(int subtypeID, int tileindex, int attribute)
{
  return subtypeInqSubEntry(subtypeID, 
			    matchAND(keyValuePair(cdiSubtypeAttributeName[SUBTYPE_ATT_TILEINDEX], tileindex),
				     keyValuePair(cdiSubtypeAttributeName[SUBTYPE_ATT_TILEATTRIBUTE], attribute)));
}


/* Construct a new subtype for a tile set. If a corresponding subtype
 * already exists, then we return this subtype ID instead. 
 *
 * See comment on subtype.c::tilesetMatchingPtr for the specification
 * of the term "corresponding" tile set.
 */
int vlistDefTileSubtype(int vlistID, subtype_t *tiles)
{
  int subtypeID = CDI_UNDEFID;

  /* loop over subtypes and search for an identical tileset */
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  int      tileset_defined = 0;
  for (int isub=0; isub<vlistptr->nsubtypes; isub++)
    {
      /* get the ID of the "isub"th subtype */
      subtypeID = vlistptr->subtypeIDs[isub];
      if (subtypeComparePtr(subtypeID, tiles) == 0)
        {
          tileset_defined = 1;
          break;
        }
    }

  /* tile set seems to be new: register at resource handler. */
  if (tileset_defined == 0)  {
    subtype_t *tiles_duplicate = NULL;
    subtypeDuplicate(tiles, &tiles_duplicate);
    subtypeID = vlistptr->subtypeIDs[vlistptr->nsubtypes++] = subtypePush(tiles_duplicate);
  }

  return subtypeID;
}



int vlistInsertTrivialTileSubtype(int vlistID)
{
  /* first, generate a subtype */
  subtype_t *subtype_ptr;
  subtypeAllocate(&subtype_ptr, SUBTYPE_TILES);

  /* create a tile set that contains only one tile/attribute pair. */
  (void) subtypeEntryInsert(subtype_ptr);

  /* register tile */
  vlist_t *vlistptr = vlist_to_pointer(vlistID);
  int subtypeID = vlistptr->subtypeIDs[vlistptr->nsubtypes++] = subtypePush(subtype_ptr);
  return subtypeID;
}




/* ------------------------------------------------------------------- */
/* NOT YET IMPLEMENTED						       */
/* ------------------------------------------------------------------- */

static int subtypeGetPackSize( void * subtype_ptr, void *context)
{  Error("Not yet implemented for subtypes!");  return 0; }

static void subtypePack( void * subtype_ptr, void * buffer, int size, int *pos, void *context)
{  Error("Not yet implemented for subtypes!"); }

static int subtypeTxCode( void )
{  Error("Not yet implemented for subtypes!");  return 0; }


