#if defined(HAVE_CONFIG_H)
#  include "config.h"
#endif

#if defined(_OPENMP)
#include <omp.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "cdo.h"
#include "remap.h"

/*****************************************************************************/

static
int isSorted(int *restrict array1, int *restrict array2, const long size)
{
  for ( long idx = 1; idx < size; ++idx )
    {
      if (  (array1[idx-1] >  array1[idx]) ||
	   ((array1[idx-1] == array1[idx]) &&
	    (array2[idx-1] >  array2[idx])) )
	{
	  return 0;
	}
    }

  return 1;
}
/*
static
void swap(int *restrict v1, int *restrict v2)
{
  const int tmp = *v1;
  *v1 = *v2;
  *v2 = tmp;
}

static
void heapify(const long pos, const long size, int *restrict arr1, int *restrict arr2, int *restrict arr3)
{
  const long leftChild = pos*2 + 1;
  
  if ( leftChild < size )
    {
      long maxIter = pos;
      if (  (arr1[leftChild] >  arr1[maxIter]) ||
	   ((arr1[leftChild] == arr1[maxIter]) &&
	    (arr2[leftChild] >  arr2[maxIter])) )
	{
	  maxIter = leftChild;
	}

      const long rightChild = leftChild + 1;
      if ( rightChild < size )
	if (  (arr1[rightChild] >  arr1[maxIter]) ||
	     ((arr1[rightChild] == arr1[maxIter]) &&
	      (arr2[rightChild] >  arr2[maxIter])) )
	  {
	    maxIter = rightChild;
	  }
 
      if ( maxIter != pos )
	{
	  swap(&arr1[maxIter], &arr1[pos]);
	  swap(&arr2[maxIter], &arr2[pos]);
	  swap(&arr3[maxIter], &arr3[pos]);
	  heapify(maxIter, size, arr1, arr2, arr3);
	}
    }
}
*/
/* remap_heapsort is 30% faster than remap_heapsort_recursiv ! */
/*
static
void remap_heapsort_recursiv(const long num_links, int *restrict arr1, int *restrict arr2, int *restrict arr3)
{
  long lvl;     // level indexes for heap sort levels

  for ( lvl = num_links/2-1; lvl >= 0; --lvl )
    {
      heapify(lvl, num_links, arr1, arr2, arr3);
    }

  for ( lvl = num_links-1; lvl > 0; --lvl )
    {
      swap(&arr1[0], &arr1[lvl]);
      swap(&arr2[0], &arr2[lvl]);
      swap(&arr3[0], &arr3[lvl]);
      heapify(0, lvl, arr1, arr2, arr3);
    }
}
*/

static
void remap_heapsort(const long num_links, int *restrict add1, int *restrict add2, int *restrict idx)
{
  int add1_tmp, add2_tmp;  /* temp for addresses during swap     */
  int idx_tmp;
  long lvl, final_lvl;     /* level indexes for heap sort levels */
  long chk_lvl1, max_lvl;
  long i;

  /*
    start at the lowest level (N/2) of the tree and shift lower 
    values to the bottom of the tree, promoting the larger numbers
  */
  for ( lvl = num_links/2-1; lvl >= 0; --lvl )
    {
      final_lvl = lvl;
      add1_tmp = add1[lvl];
      add2_tmp = add2[lvl];
      idx_tmp  = idx[lvl];

      /* Loop until proper level is found for this link, or reach bottom */

      for ( i = 0; i < num_links; ++i ) /* while ( TRUE ) */
	{
	  /* Find the largest of the two daughters */

          chk_lvl1 = 2*final_lvl+1;
          max_lvl  = 2*final_lvl+2;
          if ( chk_lvl1 == num_links-1 ) max_lvl = chk_lvl1;

          if ( (add1[chk_lvl1] >  add1[max_lvl]) ||
	      ((add1[chk_lvl1] == add1[max_lvl]) &&
               (add2[chk_lvl1] >  add2[max_lvl])) )
	    {
	      max_lvl = chk_lvl1;
	    }

          /*
	    If the parent is greater than both daughters,
	    the correct level has been found
	  */
          if ( (add1_tmp >  add1[max_lvl]) ||
              ((add1_tmp == add1[max_lvl]) &&
               (add2_tmp >  add2[max_lvl])) )
	    {
	      add1[final_lvl] = add1_tmp;
	      add2[final_lvl] = add2_tmp;
	      idx[final_lvl]  = idx_tmp;

	      break;
	    }
	  else
	    {
	      /*
		Otherwise, promote the largest daughter and push
		down one level in the tree.  If haven't reached
		the end of the tree, repeat the process.  Otherwise
		store last values and exit the loop
	      */
	      add1[final_lvl] = add1[max_lvl];
	      add2[final_lvl] = add2[max_lvl];
	      idx[final_lvl]  = idx[max_lvl];

	      final_lvl = max_lvl;
	      if ( 2*final_lvl+1 >= num_links )
		{
		  add1[final_lvl] = add1_tmp;
		  add2[final_lvl] = add2_tmp;
		  idx[final_lvl]  = idx_tmp;

		  break;
		}
	    }
	}

      if ( i == num_links ) cdoAbort("Internal problem, link 1 not found!");
    }

  /*
    Now that the heap has been sorted, strip off the top (largest)
    value and promote the values below
  */
  for ( lvl = num_links-1; lvl >= 2; --lvl )
    {
      /* Move the top value and insert it into the correct place */

      add1_tmp  = add1[lvl];
      add1[lvl] = add1[0];

      add2_tmp  = add2[lvl];
      add2[lvl] = add2[0];

      idx_tmp  = idx[lvl];
      idx[lvl] = idx[0];

      /* As above this loop sifts the tmp values down until proper level is reached */

      final_lvl = 0;

      for ( i = 0; i < num_links; ++i ) /* while ( TRUE ) */
	{
	  /* Find the largest of the two daughters */

          chk_lvl1 = 2*final_lvl+1;
          max_lvl  = 2*final_lvl+2;
          if ( max_lvl >= lvl ) max_lvl = chk_lvl1;

          if ( (add1[chk_lvl1] >  add1[max_lvl]) ||
              ((add1[chk_lvl1] == add1[max_lvl]) &&
               (add2[chk_lvl1] >  add2[max_lvl])) )
	    {
	      max_lvl = chk_lvl1;
	    }

          /*
	    If the parent is greater than both daughters,
	    the correct level has been found
	  */
          if ( (add1_tmp >  add1[max_lvl]) ||
              ((add1_tmp == add1[max_lvl]) &&
               (add2_tmp >  add2[max_lvl])) )
	    {
	      add1[final_lvl] = add1_tmp;
	      add2[final_lvl] = add2_tmp;
	      idx[final_lvl]  = idx_tmp;

	      break;
	    }
	  else
	    {
	      /*
		Otherwise, promote the largest daughter and push
		down one level in the tree.  If haven't reached
		the end of the tree, repeat the process.  Otherwise
		store last values and exit the loop
	      */
	      add1[final_lvl] = add1[max_lvl];
	      add2[final_lvl] = add2[max_lvl];
	      idx[final_lvl]  = idx[max_lvl];

	      final_lvl = max_lvl;
	      if ( 2*final_lvl+1 >= lvl )
		{
		  add1[final_lvl] = add1_tmp;
		  add2[final_lvl] = add2_tmp;
		  idx[final_lvl]  = idx_tmp;

		  break;
		}
	    }
	}

      if ( i == num_links ) cdoAbort("Internal problem, link 2 not found!");
    }

  /* Swap the last two entries */

  add1_tmp = add1[1];
  add1[1]  = add1[0];
  add1[0]  = add1_tmp;

  add2_tmp = add2[1];
  add2[1]  = add2[0];
  add2[0]  = add2_tmp;

  idx_tmp  = idx[1];
  idx[1]   = idx[0];
  idx[0]   = idx_tmp;
}


void sort_add(long num_links, long num_wts, int *restrict add1, int *restrict add2, double *restrict weights)
{
  /*
    This routine sorts address and weight arrays based on the destination address with the 
    source address as a secondary sorting criterion. The method is a standard heap sort.
  */
  /*
    Input and Output arrays:
    
       long num_links; ! num of links for this mapping
       long num_wts;   ! num of weights for this mapping
       add1,           ! destination address array [num_links]
       add2            ! source      address array
       weights         ! remapping weights [num_links*num_wts]
  */

  /* Local variables */

  int *idx;
  long i, n;
  double *wgt_tmp;

  if ( num_links <= 1 ) return;

  idx = (int*) malloc(num_links*sizeof(int));
  for ( i = 0; i < num_links; ++i ) idx[i] = i;

  remap_heapsort(num_links, add1, add2, idx);

  wgt_tmp = (double*) malloc(num_wts*num_links*sizeof(double));
  memcpy(wgt_tmp, weights, num_wts*num_links*sizeof(double));

  for ( i = 0; i < num_links; ++i )
    for ( n = 0; n < num_wts; ++n )
      weights[num_wts*i+n] = wgt_tmp[num_wts*idx[i]+n];

  free(wgt_tmp);
  free(idx);

  if ( cdoVerbose ) 
    if ( !isSorted(add1, add2, num_links) ) fprintf(stderr, ">>>> sort_add failed!!!\n");
} /* sort_add */

/*****************************************************************************/

void sort_add_orig(long num_links, long num_wts, int *restrict add1, int *restrict add2, double *restrict weights)
{
  /*
    This routine sorts address and weight arrays based on the
    destination address with the source address as a secondary
    sorting criterion. The method is a standard heap sort.
  */
  /*
    Input and Output arrays:
    
       long num_links; ! num of links for this mapping
       long num_wts;   ! num of weights for this mapping
       add1,           ! destination address array [num_links]
       add2            ! source      address array
       weights         ! remapping weights [num_links*num_wts]
  */

  /* Local variables */

  int add1_tmp, add2_tmp;  /* temp for addresses during swap     */
  long lvl, final_lvl;     /* level indexes for heap sort levels */
  long chk_lvl1, chk_lvl2, max_lvl;
  long i, n;
  double wgttmp[4];        /* temp for holding wts during swap   */

  if ( num_links <= 1 ) return;

  /*
  for ( n = 0; n < num_links; n++ )
    printf("in: %5d %5d %5d # dst_add src_add n\n", add1[n]+1, add2[n]+1, n+1);
  */
  /*
    start at the lowest level (N/2) of the tree and shift lower 
    values to the bottom of the tree, promoting the larger numbers
  */
  for ( lvl = num_links/2-1; lvl >= 0; lvl-- )
    {
      final_lvl = lvl;
      add1_tmp = add1[lvl];
      add2_tmp = add2[lvl];
      for ( n = 0; n < num_wts; n++ )
	wgttmp[n] = weights[num_wts*lvl+n];

      /* Loop until proper level is found for this link, or reach bottom */

      for ( i = 0; i < num_links; i++ ) /* while ( TRUE ) */
	{
	  /* Find the largest of the two daughters */

          chk_lvl1 = 2*final_lvl+1;
          chk_lvl2 = 2*final_lvl+2;
          if ( chk_lvl1 == num_links-1 ) chk_lvl2 = chk_lvl1;

          if ((add1[chk_lvl1] >  add1[chk_lvl2]) ||
	     ((add1[chk_lvl1] == add1[chk_lvl2]) &&
              (add2[chk_lvl1] >  add2[chk_lvl2])))
            max_lvl = chk_lvl1;
          else 
            max_lvl = chk_lvl2;

          /*
	    If the parent is greater than both daughters,
	    the correct level has been found
	  */
          if ((add1_tmp >  add1[max_lvl]) ||
             ((add1_tmp == add1[max_lvl]) &&
              (add2_tmp >  add2[max_lvl])))
	    {
	      add1[final_lvl] = add1_tmp;
	      add2[final_lvl] = add2_tmp;
	      for ( n = 0; n < num_wts; n++ )
		weights[num_wts*final_lvl+n] = wgttmp[n];

	      break;
	    }
	  else
	    {
	      /*
		Otherwise, promote the largest daughter and push
		down one level in the tree.  If haven"t reached
		the end of the tree, repeat the process.  Otherwise
		store last values and exit the loop
	      */
	      add1[final_lvl] = add1[max_lvl];
	      add2[final_lvl] = add2[max_lvl];
	      for ( n = 0; n < num_wts; n++ )
		weights[num_wts*final_lvl+n] = weights[num_wts*max_lvl+n];

	      final_lvl = max_lvl;
	      if ( 2*final_lvl+1 >= num_links )
		{
		  add1[final_lvl] = add1_tmp;
		  add2[final_lvl] = add2_tmp;
		  for ( n = 0; n < num_wts; n++ )
		    weights[num_wts*final_lvl+n] = wgttmp[n];

		  break;
		}
	    }
	}

      if ( i == num_links ) cdoAbort("Internal problem, link 1 not found!");
    }

  /*
    Now that the heap has been sorted, strip off the top (largest)
    value and promote the values below
  */
  for ( lvl = num_links-1; lvl >= 2; lvl-- )
    {
      /* Move the top value and insert it into the correct place */

      add1_tmp = add1[lvl];
      add1[lvl] = add1[0];

      add2_tmp = add2[lvl];
      add2[lvl] = add2[0];

      for ( n = 0; n < num_wts; n++ )
        wgttmp[n] = weights[num_wts*lvl+n];

      for ( n = 0; n < num_wts; n++ )
        weights[num_wts*lvl+n] = weights[n];

      /* As above this loop sifts the tmp values down until proper level is reached */

      final_lvl = 0;

      for ( i = 0; i < num_links; i++ ) /* while ( TRUE ) */
	{
	  /* Find the largest of the two daughters */

          chk_lvl1 = 2*final_lvl+1;
          chk_lvl2 = 2*final_lvl+2;
          if ( chk_lvl2 >= lvl ) chk_lvl2 = chk_lvl1;

          if ((add1[chk_lvl1] >  add1[chk_lvl2]) ||
             ((add1[chk_lvl1] == add1[chk_lvl2]) &&
              (add2[chk_lvl1] >  add2[chk_lvl2])))
            max_lvl = chk_lvl1;
          else 
            max_lvl = chk_lvl2;

          /*
	    If the parent is greater than both daughters,
	    the correct level has been found
	  */
          if ((add1_tmp >  add1[max_lvl]) ||
             ((add1_tmp == add1[max_lvl]) &&
              (add2_tmp >  add2[max_lvl])))
	    {
	      add1[final_lvl] = add1_tmp;
	      add2[final_lvl] = add2_tmp;
	      for ( n = 0; n < num_wts; n++ )
		weights[num_wts*final_lvl+n] = wgttmp[n];

	      break;
	    }
	  else
	    {
	      /*
		Otherwise, promote the largest daughter and push
		down one level in the tree.  If haven't reached
		the end of the tree, repeat the process.  Otherwise
		store last values and exit the loop
	      */
	      add1[final_lvl] = add1[max_lvl];
	      add2[final_lvl] = add2[max_lvl];
	      for ( n = 0; n < num_wts; n++ )
		weights[num_wts*final_lvl+n] = weights[num_wts*max_lvl+n];

	      final_lvl = max_lvl;
	      if ( 2*final_lvl+1 >= lvl )
		{
		  add1[final_lvl] = add1_tmp;
		  add2[final_lvl] = add2_tmp;
		  for ( n = 0; n < num_wts; n++ )
		    weights[num_wts*final_lvl+n] = wgttmp[n];

		  break;
		}
	    }
	}

      if ( i == num_links ) cdoAbort("Internal problem, link 2 not found!");
    }

  /* Swap the last two entries */

  add1_tmp = add1[1];
  add1[1]  = add1[0];
  add1[0]  = add1_tmp;

  add2_tmp = add2[1];
  add2[1]  = add2[0];
  add2[0]  = add2_tmp;

  for ( n = 0; n < num_wts; n++ ) wgttmp[n]          = weights[num_wts+n];
  for ( n = 0; n < num_wts; n++ ) weights[num_wts+n] = weights[n];
  for ( n = 0; n < num_wts; n++ ) weights[n]         = wgttmp[n];
  /*
  for ( n = 0; n < num_links; n++ )
    printf("out: %5d %5d %5d # dst_add src_add n\n", add1[n]+1, add2[n]+1, n+1);
  */
} /* sort_add_orig */


/* ******************************************************************************** 
     XXX       XXX    XXXXXXXXXX    XXXXXXXXXX       XXXXXXXXXXX   XXXXXXXXXX
    XXXX     XXXX    XXXXXXXXXX    XXXXXXXXXXX     XXXXXXXXXXX    XXXXXXXXXX
    XXXXX   XXXXX    XXX           XXX     XXX     XXX            XXX
    XXXXXX XXXXXX    XXXXXXXXX     XXXX  XXXX      XXX            XXXXXXXX
    XXX  XXX  XXX    XXXXXXXXX     XXXXXXX         XXX     XXXX   XXXXXXXX
    XXX   X   XXX    XXX           XXX XXXX        XXX      XXX   XXX
    XXX       XXX    XXXXXXXXXX    XXX   XXXX      XXXXXXXXXXXX   XXXXXXXXXX
    XXX       XXX    XXXXXXXXXX    XXX     XXXX     XXXXXXXXXX    XXXXXXXXXX


           XXXXXXXXXXX      XXXXXXXX     XXXXXXXXXX      XXXXXXXXXXXXX
          XXXXXXXXXXX      XXX    XXX    XXXXXXXXXXX     XXXXXXXXXXXXX
          XXX             XXX      XXX   XXX     XXX          XXX
          XXXXXXXXXXX     XXX      XXX   XXXX  XXXX           XXX
           XXXXXXXXXXX    XXX      XXX   XXXXXXXX             XXX
                   XXX    XXX      XXX   XXX XXXX             XXX
           XXXXXXXXXXX     XXX    XXX    XXX   XXXX           XXX
          XXXXXXXXXXX       XXXXXXXX     XXX     XXXX         XXX
********************************************************************************** */

/* MERGE SORT DEFINES */
//#define MERGE_SORT_CHUNKS          64
#define MERGE_SORT_LIMIT_SIZE      4096 //num_links/(MERGE_SORT_CHUNKS*omp_get_num_procs())


static
void merge_lists(int *nl, int *l11, int *l12, int *l21, int *l22, long *idx)
{      
  /*
    This routine writes to idx a list of indices relative to *l11 and *l12
    --> l11, l12, and l21,l22 each need to be allocated in 
        a signle block of memory
    The order is thus, that (I) l11[idx[i]]<l11[idx[i+1]]	
                        OR (II) l11[idx[i]]==l11[idx[i+1]] && l21[idx[i]]<l21[idx[i+1]]
		       where 0 <= i < nl
  */    		       
  int i1=0, i2=0, i=0, ii;
  const int n1=nl[0], n2=nl[1];

  i=0;
  while ( i2 < n2 && i1 < n1 ) 
    {
      if ( ( l11[i1] < l21[i2] ) ||
	   ( l11[i1] == l21[i2] && l12[i1] < l22[i2] ) )
	{ idx[i] = i1;    i1++; }
      else
	{ idx[i] = n1+i2; i2++; }
      i++;
    }

  for ( ii=i1; i1 < n1; ii++ ) {idx[i] = i1;    i++; i1++; }
  for ( ii=i2; i2 < n2; ii++ ) {idx[i] = n1+i2; i++; i2++; }
}

static
void sort_par(long num_links, long num_wts, int *restrict add1, int *restrict add2, 
	      double *restrict weights, int parent, int par_depth)
{
  /*
    This routine is the core of merge-sort. It does the following
     + split the address-arrays into two segments, 
     + sort each array seperately (this can be done in parallel as there
       is no data dependency)
       - the routine sort_iter, which is called for sorting the sub-arrays
         EITHER calls this routine againg, which means, that the sub-arrays
	 are further split 
	 OR     it calls sort_add, which actually sorts the sublist sequentially
     + merge the sorted arrays together
     For the merge step additional memory is needed as it cannot work in place. 
     Therefor, the merge sort algorith in this implementation uses at maximum 
     twice as much memory as the sequential sort_add.

     Parameters:
     -----------
       long num_links    | length of arrays add1 and add2 (MUST be of same length
       int *add1 *add2   | arrays with addresses, that are used as sorting criteria (ascending)
       double ** weights | weights for each address that have to be kept in the same order 
                           as add1[] and add2[]
       int parent        | the parent of this sort_par. This parameter is used to find 
                           the recursion depth and determine the actual position of the
			   sub-array within the original array 
  */


  const int nsplit = 2;                      /* (only 2 allowed) number of segments to split the data */
  int nl[nsplit];                            /* number of links in each sub-array              */
  int who_am_i,depth;                        /* current depth, depth of children and index
						to be parent in next call to sort_par          */
  int add_srt[nsplit]/*, add_end[nsplit]*/;  /* arrays for start and end index of sub array    */
  int *add1s[nsplit], *add2s[nsplit];        /* pointers to sub arrays for sort and merge step */
  int *tmp;                                  /* pointer to buffer for merging of address lists */
  double *tmp2 = NULL;                       /* pointer to buffer for merging weight lists     */
  double *wgttmp = NULL;                     /* pointer to buffer for swap weights             */
  long *idx;                                 /* index list to merge sub-arrays                 */
  long i,n,m;   

  // printf("sort_par: parent = %d numlinks = %ld\n", parent, num_links);
  if ( nsplit != 2 )
    {
      fprintf(stderr,"Error: splitting into more than two subsegments not allowed\n"
	     "       in this implementation of merge sort\n");
      exit(-1);
    }

  idx = (long*) malloc(num_links*sizeof(long));

  /* SPLIT AND SORT THE DATA FRAGMENTS */
  /*
  for ( i=0; i<nsplit; i++)
    {
      add_srt[i]= i * num_links/nsplit;
      add_end[i]= (i+1) * num_links/nsplit;
      add1s[i]  = &(add1[add_srt[i]]);
      add2s[i]  = &(add2[add_srt[i]]);
      nl[i]     = add_end[i]-add_srt[i];
    }
  */
  add_srt[0] = 0;                  add_srt[1] = num_links/nsplit;
  add1s[0]   = &add1[add_srt[0]];  add1s[1]   = &add1[add_srt[1]];
  add2s[0]   = &add2[add_srt[0]];  add2s[1]   = &add2[add_srt[1]];
  nl[0]      = num_links/nsplit;   nl[1]      = num_links-nl[0];
  //add_end[0] = nl[0];              add_end[1] = num_links;

  depth = (int) (log(parent)/log(2));

#if defined(_OPENMP)
  /* Allow for nested parallelism */
  if ( omp_in_parallel() && depth<par_depth ) 
    {
      omp_set_nested(1);            
      if ( omp_get_nested() == 0 )
	printf("Warning: OpenMP implementation seems to not support nested parallelism.\n"
	       "Maximum of CPUs used is 2 instead of %i.\n", ompNumThreads);
    }                                    
#endif

  //  printf("I am %i nl[0] %i nl[1] %i\n",parent,nl[0],nl[1]);
  //  printf("add_srt[0] %i add_Srt[1] %i\n",add_srt[0],add_srt[1]);
  //  if ( 1 )
  //      printf("\n\nSplitting thread into %i!! (I AM %i) depth %i parallel_depth %i add_srt[0]%i add_srt[1] %i\n",
  //	     nsplit,parent,depth,par_depth,add_srt[0],add_srt[1]);

#if defined(_OPENMP)
#pragma omp parallel for if ( depth < par_depth ) \
        private(i, n, m, wgttmp, who_am_i)   \
        shared(weights) num_threads(2)
#endif
  for ( i=0; i < nsplit; i++ )
    {

      who_am_i = nsplit*parent+i;
      //    my_depth = (int) (log(parent)/log(2))+1;

#if defined(_OPENMP)
      //      if ( 1 )
      //	printf("I am %i (parent %i), my_depth is: %i thread_num %i (%i) \n",
      //	       who_am_i,parent,my_depth,omp_get_thread_num()+1,omp_get_num_threads());
#endif
            
      wgttmp = (double*) malloc(num_wts*nl[i]*sizeof(double));
       
      for ( m = 0; m < nl[i]; m++ )
	for ( n = 0; n < num_wts; n++ )                      
	  wgttmp[num_wts*m+n] = weights[num_wts*(add_srt[i]+m)+n];

      sort_iter(nl[i], num_wts, add1s[i], add2s[i], wgttmp, who_am_i);

      for ( m = 0; m < nl[i]; m++ )
	for ( n = 0; n < num_wts; n++ )
	  weights[num_wts*(add_srt[i]+m)+n] = wgttmp[num_wts*m+n];

      free(wgttmp);
    }

  /* ********************************* */
  /*              TO DO                */
  /* THIS BIT NEEDS TO BE PARALLELIZED */
  /* ********************************* */
  /* Idea I: one CPU merges top-down, the other one bottom-up */
                                                              /* ********************** */
  merge_lists(nl,add1s[0],add2s[0],add1s[1],add2s[1], idx);   /* MERGE THE SEGMENTS     */
                                                              /* ********************** */
  tmp = (int*) malloc(num_links*sizeof(int));
  
#if defined(_OPENMP)
#pragma omp parallel for if ( depth < par_depth ) private(i) num_threads(2)
#endif
  for ( i = 0; i < num_links; i++ )
    tmp[i] = add1[idx[i]];
  
#if defined(_OPENMP)
#pragma omp parallel for if ( depth < par_depth ) private(i) num_threads(2)
#endif
  for ( i = 0; i < num_links; i++ )
    {
      add1[i] = tmp[i];
      tmp[i] = add2[idx[i]];
    }
  
#if defined(_OPENMP)
#pragma omp parallel for if ( depth < par_depth ) private(i) num_threads(2)
#endif
  for ( i = 0; i < num_links; i++ )
    add2[i] = tmp[i];
  
  free(tmp);
  tmp = NULL;
  
  tmp2 = (double*) malloc(num_links*num_wts*sizeof(double));
  
#if defined(_OPENMP)
#pragma omp parallel for if ( depth < par_depth ) private(i,n) num_threads(2)
#endif
  for ( i = 0; i < num_links; i++ )
    for ( n = 0; n < num_wts; n++ )
      tmp2[num_wts*i + n] = weights[num_wts*idx[i]+n];
  
#if defined(_OPENMP)
#pragma omp parallel for if ( depth < par_depth ) private(i,n) num_threads(2)
#endif
  for ( i = 0; i < num_links; i++ )
    for ( n = 0; n < num_wts; n++ )
      weights[num_wts*i+n] = tmp2[num_wts*i+n];
  
  free(tmp2);
  tmp2 = NULL;
  
  free(idx);
}


void sort_iter(long num_links, long num_wts, int *restrict add1, int *restrict add2, double *restrict weights, int parent)
{
  /*
    This routine is an interface between the parallelized (merge-sort) 
    and the sequential sorting algorithm for addresses implemented in
    the library. 
    It iterates 1 level into the binary tree if the single data chunks
    to sort are larger than the maximum size prescribed. Otherwise, it
    just sorts the chunk using the sort_add routine as implemented 
    originally. 
    Note, that even on a single CPU, the merge sort algorithm can be
    considerably faster (up to about 30% for a reasonable chunk size)
  */

  /* Parameters as in sort_par 
     additional parameters 
     int mod;      (enum TPAR_MODE) determines wether tomake use of merge sort
     int parent;   !!! CAUTION !!!
                   + determines level and position of data chunk within 
                     the original heap (level = log_2(who_am_i)) if sort_iter(...) has not
		     been called before
		   + determines number of threads to use on first call of sort_iter(...)
  */
  static int first_sort_iter_call = 1;
  static int par_depth = 1;
  static int nthreads = 1;

  if ( first_sort_iter_call )
    {
      first_sort_iter_call = 0;
      nthreads = parent;
      par_depth = (int)(log(parent)/log(2));
      parent = 1;
    }

  // fprintf(stdout, "parent %d par_depth %d num_links %ld\n", parent, par_depth, num_links);

  if ( num_links > MERGE_SORT_LIMIT_SIZE && parent <= (nthreads-1) )
    {
      sort_par(num_links, num_wts, add1, add2, weights, parent, par_depth);
      if ( cdoVerbose ) fprintf(stderr, "sort_iter: Finished iteration parent %i\n", parent);
    }
  else
    {
      // printf("sort_add: parent %d, par_depth %d num_links %ld\n", parent, par_depth, num_links);
      sort_add(num_links, num_wts, add1, add2, weights);
    }
}
