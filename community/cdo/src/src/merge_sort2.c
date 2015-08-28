#include "merge_sort2.h"

static
void merge_lists(int *nl, const double *restrict l1, const double *restrict l2, long *idx)
{      
  /*
    This routine writes to idx a list of indices relative to *l11 and *l12
    --> l1, and l2 need to be allocated in 
        a signle block of memory
    The order is thus, that (I) l1[idx[i]]<=l1[idx[i+1]]	
		       where 0 <= i < nl
  */    		       

  const int n1=nl[0], n2=nl[1];
  int i1=0, i2=0, i=0, ii=0;

  while ( i2 < n2 && i1 < n1 ) 
    {
      if ( l1[i1] < l2[i2] )
	{ idx[i] = i1;    i1++; i++; }
      else
	{ idx[i] = n1+i2; i2++; i++; }
    }

  for ( ii=i1; i1 < n1; ii++ ) {idx[i] = i1;    i++; i1++; }
  for ( ii=i2; i2 < n2; ii++ ) {idx[i] = n1+i2; i++; i2++; }

  return;
}

static
void sort_par(long num_links, double *restrict add1, int parent, int par_depth)
{
  /*
    This routine is the core of merge-sort. It does the following
     + split the address-arrays into two segments, 
     + sort each array seperately (this can be done in parallel as there
       is no data dependency)
       - the routine sort_iter, which is called for sorting the sub-arrays
         EITHER calls this routine again, which means, that the sub-arrays
	 are further split 
	 OR     it calls sort_add, which actually sorts the sublist sequentially
     + merge the sorted arrays together
     For the merge step additional memory is needed as it cannot work in place. 
     Therefor, the merge sort algorith in this implementation uses at maximum 
     twice as much memory as the sequential sort_add.

     Parameters:
     -----------
       long num_links    | length of arrays add1 and add2 (MUST be of same length
       int *add1         | arrays with addresses, that are used as sorting criteria (ascending)
       int parent        | the parent of this sort_par. This parameter is used to find 
                           the recursion depth and determine the actual position of the
			   sub-array within the original array 
  */
#define NSPLIT 2
  int nsplit = NSPLIT;                      /* (only 2 allowed) number of segments to split 
						the data */
  int nl[NSPLIT];                            /* number of links in each sub-array              */
  int who_am_i,depth,my_depth;               /* current depth, depth of children and index
						to be parent in next call to sort_par          */
  int add_srt[NSPLIT], add_end[NSPLIT];      /* arrays for start and end index of sub array    */
  double *add1s[NSPLIT];                     /* pointers to sub arrays for sort and merge step */
  double *tmp;                               /* pointer to buffer for merging of address lists */
  long *idx;                                 /* index list to merge sub-arrays                 */
  long i;   

  if ( nsplit != 2 )
    {
      cdoAbort("Error: splitting into more than two subsegments not allowed\n"
	     "       in this implementation of merge sort\n");
    }

  idx = (long*) malloc(num_links*sizeof(long));

  /* SPLIT AND SORT THE DATA FRAGMENTS */
  add_srt[0] = 0;                  add_srt[1] = num_links/nsplit;
  add1s[0]   = &add1[add_srt[0]];  add1s[1]   = &add1[add_srt[1]];
  nl[0]      = num_links/nsplit;   nl[1]      = num_links-nl[0];
  add_end[0] = nl[0];              add_end[1] = num_links;

  depth = (int) (log(parent)/log(2));

#if defined(_OPENMP)
  /* Allow for nested parallelism */
  if ( omp_in_parallel() && depth<par_depth ) 
    {
      omp_set_nested(1);            
      if ( omp_get_nested() == 0 )
	cdoWarning("openMP implementation seems to not support nested parallelism.\n"
	       "Maximum of CPUs used is 2 instead of %i.\n", omp_get_num_threads());
    }                                    
#endif

#if defined(_OPENMP)
#pragma omp parallel for if(depth<par_depth) \
  private(i,who_am_i,my_depth) \
  num_threads(2)
#endif
  for ( i=0; i < nsplit; i++ )
    {
      who_am_i = nsplit*parent+i;
      my_depth = (int) (log(parent)/log(2))+1;

#if defined(_OPENMP)
      /*      if ( 0 )
      	cdoPrint("I am %i (parent %i), my_depth is: %i thread_num %i (%i) \n",
	who_am_i,parent,my_depth,omp_get_thread_num()+1,omp_get_num_threads());
      */
#endif
            
      sort_iter_single(nl[i], add1s[i], who_am_i);

    }

  /* ********************************* */
  /*              TO DO                */
  /* THIS BIT NEEDS TO BE PARALLELIZED */
  /* ********************************* */
  /* Idea I: one CPU merges top-down, the other one bottom-up */
  /*         ideally they should meet in the middle :)        */
  {
    //    uint64_t end,start;
    //    start = mach_absolute_time();                         /* ********************** */
    merge_lists(nl,add1s[0],add1s[1], idx);                     /* MERGE THE SEGMENTS     */
    //    end = mach_absolute_time();                           /* ********************** */
    //    merge_time += end-start;
  }

  tmp = (double*) malloc(num_links*sizeof(double));

#if defined(_OPENMP)
#pragma omp parallel for if ( depth < par_depth /* && num_links > 4096*/ ) \
    private(i) num_threads(2) schedule(static,1024)
#endif
  for ( i=0; i < num_links; i++ )
    tmp[i] = add1[idx[i]];
  
  memcpy(add1,tmp,num_links*sizeof(double));
  
  free(tmp);
  free(idx);

  tmp=NULL;

  return;
}

static
void sort_add(long num_links, double *restrict add1)
{
  /*
    This routine sorts address and weight arrays based on the
    destination address with the source address as a secondary
    sorting criterion. The method is a standard heap sort.

    Input and Output arrays:
    
       long num_links; ! num of links for this mapping
       double *add1    ! destination address array [num_links]
  */

  /* Local variables */

  double add1_tmp;         /* temp for addresses during swap     */
  long lvl, final_lvl;     /* level indexes for heap sort levels */
  long chk_lvl1, chk_lvl2, max_lvl;
  long i;

  if ( num_links <= 1 ) return;

  /*
    start at the lowest level (N/2) of the tree and shift lower 
    values to the bottom of the tree, promoting the larger numbers
  */
  for ( lvl = num_links/2-1; lvl >= 0; lvl-- )
    {
      final_lvl = lvl;
      add1_tmp = add1[lvl];

      /* Loop until proper level is found for this link, or reach bottom */

      for ( i = 0; i < num_links; i++ ) /* while ( TRUE ) */
	{
	  /* Find the largest of the two daughters */

          chk_lvl1 = 2*final_lvl+1;
          chk_lvl2 = 2*final_lvl+2;
          if ( chk_lvl1 == num_links-1 ) chk_lvl2 = chk_lvl1;

          if (add1[chk_lvl1] >  add1[chk_lvl2])
            max_lvl = chk_lvl1;
          else 
            max_lvl = chk_lvl2;

          /*
	    If the parent is greater than both daughters,
	    the correct level has been found
	  */
          if (add1_tmp >  add1[max_lvl])
	    {
	      add1[final_lvl] = add1_tmp;
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
	      final_lvl = max_lvl;
	      if ( 2*final_lvl+1 >= num_links )
		{
		  add1[final_lvl] = add1_tmp;
		  break;
		}
	    }
	}

      if ( i == num_links ) {
	cdoAbort("Internal problem; link 1 not found!\n");
	exit(1);
      }
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

      /* As above this loop sifts the tmp values down until proper level is reached */

      final_lvl = 0;

      for ( i = 0; i < num_links; i++ ) /* while ( TRUE ) */
	{
	  /* Find the largest of the two daughters */

          chk_lvl1 = 2*final_lvl+1;
          chk_lvl2 = 2*final_lvl+2;
          if ( chk_lvl2 >= lvl ) chk_lvl2 = chk_lvl1;

          if (add1[chk_lvl1] >  add1[chk_lvl2]) 
            max_lvl = chk_lvl1;
          else 
            max_lvl = chk_lvl2;

          /*
	    If the parent is greater than both daughters,
	    the correct level has been found
	  */
          if (add1_tmp >  add1[max_lvl]) 
	    {
	      add1[final_lvl] = add1_tmp;
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

	      final_lvl = max_lvl;
	      if ( 2*final_lvl+1 >= lvl )
		{
		  add1[final_lvl] = add1_tmp;
		  break;
		}
	    }
	}

      if ( i == num_links ) {
	cdoAbort("Internal problem; link 2 not found!\n");
      }
    }

  /* Swap the last two entries */
  add1_tmp = add1[1];
  add1[1]  = add1[0];
  add1[0]  = add1_tmp;

} /* sort_add */


void sort_iter_single(long num_links, double *restrict add1, int parent)
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

  static int MERGE_SORT_LIMIT_SIZE = 16384; 
  static int first_sort_iter_call = 1;
  static double merge_time;

  int par_depth = 1;
              
  if ( first_sort_iter_call )
    {
      first_sort_iter_call = 0; 
      par_depth = (int)(log(parent)/log(2));
      MERGE_SORT_LIMIT_SIZE = num_links/parent/2 ;
      if ( MERGE_SORT_LIMIT_SIZE < 4096 ) MERGE_SORT_LIMIT_SIZE = 4096;
      parent = 1;
    }
 
  if ( num_links > MERGE_SORT_LIMIT_SIZE )
    sort_par(num_links, add1, parent, par_depth);
  else
    sort_add(num_links,add1);

  if ( parent == 1 ) {
    first_sort_iter_call = 1;
    //    mach_timebase_info_data_t info = {0,0};
    //    mach_timebase_info(&info);
    //    merge_time = merge_time * (info.numer / info.denom)/1000./num_links;
    //    fprintf(stderr,"%12.8g ",merge_time);
    
    merge_time = 0;
  }

  return;
}
