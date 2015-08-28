#if defined(HAVE_CONFIG_H)
#  include "config.h" /* HAVE_LIBMAGICS */
#endif

#include<limits.h>  /* TEMPORARY FIX, UNTIL NEXT MAGICS LIBRARY RELEASE */ 

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "pstream.h"

#if defined(HAVE_LIBMAGICS)
#include "magics_api.h"
#endif


#if defined(HAVE_LIBXML2)

#include<libxml/parser.h>
#include<libxml/tree.h>
#include "template_parser.h"
#include "magics_template_parser.h"
#include "results_template_parser.h"
#include <ctype.h>


extern xmlNode  *magics_node;

#endif

#define DBG 0


char *line_colours[] = {     "red", "green", "blue", "yellow", "cyan", "magenta",
			     "avocado","beige", "brick", "brown", "burgundy",
			     "charcoal", "chestnut", "coral", "cream", 
			     "evergreen", "gold", 
			     "khaki", "kellygreen", "lavender",
			     "mustard", "navy", "ochre", "olive",
			     "peach", "pink", "rose", "rust", "sky",
			     "tan", "tangerine","turquoise",
			     "violet", "reddishpurple",
			     "purplered", "purplishred",
			     "orangishred", "redorange", "reddishorange",
			     "orange", "yellowishorange",
			     "orangeyellow", "orangishyellow", 
			     "greenishyellow", "yellowgreen",
			     "yellowishgreen", "bluishgreen",
			     "bluegreen", "greenishblue",
			     "purplishblue", "bluepurple",
			     "bluishpurple", "purple",
			};

char  *graph_params[] = {"ymin","ymax","sigma","stat","obsv","device"};

int graph_param_count = sizeof(graph_params)/sizeof(char*);
int num_colours = sizeof( line_colours )/sizeof( char* );

void VerifyGraphParameters( int num_param, char **param_names );
int compareDateOrTimeStr( char *datetimestr1, char *datetimestr2, char *sep_char );

extern int checkdevice();
extern int IsNumeric();
extern void StrToUpperCase();
extern int StringSplitWithSeperator();

extern char *DEVICE;
extern char *DEVICE_TABLE;
extern int DEVICE_COUNT;


static
void maggraph(const char *plotfile, const char *varname,const char *varunits, long nfiles, long *nts, int **vdate, int **vtime, double **datatab, int nparam, char **params)
{
  
  char *lines[1];
  char *temp_str;
  char **split_str = NULL;
  char *sep_char = "=";
  char **date_time_str[nfiles];
  char min_date_time_str[1024], max_date_time_str[1024];
  int  min_index, max_index;
  char vdatestr[32], vtimestr[32], legend_text_data[256];
  char vdatestr1[32], vtimestr1[32];
  char vdatestr2[32], vtimestr2[32];
  int num_sigma = 2;
  int stat = FALSE, obsv = FALSE;
  int split_str_count;
  int file_begin = 0;
  int count ;
  int num_years = 0, num_months = 0, num_days = 0;
  int ret;
  long tsID, fileID, i, ntime_steps;
  double *date_time;
  double min_val = 1.0e+200, max_val = -1.0e+200;
  double *mean_val, *std_dev_val;
  double *spread_min, *spread_max;
  double y_min_val = 1.0e+200, y_max_val = -1.0e+200;
  
  if( DBG )
    {
      fprintf(stderr, "Num params %d\n", nparam);
  
      for( i = 0; i< nparam; i++ )
	fprintf(stderr, "Param %s\n", params[i]);
    }
  
  for( i = 0; i < nparam; ++i )
    {
      split_str_count = 0;
      sep_char = "=";
      split_str_count = StringSplitWithSeperator( params[i], sep_char, &split_str );
      
      if( !strcmp( split_str[0],"obsv" ) ) 
	{  
	  temp_str = strdup( split_str[1] );    
	  StrToUpperCase( temp_str );
	  if( !strcmp( temp_str, "TRUE" ) )
	    {
	      obsv = TRUE;
	      file_begin = 1;
	      if( DBG )
		fprintf( stderr,"OBSV TRUE\n" );
	    }
	}
	
      if( !strcmp( split_str[0],"stat" ) ) 
	{  
	  temp_str = strdup( split_str[1] );    
	  StrToUpperCase( temp_str );
	  
	  if( !strcmp( temp_str, "TRUE" ) )
	    {
	      stat = TRUE;
	      if( DBG )
		fprintf(stderr,"STAT TRUE\n");
	    }
	}
	
      if( !strcmp( split_str[0],"ymin" ) )
	{
	  y_min_val = atof( split_str[1] );
	  if( DBG )
	    fprintf(stderr,"Y min Val %g\n",y_min_val);
	}
	
      if( !strcmp( split_str[0],"ymax" ) )
	{
	  y_max_val = atof( split_str[1] );
	  if( DBG )
	    fprintf(stderr,"Y max Val %g\n",y_max_val);
	}
	
      if( !strcmp( split_str[0],"sigma" ) )
	{
	  num_sigma = atof( split_str[1] );
	  if( DBG )
	    fprintf(stderr,"SIGMA %d\n",num_sigma);
	}   
	
      if( !strcmp( split_str[0],"device" ) ) 
	{  
	  temp_str = strdup( split_str[1] );    
	  StrToUpperCase( temp_str );
	  DEVICE = temp_str;
	  if( DBG )
	    fprintf( stderr,"DEVICE %s\n",DEVICE );
	  
	  mag_setc ("output_format", DEVICE );
	}
	
      free( split_str );
    }

  if ( DBG )
    {
      ntime_steps = nts[0];
      fprintf(stderr," %ld %ld\n", nfiles, ntime_steps );
      fprintf(stderr,"STAT  %d\n", stat );
    }
    
  if ( stat == TRUE )
    {
        ntime_steps = nts[0];
	/* First date & time of first file */
	date2str( vdate[0][0],          vdatestr1, sizeof( vdatestr ) );
	date2str( vdate[0][ nts[0]-1 ], vdatestr2, sizeof( vdatestr ) );	
	
	/* Last date & time of first file */
	time2str( vtime[0][0],          vtimestr1, sizeof( vtimestr ) );
	time2str( vtime[0][ nts[0]-1 ], vtimestr2, sizeof( vtimestr ) );
	
	
	for ( fileID = 1; fileID < nfiles; fileID++ )
	  {
	    if ( nts[ fileID ] != ntime_steps )
	      {
		cdoWarning("  Unequal number of time steps! Statistics disabled.");
		stat = FALSE;
		break;
	      }
	    
	    /* First date & time of the present file */
	    date2str( vdate[ fileID ][0], vdatestr, sizeof( vdatestr ) );
	    sep_char = "-";
	    ret = compareDateOrTimeStr( vdatestr, vdatestr1, sep_char );
	    if( ret )
	      {
		cdoWarning("  Incosistent start date! Statistics disabled.");
		stat = FALSE;
		break;
	      }
	      
	    /* First time of the present file */	
	    time2str( vtime[ fileID ][0], vtimestr, sizeof( vtimestr ) );  
	    sep_char = ":";
	    ret = compareDateOrTimeStr( vtimestr, vtimestr1, sep_char );
	    if( ret )
	      {
		cdoWarning("  Incosistent start time! Statistics disabled.");
		stat = FALSE;
		break;
	      }
	      
	     /* Last date of the present file */
	    date2str( vdate[ fileID ][ nts[ fileID ]-1 ], vdatestr, sizeof( vdatestr ) );
	    sep_char = "-";
	    ret = compareDateOrTimeStr( vdatestr, vdatestr2, sep_char );
	    if( ret )
	      {
		cdoWarning("  Incosistent end date! Statistics disabled.");
		stat = FALSE;
		break;
	      }
	      
	    /* Last time of the present file */	
	    time2str( vtime[ fileID ][ nts[ fileID ]-1 ], vtimestr, sizeof( vtimestr ) );  
	    sep_char = ":";
	    ret = compareDateOrTimeStr( vtimestr, vtimestr2, sep_char );
	    if( ret )
	      {
		cdoWarning("  Incosistent end time! Statistics disabled.");
		stat = FALSE;
		break;
	      } 
	  }
    }
    
  if ( DBG )
    {
      fprintf(stderr,"STAT  %d\n", stat );
    }
    
  
  if ( stat == TRUE )
    {
	/* if all files are of same number of steps, only one date_time_str array is being used */
        date_time_str[0] = (char **) malloc( ntime_steps*sizeof(char *)); 
	
	date_time   = (double*) malloc( ntime_steps*sizeof(double));
	mean_val    = (double*) malloc( ntime_steps*sizeof(double));
	std_dev_val = (double*) malloc( ntime_steps*sizeof(double));
	spread_min  = (double*) malloc( ntime_steps*sizeof(double));
	spread_max  = (double*) malloc( ntime_steps*sizeof(double));
	
	for ( tsID = 0; tsID < ntime_steps; ++tsID )
	  {
	    date_time[tsID] = tsID+1;
	    date2str(vdate[0][tsID], vdatestr, sizeof(vdatestr));
	    time2str(vtime[0][tsID], vtimestr, sizeof(vtimestr));
	    date_time_str[0][tsID] = (char*) malloc(256);
	    sprintf(date_time_str[0][tsID], "%s %s", vdatestr, vtimestr);
	    mean_val[tsID] = 0.;
	    std_dev_val[tsID] = 0.;

	    if( DBG )
	      {
		fprintf(stderr,"%ld: %s\n", tsID, date_time_str[0][tsID]);
		fprintf(stderr,"%6d %6d", vdate[0][tsID], vtime[0][tsID]);
	      }
	
	    for ( fileID = 0; fileID < nfiles; ++fileID )
	      {
		if( DBG )
		  fprintf( stderr,"%ld\n", fileID );
		
		if( datatab[fileID][tsID] < min_val )
		  min_val = datatab[ fileID ][ tsID ];	
		if( datatab[fileID][tsID] > max_val )
		  max_val = datatab[ fileID ][ tsID ];	
	  
		mean_val[tsID] += datatab[fileID][tsID];
		std_dev_val[tsID] = 0.;
		spread_min[tsID] = 0.;
		spread_max[tsID] = 0.;

		if( DBG )
		  {
		    fprintf(stderr," %6g", datatab[fileID][tsID]);
		    fprintf(stderr,"\n");
		  }
	      }
	  }
	  
	for ( tsID = 0; tsID < ntime_steps; ++tsID )
	  {
	    mean_val[tsID] /= ( double )nfiles;
	    spread_min[tsID] = mean_val[tsID];
	    spread_max[tsID] = mean_val[tsID];

	    for ( fileID = 0; fileID < nfiles; ++fileID )
	      {
		std_dev_val[tsID] += ( datatab[fileID][tsID]-mean_val[tsID] ) * ( datatab[fileID][tsID]-mean_val[tsID] );
	      }
	    std_dev_val[tsID] /= ( double )nfiles;
	    std_dev_val[tsID] = pow( std_dev_val[tsID], 0.5 ); 
      
	    if( DBG )
	      fprintf(stderr," Mean : %g Std Dev: %g\n",mean_val[tsID],std_dev_val[tsID] ); 

	    spread_min[tsID] = mean_val[tsID] - num_sigma * std_dev_val[tsID];
	    spread_max[tsID] = mean_val[tsID] + num_sigma * std_dev_val[tsID];
	
	    if( DBG )
	      fprintf(stderr," Min : %g Max: %g\n",spread_min[tsID],spread_max[tsID] ); 
	  }

	for ( tsID = 0; tsID < ntime_steps; ++tsID )
	  {
	      if( spread_min[tsID] < min_val )
		  min_val = spread_min[ tsID ];	
	      if( spread_max[tsID] > max_val )
		  max_val = spread_max[ tsID ];	
	  }
	  
	if( DBG )
	  {
	    fprintf(stderr," %6g %6g\n", min_val, max_val );
	    fprintf(stderr," %s %s\n", date_time_str[0][0], date_time_str[0][ ntime_steps-1 ] );
	    fprintf(stderr,"\n");
	  }

	strcpy( min_date_time_str,date_time_str[0][0] );
	strcpy( max_date_time_str,date_time_str[0][ ntime_steps - 1 ] );
    }
    else
    {
      /* Find the min_date_time_str from the min's of nfiles
         Find the max_date_time_str from the max's of nfiles
         Construct the date_time_str array
      */
      
      if ( DBG )
	  fprintf(stderr,"STAT  %d\n", stat );
	
      for ( fileID = 0; fileID < nfiles; fileID++ )
	{
	  if ( DBG )
	    fprintf(stderr,"FILE  %ld\n", fileID );
	  date_time             = (double*) malloc( nts[fileID]*sizeof(double));
	  date_time_str[fileID] = (char **) malloc( nts[fileID]*sizeof(char *));
	  
	  for ( tsID = 0; tsID <  nts[fileID]; ++tsID )
	    {
	      date_time[tsID] = tsID+1;
	      date2str(vdate[fileID][tsID], vdatestr, sizeof(vdatestr));
	      time2str(vtime[fileID][tsID], vtimestr, sizeof(vtimestr));
	      
	      date_time_str[fileID][tsID] = (char*) malloc(256);
	      sprintf(date_time_str[fileID][tsID], "%s %s", vdatestr, vtimestr);
	      if ( DBG )
		fprintf( stderr,"%s %s %s\n", vdatestr, vtimestr, date_time_str[fileID][tsID] );
	      
	      if( datatab[fileID][tsID] < min_val )
	        min_val = datatab[ fileID ][ tsID ];	
	      if( datatab[fileID][tsID] > max_val )
	        max_val = datatab[ fileID ][ tsID ];	
	    }
	  free( date_time );
	  
	  if( fileID == 0 )
	    {
	      if ( DBG )
		fprintf( stderr,"\n %s %s\n", date_time_str[ fileID ][0], date_time_str[ fileID ][ nts[0]-1 ] );
	      min_index = 0;
	      max_index = 0;
	    }
	  else
	    {
	      if ( DBG )
		fprintf( stderr,"compareDateOrTimeStr %s\n", date_time_str[ fileID ][0] );
	      date2str( vdate[ min_index ][0], vdatestr1, sizeof( vdatestr ) );
	      date2str( vdate[ fileID ][0]   , vdatestr2, sizeof( vdatestr ) );	
	      sep_char = "-";
	      ret = compareDateOrTimeStr( vdatestr1, vdatestr2, sep_char );
	      if ( ret == -999 )
		cdoAbort("Error in input Date Time");
	      else if( ret == 1 )
		min_index = fileID;
	      else if( !ret )
	      {
		time2str( vtime[ min_index ][0], vtimestr1, sizeof( vtimestr ) );
		time2str( vtime[ fileID ][0]   , vtimestr2, sizeof( vtimestr ) );
		sep_char = ":";
		ret = compareDateOrTimeStr( vtimestr1, vtimestr2, sep_char );
					   
		if ( ret == -999 )
		  cdoAbort("Error in input Date Time");
		else if( ret == 1 )
		  min_index = fileID;			      
					   
	      }
	      if ( DBG )
		fprintf( stderr,"Min File ID %d\n",min_index);
	      
	      
	      if ( DBG )
		fprintf( stderr,"compareDateOrTimeStr  %s\n", date_time_str[ fileID ][ nts[ fileID ]-1 ] );
	      
	      date2str( vdate[ max_index ][ nts[ max_index ]-1 ], vdatestr1, sizeof( vdatestr ) );
	      date2str( vdate[ fileID ][ nts[ fileID ]-1 ]   , vdatestr2, sizeof( vdatestr ) );
	      sep_char = "-";
	      ret = compareDateOrTimeStr( vdatestr1, vdatestr2, sep_char );
					 
	      if ( ret == -999 )
		cdoAbort( "Error in input Date Time" );
	      else if( ret == -1 )
		max_index = fileID;
	      else if( !ret )
	      {
		time2str( vtime[ max_index ][ nts[ max_index ]-1 ], vtimestr1, sizeof( vtimestr ) );
		time2str( vtime[ fileID ][ nts[ fileID ]-1 ]   , vtimestr2, sizeof( vtimestr ) );
		sep_char = ":";
		ret = compareDateOrTimeStr( vtimestr1, vtimestr2, sep_char );
					   
		if ( ret == -999 )
		  cdoAbort("Error in input Date Time");
		else if( ret == -1 )
		  max_index = fileID;			      
	      }

            if( DBG )
	      fprintf( stderr,"Max File ID %d\n",max_index);

	    }
	}
	
	strcpy( min_date_time_str, date_time_str[ min_index ][0] );
	strcpy( max_date_time_str, date_time_str[ max_index ][ nts[ max_index ]-1 ] );
	if ( DBG )
	  fprintf( stderr,"%s %s\n",min_date_time_str, max_date_time_str );
    }
    
    if ( DBG )
      fprintf( stderr,"%s %s\n",min_date_time_str,max_date_time_str );

    split_str_count = 0;
    sep_char = "-";
    split_str_count = StringSplitWithSeperator( max_date_time_str, sep_char, &split_str );
    num_years  = atoi( split_str[0] );
    num_months = atoi( split_str[1] );
    num_days   = atoi( split_str[2] );
    free( split_str  );
    
    split_str_count = StringSplitWithSeperator( min_date_time_str, sep_char, &split_str );
    num_years -= atoi( split_str[0] );

    if( num_years <= 1 )
      {
	if( num_years == 1 )
	  num_months += ( 12 - atoi( split_str[1] ) );
	else
	  num_months -= ( atoi( split_str[1] ) );
	
	if( !num_months )
	  num_days -= atoi( split_str[2] );
	else if( num_months == 1 )
	  num_days += ( 31- atoi( split_str[2] ) );
      }
    free( split_str );
    
    if( DBG )
      fprintf(stderr," %d %d\n", num_years, num_months );

  /* 
	1. Loop over the Files
	2. Loop over the number of time steps 
	3. Set the attributes for the magics data and plot
  */  
   
#if defined(HAVE_LIBMAGICS)


  /* magics_template_parser( magics_node ); */

  mag_setc("output_name", plotfile);
  mag_setc("subpage_map_projection", "cartesian"); 
  mag_setr("subpage_y_length", 14.);
  mag_setr("subpage_y_position", 1.5);


  /* Horizontal Axis attributes */
  mag_setc("axis_orientation","horizontal");
  mag_setc("axis_grid", "on");
  mag_setc("axis_grid_colour", "grey");
  mag_seti("axis_grid_thickness", 1);
  mag_setc("axis_grid_line_style", "dot");
  mag_setc("axis_type", "date");
  
  if( num_years > 1 )
    mag_setc("axis_date_type", "years");
  else if( num_years <= 1 )
    {
      if( num_months > 1 )
	mag_setc("axis_date_type", "months");
      else
	{
	  if( num_months == 1 )
	    mag_setc("axis_date_type", "days");
	  else
	    {
	      if( num_days )
		mag_setc("axis_date_type", "days");
	      else
		mag_setc("axis_date_type", "hours");
	    }
	}
    }
  
  
  mag_setc("axis_date_min_value", min_date_time_str);
  mag_setc("axis_date_max_value", max_date_time_str);
  mag_setc("axis_title_text","Time");
  mag_setc("axis_title_orientation","horizontal");

  mag_seti("axis_tick_label_frequency", 2);
  mag_setr("axis_years_label_height", 0.4);

  mag_axis();

  /* Vertical Axis attributes */
  mag_setc("axis_orientation", "vertical");
  mag_setc("axis_grid", "on");
  mag_setc("axis_type", "regular");
  mag_setc("axis_grid_colour", "grey");
  mag_seti("axis_grid_thickness", 1);
  mag_setc("axis_grid_line_style", "dot");

  /*  To redefine the y- axis scale based on user input in .xml file */

  /* min & max values from the input data files */
  mag_setr("axis_min_value", min_val);
  mag_setr("axis_max_value", max_val);
  
  /* min & max values specified by the user in the command line args */
  if( y_min_val < 1.0e+200 )
    mag_setr("axis_min_value", y_min_val);
  
  if( y_max_val > -1.0e+200)
    mag_setr("axis_max_value", y_max_val);
  
  mag_setc("axis_title_text",varname);
  
  mag_setc("axis_title_orientation","vertical");

  mag_seti("axis_tick_label_frequency", 2);
  mag_setr("axis_tick_label_height", 0.5);

  mag_axis();
 

  /* Legend */
  mag_setc("legend", "on");
  mag_setc("legend_text_colour", "black");

  mag_setc("graph_symbol","off");
  mag_seti("graph_line_thickness", 8 );
  
  if( DBG )
    fprintf(stderr, "FILE BEGIN %d\n", file_begin );
  
  for ( i = file_begin; i < nfiles; ++i )
    {
      count = i; 
      if( obsv == TRUE )
	count = i -1;
      if( DBG )
	fprintf(stderr, "Current File %ld\n", i );
      /*sprintf(legend_text_data, "ens_%d", count + 1);*/
      sprintf(legend_text_data, "data_%d", count + 1);
      mag_setc("graph_line_colour", line_colours[ count%num_colours ]);
      mag_setc("legend_user_text", legend_text_data);
      if( stat == TRUE )
	mag_set1c("graph_curve_date_x_values",(const char**)date_time_str[0], ntime_steps);
      else
	mag_set1c("graph_curve_date_x_values",(const char**)date_time_str[i], nts[i]);

      /* TEMPORARY FIX, UNITL NEW MAGICS LIBRARY RELEASE *  begin**/
      mag_setr("graph_x_suppress_below",LLONG_MIN);
      mag_setr("graph_x_suppress_above",LLONG_MAX);
      /* TEMPORARY FIX, UNITL NEW MAGICS LIBRARY RELEASE *  end  **/

      mag_set1r("graph_curve_y_values", datatab[i], nts[i]);
      mag_graph ();
    }
      
  if( obsv == TRUE )
    {
      mag_setc("graph_line_colour", "black");
      sprintf(legend_text_data, "%s","Obsv" );
      mag_setc("legend_user_text", legend_text_data);
      mag_set1c("graph_curve_date_x_values",(const char**)date_time_str[0], nts[0]);

      /* TEMPORARY FIX, UNITL NEW MAGICS LIBRARY RELEASE *  begin**/
      mag_setr("graph_x_suppress_below",LLONG_MIN);
      mag_setr("graph_x_suppress_above",LLONG_MAX);
      /* TEMPORARY FIX, UNITL NEW MAGICS LIBRARY RELEASE *  end  **/

      mag_set1r("graph_curve_y_values", datatab[0], nts[0]);
      mag_setc("graph_line_style", "dot" );
      mag_seti("graph_line_thickness", 10 );
      mag_graph ();
    }
    
  if( DBG )
	fprintf(stderr, "NTIME STEPS %ld\n", ntime_steps ); 

  if( stat == TRUE )
    {
      if( DBG )
	fprintf(stderr, "NTIME STEPS %ld\n", ntime_steps );
      
      mag_seti("graph_line_thickness", 8 );
      mag_setc("graph_line_colour", "grey" );
      mag_setc("graph_line_style", "dash" );
      mag_set1c("graph_curve_date_x_values", (const char**)date_time_str[0], ntime_steps);

      /* TEMPORARY FIX, UNITL NEW MAGICS LIBRARY RELEASE *  begin**/
      mag_setr("graph_x_suppress_below",LLONG_MIN);
      mag_setr("graph_x_suppress_above",LLONG_MAX);
      /* TEMPORARY FIX, UNITL NEW MAGICS LIBRARY RELEASE *  end  **/

      mag_set1r("graph_curve_y_values",mean_val, ntime_steps);
      sprintf(legend_text_data, "Mean");
      mag_setc("legend_user_text", legend_text_data);
      mag_graph ();

      mag_reset("graph_type");
      mag_setc("graph_type", "area");
      mag_seti("graph_line_thickness", 1 );
      mag_setc("graph_shade_style", "dot");
      mag_setr("graph_shade_dot_size",1.);
      mag_set1c("graph_curve2_date_x_values", (const char**)date_time_str[0], ntime_steps);
      mag_set1r("graph_curve2_y_values",spread_max, ntime_steps);
      mag_set1c("graph_curve_date_x_values", (const char**)date_time_str[0], ntime_steps);
      mag_set1r("graph_curve_y_values",spread_min, ntime_steps);
      mag_setc("graph_shade_colour", "grey");
      sprintf(legend_text_data, "%dSigma", num_sigma);
      mag_setc("legend_user_text", legend_text_data);

      /* TEMPORARY FIX, UNITL NEW MAGICS LIBRARY RELEASE *  begin**/
      mag_setr("graph_x_suppress_below",LLONG_MIN);
      mag_setr("graph_x_suppress_above",LLONG_MAX);
      /* TEMPORARY FIX, UNITL NEW MAGICS LIBRARY RELEASE *  end  **/

      mag_graph ();
    }
  
  
  lines[0] = (char*) malloc(1024);
  /* To be obtained from Meta Data */
  /*sprintf( lines[0],"%s","ExpID : " );*/ 
  /*sprintf( lines[0],"%sxxxx  Variable : %s[%s]",lines[0], varname, varunits );*/
  // sprintf( lines[0],"Variable : %s[%s]",varname, varunits );
  // sprintf( lines[0],"%s  Date : %s --%s",lines[0], min_date_time_str, max_date_time_str );
  sprintf( lines[0],"Variable : %s[%s]  Date : %s --%s",varname, varunits, min_date_time_str, max_date_time_str );
  mag_set1c( "text_lines", (const char**)lines, 1 );
  
  mag_setc("text_html", "true");
  mag_setc("text_colour", "black");
  mag_setr("text_font_size", 0.6);
  mag_setc("text_mode", "positional");
  mag_setr("text_box_x_position", 1.5);
  mag_setr("text_box_y_position", 16.5);
  mag_setr("text_box_x_length", 20.);
  mag_setr("text_box_y_length", 2.5);
  mag_setc("text_border", "off");
  mag_setc("text_justification", "left");
  mag_text();

  if ( stat == TRUE )
    {
      free( date_time );
      free( mean_val );
      free( std_dev_val );
      free( spread_min );
      free( spread_max );
    }
  
  
  if( DBG )
    fprintf(stderr, "%s\n",lines[0]);

#endif

}

int compareDateOrTimeStr( char *datetimestr1, char *datetimestr2, char *sep_char )
{
  
  int    split_str_count1, split_str_count2;
  int	 i,flag[3]; /*  '3' since, three fields are expected in the input strings */
  char   **split_str1 = NULL;
  char   **split_str2 = NULL;
   
  if( DBG )
    fprintf(stderr,"Inside compareDateOrTimeStr %s %s\n",datetimestr1,datetimestr2);
  split_str_count1 = StringSplitWithSeperator( datetimestr1, sep_char, &split_str1 );
  
  if( split_str_count1 )
    split_str_count2 = StringSplitWithSeperator( datetimestr2, sep_char, &split_str2 );
  else
    {
      free( split_str1 );
      return -999;
    }
  
  if( split_str_count2 && !( split_str_count1 - split_str_count2 ) )
    {
      flag[0] = atoi( split_str1[0] ) - atoi( split_str2[0] ) ;
      flag[1] = atoi( split_str1[1] ) - atoi( split_str2[1] ) ;
      flag[2] = atoi( split_str1[2] ) - atoi( split_str2[2] ) ;
    }
  else
    {
      free( split_str1 );
      free( split_str2 );
      return -999;
    }
  
  free( split_str1 );
  free( split_str2 );
  
  for ( i = 0;i < 3 ; i++ )
    {
	if( flag[i] > 0 )
	    return 1;
	else if(  flag[i] < 0 )
	  return -1;
	else 
	  continue;
    }
    return 0;
}


#if defined(HAVE_LIBMAGICS)

static
void init_MAGICS( )

{
  setenv( "MAGPLUS_QUIET","1",1 ); /* To suppress magics messages */
  mag_open();

/* Some standard parameters affectng the magics environment, moved from the xml file  ** begin ** */
  mag_setc ("page_id_line","off");
/* Some standard parameters affectng the magics environment, moved from the xml file  ** end ** */

}


static
void quit_MAGICS( )

{

  mag_close ();
  if( DBG )
    fprintf( stdout,"Exiting From MAGICS\n" );

}

#endif

#define NINC_ALLOC 1024

void *Maggraph(void *argument)
{
  const char *ofilename;
  char varname[CDI_MAX_NAME], units[CDI_MAX_NAME];
  char **pnames = NULL;
  int varID, levelID;
  int gridID;
  int nrecs;
  int tsID;
  int streamID;
  int vlistID, vlistID0 = -1;
  int nmiss;
  int taxisID;
  int **vdate = NULL, **vtime = NULL;
  int fileID, nfiles;
  long *nts, nts_alloc;
  int nparam = 0;
  double **datatab = NULL;
  double val;
  int i;
  
  cdoInitialize(argument);

  nparam = operatorArgc();
  pnames = operatorArgv();
  
  if( nparam )
    VerifyGraphParameters(nparam,pnames);
  
  nfiles = cdoStreamCnt() - 1;
  ofilename = cdoStreamName(nfiles)->args;
  
  if( DBG )
    {
       fprintf( stderr," Num of files %d\n",nfiles );
       fprintf( stderr," files %s\n",ofilename );
    }
	
  datatab = (double **) malloc(nfiles*sizeof(double *));
  vdate   = (int **) malloc(nfiles*sizeof(int *));
  vtime   = (int **) malloc(nfiles*sizeof(int *));
  nts     = (long*) malloc(nfiles*sizeof(long));
  
  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      datatab[fileID] = NULL;
      vdate[fileID]   = NULL;
      vtime[fileID]   = NULL;
      nts[fileID]     = 0;
    }

  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      
     
      if( DBG )
        fprintf( stderr," file %d is %s\n", fileID, cdoStreamName(fileID)->args );
      streamID = streamOpenRead(cdoStreamName(fileID));

      vlistID = streamInqVlist(streamID);
      taxisID = vlistInqTaxis(vlistID);

      vlistInqVarUnits(vlistID, 0, units);
      if( DBG )
	fprintf(stderr," %s\n", units );
      if ( fileID == 0 )
	{
	  vlistInqVarName(vlistID, 0, varname);
	  
	  
	  gridID = vlistInqVarGrid(vlistID, 0);

	  if ( gridInqSize(gridID) != 1 ) cdoAbort("Variable has more than one grid point!");

	  vlistID0 = vlistDuplicate(vlistID);
	}
      else
	{
	  vlistCompare(vlistID0, vlistID, CMP_ALL);
	}

      tsID = 0;
      nts_alloc = 0;
      while ( (nrecs = streamInqTimestep(streamID, tsID)) )
	{
	  if ( nrecs != 1 ) cdoAbort("Input streams have more than one record!\n");
	  
	  if ( tsID == 0 )
	    {
	      nts_alloc += NINC_ALLOC;
	      datatab[ fileID ] = (double*) malloc( nts_alloc*sizeof(double));
	      vdate[ fileID ]   = (int*) malloc(  nts_alloc*sizeof(int));
	      vtime[ fileID ]   = (int*) malloc(  nts_alloc*sizeof(int));
	    }
		
	  nts[ fileID ]++;

	  if ( nts[ fileID ] > nts_alloc )
	    {
	      nts_alloc += NINC_ALLOC;
	      datatab[ fileID ] = (double*) realloc(datatab[fileID], nts_alloc*sizeof(double));
	      vdate[ fileID ]   = (int*) realloc(vdate[fileID], nts_alloc*sizeof(int));
	      vtime[ fileID ]   = (int*) realloc(vtime[fileID], nts_alloc*sizeof(int));
	    }
	  
	  streamInqRecord( streamID, &varID, &levelID );
	  streamReadRecord( streamID, &val, &nmiss );	
	  datatab[ fileID ][ tsID ] = val;
	  vdate[ fileID ][ tsID ] = taxisInqVdate(taxisID);
	  vtime[ fileID ][ tsID ] = taxisInqVtime(taxisID);

          if( DBG )
	    fprintf(stderr, "%f %f\n", datatab[ fileID ][ tsID ],val ); 
	  tsID++;
	}
      streamClose(streamID);
    }
  
#if defined(HAVE_LIBXML2)
  /* HARDCODED THE FILE NAME .. TO BE SENT AS COMMAND LINE ARGUMENT FOR THE MAGICS OPERATOR */
  /*
  init_XMLtemplate_parser( Filename );
  updatemagics_and_results_nodes( );
  */
#endif


#if defined(HAVE_LIBMAGICS)
  init_MAGICS( );
#endif

  cdoPrint(" Creating PLOT for %s", varname);
  if( DBG )
    {
      fprintf(stderr, "Num params %d\n", nparam);
  
      for( i = 0; i< nparam; i++ )
	fprintf(stderr, "Param %s\n", pnames[i]);
    }
  maggraph(ofilename, varname, units, nfiles, nts, vdate, vtime, datatab, nparam, pnames);

#if defined(HAVE_LIBXML2)
  /* quit_XMLtemplate_parser( ); */
#endif

#if defined(HAVE_LIBMAGICS)
  quit_MAGICS( );
#endif

  if ( vlistID0 != -1 ) vlistDestroy(vlistID0);

  for ( fileID = 0; fileID < nfiles; fileID++ )
    {
      if ( datatab[fileID] ) free(datatab[fileID]);
    }

  free(datatab);


  if ( vdate ) free(vdate);
  if ( vtime ) free(vtime);

  cdoFinish();

  return (0);
}


void VerifyGraphParameters( int num_param, char **param_names )

{
  int i, j;
  int  found = FALSE, syntax = TRUE, halt_flag = FALSE, split_str_count;
  char **split_str = NULL;
  char *sep_char = "=";
  char *temp_str;
  
  for ( i = 0; i < num_param; ++i )
    {
      split_str_count = 0;
      found = FALSE;
      syntax = TRUE;
      split_str_count = StringSplitWithSeperator( param_names[i], sep_char, &split_str );
      if( split_str_count > 1 ) 
	{
	  for ( j = 0; j < graph_param_count; ++j )
	    {
	      if( !strcmp( split_str[0], graph_params[j] ) )
		{
		  found = TRUE;
		  if( !strcmp( split_str[0],"obsv" ) ||  !strcmp( split_str[0],"stat" ) )
		    {  
		      if( IsNumeric( split_str[1] ) )
			syntax = FALSE;
		      else 
			{			
			  temp_str = strdup( split_str[1] );    
			  StrToUpperCase( temp_str );
			  if( strcmp( temp_str,"TRUE" ) && strcmp( temp_str,"FALSE" ) )
			    syntax = FALSE;			      
			}
		    }	 
		      
		  if( !strcmp( split_str[0],"ymin" ) ||  !strcmp( split_str[0],"ymax" ) || !strcmp( split_str[0],"sigma" )  )
		    {
		      if( !IsNumeric( split_str[1] ) )
			syntax = FALSE;       
		    }
		    
		    
      		  if( !strcmp( split_str[0],"device" ) )
		    {
		      if( IsNumeric( split_str[1] ) )
			syntax = FALSE;       
		      else 
			{
			  if( !strcmp( split_str[0],"device" ) )
			    {
			      if( DBG )
				fprintf( stderr,"Parameter value '%s'\n",split_str[1] );
			      if( checkdevice( split_str[1] ) )
				syntax = FALSE;

                              /* Graph not supported in google earth format */
                              if( !strcmp( split_str[1],"GIF_ANIMATION" ) || !strcmp( split_str[1],"gif_animation" ))
                                {
                                   syntax = FALSE;
	                           fprintf( stderr,"Animation not supported for Graph!\n");
                                   if( DBG )
                                     fprintf( stderr,"Parameter value '%s'\n",split_str[1] );
                                }
                              if( !strcmp( split_str[1],"KML" ) || !strcmp( split_str[1],"kml" ) )
                                {
                                   syntax = FALSE;
	                           fprintf( stderr," 'kml' format not supported for  Graph!\n");
                                   if( DBG )
                                     fprintf( stderr,"Parameter value '%s'\n",split_str[1] );
                                }
			    }
			}
		    }

/*		    
		  if( !strcmp( split_str[0],"xml" ) )
		    {
		      if( ( fp = fopen( split_str[1],"r") ) == NULL )
			{
			  fprintf( stderr,"Input XML File not found in specified path '%s'\n", split_str[1] );
			  halt_flag = TRUE;
			}
		      else
			{
#if defined(HAVE_LIBXML2)
			  // HARDCODED THE FILE NAME .. TO BE SENT AS COMMAND LINE ARGUMENT FOR THE MAGICS OPERATOR 
			  fclose(fp);
			  init_XMLtemplate_parser( split_str[1] );
			  updatemagics_and_results_nodes( );
#endif			
			}
		    }
*/
		}
	    }
	}
      else
	{
	  syntax = FALSE;
	}
	
      if( found == FALSE )
	{
	  halt_flag = TRUE;
	  fprintf( stderr,"Unknown parameter  '%s'!\n", param_names[i] );
	} 
      if( found == TRUE && syntax == FALSE )
	{
	  halt_flag = TRUE;
	  fprintf( stderr,"Invalid parameter specification  '%s'!\n", param_names[i] );
	}
      free( split_str );
    }
      
    if( halt_flag == TRUE )
    {
      exit(0);
    }
    
}
