#if defined(HAVE_CONFIG_H)
#  include "config.h" /* HAVE_LIBMAGICS */
#endif

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "grid.h"
#include "pstream.h"
#include "magics_api.h"

#include<libxml/parser.h>
#include<libxml/tree.h>
#include "template_parser.h"
#include "magics_template_parser.h"
#include "results_template_parser.h"

xmlDoc *param_doc = NULL;
xmlNode *root_node = NULL, *magics_node = NULL, *results_node = NULL;

#define DBG 0


int CONTOUR, SHADED, GRFILL;

char  *contour_params[] = {"min","max","count","interval","list","colour","thickness","style","RGB","device", "step_freq","file_split"};
int contour_param_count = sizeof(contour_params)/sizeof(char*);

char  *shaded_params[] = {"min","max","count","interval","list","colour_min","colour_max","colourtable","RGB","colour_triad","device","step_freq","file_split"};
int shaded_param_count = sizeof(shaded_params)/sizeof(char*);

char  *grfill_params[] = {"min","max","count","interval","list","colour_min","colour_max","colourtable","resolution","RGB","colour_triad","device","step_freq","file_split"};
int grfill_param_count = sizeof(grfill_params)/sizeof(char*);

char  *STD_COLOUR_TABLE[] = {"red", "green", "blue", "yellow", "cyan", "magenta", "black", "avocado",
			     "beige", "brick", "brown", "burgundy",
			     "charcoal", "chestnut", "coral", "cream", 
			     "evergreen", "gold", "grey", 
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
			     "bluishpurple", "purple", "white"
			    };


char **USR_COLOUR_TABLE = NULL;

int  STD_COLOUR_COUNT = sizeof( STD_COLOUR_TABLE )/sizeof( char* );
int  USR_COLOUR_COUNT =0;


char *STYLE_TABLE[] = { "SOLID","DASH","DOT","CHAIN_DASH","CHAIN_DOT"};
int STYLE_COUNT = sizeof( STYLE_TABLE )/ sizeof( char *);

char *DEVICE_TABLE[] = { "PS","EPS","PDF","PNG","GIF","GIF_ANIMATION","JPEG","SVG","KML"};
int DEVICE_COUNT = sizeof( DEVICE_TABLE )/ sizeof( char *);

int ANIM_FLAG = 0, STEP_FREQ = 0; /* '0' for static images like jpeg,ps, etc.. , '1' for animation formats */


int checkcolour( char *colour_in );
int ReadColourTable ( char *filepath );
int checkstyle( char *style_in );
int checkdevice( char *device_in );
void VerifyPlotParameters( int num_param, char **param_names, int opID );

extern int IsNumeric();
extern void StrToUpperCase();
extern void StrToLowerCase();
extern int StringSplitWithSeperator();
extern void StrReplaceChar( );

 /* Magics default values */
int COUNT = 10, isRGB = FALSE,   THICKNESS = 1, NUM_LEVELS = 0, FILE_SPLIT = FALSE;
double YMIN = 1.0e+200, YMAX = -1.0e+200, INTERVAL = 8.0, RESOLUTION = 10.0f, *LEV_LIST = NULL ;
char *COLOUR = NULL, *COLOUR_MIN = NULL, *COLOUR_MAX = NULL, *STYLE = NULL, *DEVICE = NULL, *COLOUR_TRIAD = NULL;


static
void magplot( const char *plotfile, int operatorID, const char *varname, const char *units, long nlon, long nlat, double *grid_center_lon, double *grid_center_lat, double *array,  int nparam, char **params, char *datetime )

{
  long i;
  double dlon = 0, dlat = 0;
  char plotfilename[4096];
  char *titlename;
  int j, split_str_count, split_str_count1;
  char *sep_char = "=";
  char **split_str = NULL, **split_str1 = NULL;
  char *temp_str = NULL;
  char  orig_char = ';', rep_char = ',';
  char tempname[256];
  
  
  if( DBG )
    {
      fprintf(stderr, "Num params %d\n", nparam);
  
      for( i = 0; i< nparam; i++ )
	fprintf(stderr, "Param %s\n", params[i]);
      fflush( stderr );
  
      for( i = 0; i < nparam; ++i )
        {
           split_str_count = 0;
           sep_char = "=";
           split_str_count = StringSplitWithSeperator( params[i], sep_char, &split_str );
	
           if( !strcmp( split_str[0],"min" ) )
	     fprintf(stderr," Min Val %g\n",YMIN );
	
           if( !strcmp( split_str[0],"max" ) )
	     fprintf(stderr,"Max Val %g\n",YMAX );
	
           if( !strcmp( split_str[0],"resolution" ) )
	     fprintf( stderr,"RESOLUTION %g\n",RESOLUTION );
	
           if( !strcmp( split_str[0],"colour" ) ) 
	     fprintf(stderr,"COLOUR %s\n",COLOUR );
	
           if( !strcmp( split_str[0],"colour_min" ) ) 
	     fprintf(stderr,"COLOUR %s\n",COLOUR_MIN );
	
           if( !strcmp( split_str[0],"colour_max" ) ) 
	     fprintf(stderr,"COLOUR %s\n",COLOUR_MAX );
	
           if( !strcmp( split_str[0],"interval" ) )
	     fprintf( stderr,"INTERVAL %f\n",INTERVAL );
	
           if( !strcmp( split_str[0],"count" ) )
	     fprintf( stderr,"COUNT %d\n",COUNT );
	
           if( !strcmp( split_str[0],"list" ) ) 
	     {
	        for( j = 0; j < split_str_count1; j++ )
	           fprintf( stderr,"LIST %f\n",LEV_LIST[j] ); 
	     }
	
           if( !strcmp( split_str[0],"thickness" ) )
	     fprintf( stderr,"THICKNESS %d\n",THICKNESS );
	
           if( !strcmp( split_str[0],"style" ) )
	     fprintf( stderr,"STYLE %s\n",STYLE );
	
           if( !strcmp( split_str[0],"device" ) )
	     fprintf( stderr,"DEVICE %s\n",DEVICE );

           if( !strcmp( split_str[0],"step_freq" ) )
	     fprintf( stderr,"STEP_FREQ %d\n",STEP_FREQ );

           free( split_str );
        }
    }

  if ( nlon > 1 )
    {
      for ( i = 1; i < nlon; ++i ) dlon += (grid_center_lon[i] - grid_center_lon[i-1]);
      dlon /= (nlon-1);
    }
  if ( nlat > 1 )
    {
      for ( i = 1; i < nlat; ++i ) dlat += (grid_center_lat[nlon*i] - grid_center_lat[nlon*(i-1)]);
      dlat /= (nlat-1);
    }

  sprintf( plotfilename, "%s [%s] %s", varname, units, datetime );
  titlename = strdup( plotfilename );
  sprintf( plotfilename, "%s_%s", plotfile, varname );

  mag_setc ("output_name",      plotfilename);
  mag_new( "page");

  /* Set the input data arrays to magics++ */
   
  mag_set2r("input_field", array, nlon, nlat);

  /*
  	mag_setc("input_field_organization", "REGULAR");
  	mag_set2r("input_field_latitudes", grid_center_lat, nlon, nlat);
  	mag_set2r("input_field_longitudes", grid_center_lon, nlon, nlat);
  */

  mag_setr("input_field_initial_latitude", grid_center_lat[0]);
  mag_setr("input_field_latitude_step", dlat);

  mag_setr("input_field_initial_longitude", grid_center_lon[0]);
  mag_setr("input_field_longitude_step", dlon);
 
  /* magics_template_parser( magics_node ); */
  /* results_template_parser(results_node, varname ); */


  /* set up the coastline attributes */
  /* mag_setc ("map_coastline_colour", "khaki"); */
  /* mag_setc ("map_grid_colour",      "grey");  */ 


  /* Parameters common to all operators */
  if( DEVICE )                                
    {
      mag_setc ("output_format", DEVICE );
    }

  mag_seti ("map_label_latitude_frequency",2);
  mag_seti ("map_label_longitude_frequency",2);
  /*mag_setr ("map_label_height",0.5);*/
  mag_setr ("map_label_height",0.4);


  /* define the contouring parameters */
  if ( operatorID == SHADED )
    {

      mag_setc ( "contour", "off" );
      mag_setc ( "contour_shade", "on" );
      mag_setc ( "contour_shade_method", "area_fill" );
      mag_setc ( "contour_label", "off" );
      
      if( YMIN < 1.0e+200  )
        {
	   mag_setr( "contour_shade_min_level", YMIN );
	   mag_setr( "contour_min_level", YMIN );
        }
      
      if( YMAX > -1.0e+200 )
        {
	   mag_setr( "contour_shade_max_level", YMAX );
	   mag_setr( "contour_max_level", YMAX );
        }
      
      if( COLOUR_MIN )
	mag_setc( "contour_shade_min_level_colour", COLOUR_MIN );
      
      if( COLOUR_MAX )
	mag_setc( "contour_shade_max_level_colour", COLOUR_MAX );

      if( INTERVAL != 8.0f )
	{
	  mag_setc( "contour_level_selection_type", "INTERVAL" );
	  mag_setr( "contour_interval", INTERVAL );
	}
	
      if( COUNT != 10 )
	{
	  mag_setc( "contour_level_selection_type", "COUNT" );
	  mag_seti( "contour_level_count", COUNT );
	}
	
      if( NUM_LEVELS  )
	{
	  mag_setc( "contour_level_selection_type", "LEVEL_LIST" );
	  mag_set1r( "contour_level_list", LEV_LIST, NUM_LEVELS );
	}
	
      if( USR_COLOUR_COUNT ) 
	{
	  mag_setc( "contour_shade_colour_method", "LIST" );
	  mag_set1c( "contour_shade_colour_list",( const char **)USR_COLOUR_TABLE, USR_COLOUR_COUNT ); 
	}
	
      if( COLOUR_TRIAD )                                
	{
	  mag_setc( "contour_shade_colour_direction", COLOUR_TRIAD );
	}
	
      
      /* Adjust Set The page slightly to fit the legend */
      mag_setr ( "subpage_x_length", 24. );
      mag_setr ( "subpage_y_length", 30. );

      /* Legend Settings */
      mag_setc ( "legend", "on" );
      mag_setc ( "legend_display_type", "continuous" );
      mag_setc ( "legend_entry_plot_direction", "column" );
      mag_setc ( "legend_box_mode", "positional" );
      mag_setr ( "legend_box_x_position", 26.5 );
      mag_setr ( "legend_box_y_position", 0.39 );
      mag_setr ( "legend_box_x_length", 2.0 );
      mag_setr ( "legend_box_y_length", 12.69 );

      if( DBG )
        {
           mag_enqc ( "output_name", &tempname );
           fprintf( stderr, " SHADED Done %s!\n",tempname );
           fprintf( stderr, " SHADED Done!\n" );
        }
    }
  else if ( operatorID == CONTOUR )
    {

      mag_setc ("contour",                  "on");
      mag_setc ("contour_shade",            "off");
      mag_setc ("contour_label",            "on");
      mag_setc ("contour_highlight",        "off");
      
    
      if( YMIN < 1.0e+200  )
	mag_setr( "contour_min_level", YMIN );

      if( YMAX > -1.0e+200 )
	mag_setr( "contour_max_level", YMAX );

      
      if( COLOUR )
	mag_setc( "contour_line_colour", COLOUR );
      
      
      if( INTERVAL != 8.0f )
	{
	  mag_setc( "contour_level_selection_type", "INTERVAL" );
	  mag_setr( "contour_interval", INTERVAL );
	}
	
      if( COUNT != 10 )
	{
	  mag_setc( "contour_level_selection_type", "COUNT" );
	  mag_seti( "contour_level_count", COUNT );
	}
	
      if( NUM_LEVELS  )
	{
	  mag_setc( "contour_level_selection_type", "LEVEL_LIST" );
	  mag_set1r( "contour_level_list", LEV_LIST, NUM_LEVELS );
	}
	
      if( THICKNESS != 1 )
	mag_seti( "contour_line_thickness", THICKNESS );
      
      if( STYLE )
      	  mag_setc( "contour_line_style", STYLE );
      
      if( DBG )
        fprintf( stderr, " CONTOUR Done!\n" );
    }
  else if ( operatorID == GRFILL )
    {

      mag_setc ( "contour", "off" );
      mag_setc ( "contour_shade", "on" );

      mag_setc ( "contour_shade_technique", "cell_shading" );

      mag_setc ( "contour_shade_method", "area_fill" );
      mag_setc ( "contour_label", "off" );
      
      if( YMIN < 1.0e+200  )
        {
	   mag_setr( "contour_shade_min_level", YMIN );
	   mag_setr( "contour_min_level", YMIN );
        }

      if( YMAX > -1.0e+200 )
        {
	   mag_setr( "contour_shade_max_level", YMAX );
	   mag_setr( "contour_max_level", YMAX );
        }

      /*
      if( YMIN < 1.0e+200  )
	mag_setr( "contour_shade_min_level", YMIN );
      
      if( YMAX > -1.0e+200 )
	mag_setr( "contour_shade_max_level", YMAX );
      */
      
      if( COLOUR_MIN )
	mag_setc( "contour_shade_min_level_colour", COLOUR_MIN );
      
      if( COLOUR_MAX )
	mag_setc( "contour_shade_max_level_colour", COLOUR_MAX );
      
      if( INTERVAL != 8.0f )
	{
	  mag_setc( "contour_level_selection_type", "INTERVAL" );
	  mag_setr( "contour_interval", INTERVAL );
	}
	
      if( COUNT != 10 )
	{
	  mag_setc( "contour_level_selection_type", "COUNT" );
	  mag_seti( "contour_level_count", COUNT );
	}
	
      if( NUM_LEVELS  )
	{
	  mag_setc( "contour_level_selection_type", "LEVEL_LIST" );
	  mag_set1r( "contour_level_list", LEV_LIST, NUM_LEVELS );
	}
	
      if( USR_COLOUR_COUNT )
	{
	  mag_setc( "contour_shade_colour_method", "LIST" );
	  mag_set1c( "contour_shade_colour_list",( const char ** ) USR_COLOUR_TABLE, USR_COLOUR_COUNT ); 
	}
	
      if( RESOLUTION != 10.0f)
	mag_setr( "contour_shade_cell_resolution", RESOLUTION );
      
      if( COLOUR_TRIAD )                                
	  mag_setc( "contour_shade_colour_direction", COLOUR_TRIAD );

      /* Adjust Set The page slightly to fit the legend */
      mag_setr ( "subpage_x_length", 24. );
      mag_setr ( "subpage_y_length", 30. );

      /* Legend Settings */
      mag_setc ( "legend", "on" );
      mag_setc ( "legend_display_type", "continuous" );
      mag_setc ( "legend_entry_plot_direction", "column" );
      mag_setc ( "legend_box_mode", "positional" );
      mag_setr ( "legend_box_x_position", 26.5 );
      mag_setr ( "legend_box_y_position", 0.39 );
      mag_setr ( "legend_box_x_length", 2.0 );
      mag_setr ( "legend_box_y_length", 12.69 );

      if( DBG )
        fprintf( stderr, " GrFILL Done!\n");
    }

  /* plot the title text and the coastlines */
  mag_cont ();
  mag_coast ();


  mag_set1c("text_lines", (const char **) &titlename, 1);
  mag_setc("text_colour", "black");

/*
  mag_setr("text_font_size", 0.6);
  mag_setc("text_mode", "positional");
  mag_setr("text_box_x_position", 1.5);
  mag_setr("text_box_y_position", 16.5);
  mag_setr("text_box_x_length", 20.);
  mag_setr("text_box_y_length", 2.5);
  mag_setc("text_border", "off");
*/

  mag_setc("text_justification", "left");
  mag_text();

  if( LEV_LIST ) 
    free( LEV_LIST );


}


static
void init_MAGICS( )
{
  setenv( "MAGPLUS_QUIET","1",1 ); /* To suppress magics messages */

  mag_open();
/* Some standard parameters affectng the magics environment, moved from the xml file  ** begin ** */
  mag_setc ("page_id_line","off");
  mag_setc(  "output_name_first_page_number", "off" );
  if( FILE_SPLIT == TRUE )
    mag_setc(  "output_ps_split" , "on" );
}

static
void quit_MAGICS( )
{

  mag_close ();
  if( DBG )
    fprintf( stderr,"Exiting From MAGICS\n" );

}

void *Magplot(void *argument)
{
  int operatorID;
  int varID, recID;
  int gridsize;
  int gridID;
  int nrecs;
  int levelID;
  int tsID;
  int streamID;
  int vlistID;
  int nmiss;
  int nlon, nlat;
  int nlev;
  int zaxisID, taxisID;
  int vdate, vtime;
  int nparam = 0;
  int i;
  char **pnames = NULL;
  char varname[CDI_MAX_NAME];
  double missval;
  double *array = NULL;
  double *grid_center_lat = NULL, *grid_center_lon = NULL;
  char units[CDI_MAX_NAME];
  char vdatestr[32], vtimestr[32], datetimestr[64];


  cdoInitialize(argument);
  
  
  nparam = operatorArgc();
  pnames = operatorArgv();
  
  CONTOUR = cdoOperatorAdd("contour", 0, 0, NULL);
  SHADED  = cdoOperatorAdd("shaded", 0, 0, NULL);
  GRFILL  = cdoOperatorAdd("grfill", 0, 0, NULL);

  operatorID = cdoOperatorID();
  
  if( nparam )
    {
      if( DBG )
	{
	  for( i = 0; i < nparam; i++ )
	    fprintf( stderr,"Param %d is %s!\n",i+1, pnames[i] );
	}
      
      VerifyPlotParameters( nparam, pnames, operatorID );
    }

  streamID = streamOpenRead(cdoStreamName(0));

  vlistID = streamInqVlist(streamID);
  taxisID = vlistInqTaxis(vlistID);

  varID = 0;
  gridID  = vlistInqVarGrid(vlistID, varID);
  zaxisID = vlistInqVarZaxis(vlistID, varID);
  missval = vlistInqVarMissval(vlistID, varID);

  if ( gridInqType(gridID) == GRID_GME          ) cdoAbort("GME grid unspported!");
  if ( gridInqType(gridID) == GRID_UNSTRUCTURED ) cdoAbort("Unstructured grid unspported!");

  if ( gridInqType(gridID) != GRID_CURVILINEAR )
    gridID = gridToCurvilinear(gridID, 1);

  gridsize = gridInqSize(gridID);
  nlon     = gridInqXsize(gridID);
  nlat     = gridInqYsize(gridID);
  nlev     = zaxisInqSize(zaxisID);

  array           = (double*) malloc(gridsize*sizeof(double));
  grid_center_lat = (double*) malloc(gridsize*sizeof(double));
  grid_center_lon = (double*) malloc(gridsize*sizeof(double));

  gridInqYvals(gridID, grid_center_lat);
  gridInqXvals(gridID, grid_center_lon);

  /* Convert lat/lon units if required */
  gridInqXunits(gridID, units);
  grid_to_degree(units, gridsize, grid_center_lon, "grid center lon");
  gridInqYunits(gridID, units);
  grid_to_degree(units, gridsize, grid_center_lat, "grid center lat");
					
  tsID = 0;

  /* HARDCODED THE FILE NAME .. TO BE SENT AS COMMAND LINE ARGUMENT FOR THE MAGICS OPERATOR */
  /*
     init_XMLtemplate_parser( Filename );
     updatemagics_and_results_nodes( );
  */


  init_MAGICS( );

  while ( (nrecs = streamInqTimestep(streamID, tsID)) )
    {
      if( ANIM_FLAG )
        {
      	  if( nrecs > 1 )
	    {
	      cdoWarning("File has more than one variable! Animation creation not possible!!! \n");
	      break;
            }
      	  if( tsID % STEP_FREQ )
	    {
                tsID++;
		continue;
            }
	}
      else 	
        {
          if( STEP_FREQ )
	    {
          	if( tsID % STEP_FREQ )
	    	  {
                     tsID++;
	             cdoWarning("NOT PLOTTING STEP %d!!!\n",tsID);
	             continue;
	    	  }
            }
         else 
            {
		if( tsID )
		  {
	   		cdoWarning("File variables have values at more than one time step! Images created for first time step!!!");
           		cdoWarning("To plot steps at a particular interval, set 'step_freq' to the frequency of the steps to be plotted!!!");
           		cdoWarning("To plot steps at random interval, set 'step_freq' to '1' and select the steps using the selection operators!!!");
       	   		break;
		  }
	    }
        }
      
      vdate = taxisInqVdate(taxisID);
      vtime = taxisInqVtime(taxisID);
	      
      date2str(vdate, vdatestr, sizeof(vdatestr));
      time2str(vtime, vtimestr, sizeof(vtimestr));
      sprintf( datetimestr, "%s %s", vdatestr, vtimestr );
      if( DBG )
        fprintf( stderr,"Date %s Time %s\n",vdatestr, vtimestr );

      for ( recID = 0; recID < nrecs; recID++ )
	{
	  streamInqRecord(streamID, &varID, &levelID);
	  streamReadRecord(streamID, array, &nmiss);
	  vlistInqVarName(vlistID, varID, varname);
	  vlistInqVarUnits(vlistID, varID, units);

	  if ( operatorID == SHADED || operatorID == CONTOUR || operatorID == GRFILL )
          {
                if( DBG )
                  {
                     if( operatorID == SHADED )
                       fprintf( stderr," Creating SHADED PLOT for %s\n",varname );
                     else if( operatorID == CONTOUR )
                       fprintf( stderr," Creating CONTOUR PLOT for %s\n",varname );
                     else if( operatorID == GRFILL )
                       fprintf( stderr," Creating GRFILL PLOT for %s\n",varname );
                  }

                if( DBG )
                  fprintf( stderr,"Plot %d\n",varID );
	  	magplot(cdoStreamName(1)->args, operatorID, varname, units, nlon, nlat, grid_center_lon, grid_center_lat, array, nparam, pnames, datetimestr );
          }
	  else
	  	fprintf(stderr,"operator not implemented\n");
	}

      if( DBG )
        fprintf( stderr,"TimeStep %d\n",tsID );

       
      tsID++;
      /*
      if( !STEP_FREQ  && tsID )
        {
	   cdoWarning("File variables have values at more than one time step! Images created for first time step!!!");
           cdoWarning("To plot steps at a particular interval, set 'step_freq' to the frequency of the steps to be plotted!!!");
           cdoWarning("To plot steps at random interval, set 'step_freq' to '1' and select the steps using the selection operators!!!");
       	   break;
	}
      else
        {
      	   tsID++;
           if( DBG )
             fprintf( stderr,"TimeStep %d\n",tsID );
	}
      */
    }

  if( ANIM_FLAG )
    {
      if( FILE_SPLIT == TRUE  ) 
        cdoWarning("File split parameter ignored!!!");
    }
  quit_MAGICS( );

  streamClose(streamID);

  if ( array  ) free(array);
  if ( grid_center_lon ) free(grid_center_lon);
  if ( grid_center_lat ) free(grid_center_lat);

/*   quit_XMLtemplate_parser( ); */

  cdoFinish();

  return (0);

}


void VerifyPlotParameters( int num_param, char **param_names, int opID )

{
  int i, j, k;
  int found = FALSE, syntax = TRUE, halt_flag = FALSE, /* file_found = TRUE, */ split_str_count;
  int param_count;
  char **params;
  char **split_str = NULL, **split_str1 = NULL;
  char *sep_char = "=";
  char *temp_str;
  char  orig_char = ';', rep_char = ',';
  FILE *fp;

/*  
  char  *contour_params[] = {"ymin","ymax","count","interval","list","colour","thickness","style"};
  char  *shaded_params[]  = {"ymin","ymax","count","interval","list","colour_min","colour_max","colortable","step_freq"};
  char  *grfill_params[]  = {"ymin","ymax","count","interval","list","colour_min","colour_max","colortable","resolution"};
*/


  for ( i = 0; i < num_param; ++i )
    {
      split_str_count = 0;
      found = FALSE;
      syntax = TRUE;
      split_str_count = StringSplitWithSeperator( param_names[i], sep_char, &split_str );
      
      if( DBG )
	fprintf( stderr, "Verifying params!\n");
      
      if( split_str_count > 1 ) 
	{
	  
	  if( opID == CONTOUR )
	    {
	      param_count = contour_param_count;
	      params = contour_params;
	    }
	  else if( opID == SHADED )
	    {
	      param_count = shaded_param_count;
	      params = shaded_params;
	    }
	  else if( opID == GRFILL )
	    {
	      param_count = grfill_param_count;
	      params = grfill_params;
	    }
	  
	  for ( j = 0; j < param_count; ++j )
	    {
	      if( !strcmp( split_str[0], params[j] ) )
		{
		  found = TRUE;
		  if( !strcmp( split_str[0],"colour" )     || !strcmp( split_str[0],"style" )       ||
		      !strcmp( split_str[0],"colour_min" ) || !strcmp( split_str[0],"colour_max" )  ||
		      !strcmp( split_str[0],"RGB" )        || !strcmp( split_str[0],"colour_triad" )||
		      !strcmp( split_str[0],"device")      || !strcmp( split_str[0],"file_split" )
		    )
		    {
		      if( IsNumeric( split_str[1] ) )
			syntax = FALSE;
		      else
			{
			  if( !strcmp( split_str[0],"RGB" ) || !strcmp( split_str[0],"file_split") )
			    {
			      temp_str = strdup( split_str[1] );    
			      StrToUpperCase( temp_str );
			      if( strcmp( temp_str,"TRUE" ) && strcmp( temp_str,"FALSE" ) )
				syntax = FALSE;      
			      else
				{
                                  if( !strcmp( split_str[0],"RGB" ) )
                                    {
				      if( !strcmp( temp_str,"TRUE" ) )
				        isRGB = TRUE;
				      else
				        isRGB = FALSE;
				    }
                                  else if( !strcmp( split_str[0],"file_split" ) )
                                    {
				      if( !strcmp( temp_str,"TRUE" ) )
				        FILE_SPLIT = TRUE;
				      else
				        FILE_SPLIT = FALSE;
				    }
				}
			    }
			  else if( !strcmp( split_str[0],"style" ) )
			    {
			      if( checkstyle( split_str[1] ) )
				syntax = FALSE;
			    }
			  else if( !strcmp( split_str[0],"colour" ) || !strcmp( split_str[0],"colour_min" ) || !strcmp( split_str[0],"colour_max" ) )
			    {
			      
			      if( checkcolour( split_str[1] ) )
				syntax = FALSE;
                              else
                                {
      				   if( !strcmp( split_str[0],"colour" ) ) 
	                             {  
	                               temp_str = strdup( split_str[1] );  
	                               if( !isRGB )
	                                 StrToLowerCase( temp_str );
	                               else
	                                 {
	                                    StrToUpperCase( temp_str );
	                                    StrReplaceChar( temp_str, orig_char, rep_char ); /* replace ';' in RGB format to ',' */
	                                 }
                                       COLOUR = temp_str;
	                               if( DBG )
	                                 fprintf(stderr,"COLOUR %s\n",COLOUR );
	                             }
                                   if( !strcmp( split_str[0],"colour_min" ) ) 
	                             {  
	                               temp_str = strdup( split_str[1] );    
                                       if( !isRGB )
	                                 StrToLowerCase( temp_str );
	                               else
	                                 {
	                                   StrToUpperCase( temp_str );
	                                   StrReplaceChar( temp_str, orig_char, rep_char ); /* replace ';' in RGB format to ',' */
	                                 }
	                               COLOUR_MIN = temp_str;
	                               if( DBG )
	                                 fprintf(stderr,"COLOUR %s\n",COLOUR_MIN );
	                             }
                                   if( !strcmp( split_str[0],"colour_max" ) ) 
	                             {
	                               temp_str = strdup( split_str[1] );    
                                       if( !isRGB )
	                                 StrToLowerCase( temp_str );
	                               else
	                                 {
	                                   StrToUpperCase( temp_str );
	                                   StrReplaceChar( temp_str, orig_char, rep_char ); /* replace ';' in RGB format to ',' */
	                                 }
	                               COLOUR_MAX = temp_str;
	                               if( DBG )
	                                 fprintf(stderr,"COLOUR %s\n",COLOUR_MAX );
	                             }
			        }
			    }
			  else if( !strcmp( split_str[0],"device" ) )
			    {
			      if( checkdevice( split_str[1] ) )
  				    syntax = FALSE;
			    }
			  else if( !strcmp( split_str[0],"colour_triad" ) )
			    {
			      temp_str = strdup( split_str[1] );    
			      StrToUpperCase( temp_str );
			      if( strcmp( temp_str,"CW" ) && strcmp( temp_str,"ACW" ) )
				syntax = FALSE;      
			      else
				{
				   if( DBG )
				     fprintf( stderr, "TRIAD check  %s!\n",temp_str);
				   if( !strcmp( temp_str,"CW" ) )
				     COLOUR_TRIAD = "clockwise";
				   else
				     COLOUR_TRIAD = "anti_clockwise";
				}
			    }
			}
		    }
		      
		  if( !strcmp( split_str[0],"min" )      ||  !strcmp( split_str[0],"max" )     ||
		      !strcmp( split_str[0],"count" )     ||  !strcmp( split_str[0],"interval" ) ||
		      !strcmp( split_str[0],"thickness" ) ||  !strcmp( split_str[0],"resolution" ) || 
		      !strcmp( split_str[0],"step_freq" )
                    )
		    {
		      if( !IsNumeric( split_str[1] ) )
			syntax = FALSE;
		      else
			{
                           if( !strcmp( split_str[0],"min" ) )
	                     {
	                        YMIN = atof( split_str[1] );
			     }
                           if( !strcmp( split_str[0],"max" ) )
	                     {
	  		        YMAX = atof( split_str[1] );
			     }
		           if( !strcmp( split_str[0],"count" ) )
			     {
	                        COUNT = atoi( split_str[1] );
			     }
                           if( !strcmp( split_str[0],"interval" ) )
	                     {
	  		        INTERVAL = atof( split_str[1] );
	                     }
		           if( !strcmp( split_str[0],"thickness" ) )
			     {
	  		        THICKNESS = atoi( split_str[1] );
			     }
                           if( !strcmp( split_str[0],"resolution" ) )
	                     {
	                        RESOLUTION = atoi( split_str[1] );
	                     }	
		           if( !strcmp( split_str[0],"step_freq" ) )
			     {
	  	                STEP_FREQ = atoi( split_str[1] );
			     }
	                }
		    }
		    
		  if( !strcmp( split_str[0],"colourtable" ) )
		    {
		      if( ( fp = fopen( split_str[1],"r") ) == NULL )
			{
			  fprintf( stderr,"Input Color Table File not found in specified path '%s'\n", split_str[1] );
			  halt_flag = TRUE;
			}
		      else
			{
			  ReadColourTable ( split_str[1] );
			}
		    }
		    
		  if( !strcmp( split_str[0],"list" ) )
		    {
		      sep_char = ";";
		      split_str_count = StringSplitWithSeperator( split_str[1], sep_char, &split_str1 );
		      if( !split_str_count )
			{
			  syntax = FALSE;
			}
		      else
		        {
			  for( k = 0; k < split_str_count; k++ )
			    {
			      if( !IsNumeric( split_str1[k] ) )
			        syntax = FALSE;
			    }
			  if( syntax == TRUE )
			    {
	                       NUM_LEVELS = split_str_count;
	                       LEV_LIST = (double*) malloc( sizeof( double ) * split_str_count );
			       for( k = 0; k < split_str_count; k++ )
		                 {
		                    LEV_LIST[k] = atof( split_str1[k] );
		                 }
	                       free( split_str1 );
	                    }
			  }
		      }
		  sep_char = "=";
		} /*** if( !strcmp( split_str[0], params[j] ) )  ***/
	    } /*** Loop over param count ***/
	} /*** ( split_str_count > 1 ) ***/
      else
	{
	  syntax = FALSE;
	}
	
      if( found == FALSE )
	{
	  halt_flag = TRUE;
	  fprintf( stderr,"Invalid parameter  '%s'\n", param_names[i] );
	} 
      if( found == TRUE && syntax == FALSE )
	{
	  halt_flag = TRUE;
	  fprintf( stderr,"Invalid parameter specification  '%s'\n", param_names[i] );
	}
	
      if( split_str ) 	  
	free( split_str );
    } /*** Loop over params ****/
      
    if( halt_flag == TRUE )
      {
        exit(0);
      }
    
}


int checkcolour( char *colour_in )

{

    int i, n, found = FALSE;
    int split_str_count;
    char *sep_char =",";
    char **split_str = NULL;
    float  rgb_values[3];
    char temp[256];
    char *ref;
   
    ref = colour_in;
    
    if( isRGB )
      {
	if( strchr( colour_in,';') == NULL || strstr( colour_in,"RGB(") == NULL )
	  {
	    cdoWarning( "Found 'RGB=true',Specify Colour in 'RGB(r;g;b)' ( where r,g,b in [0.0,1.0] ) format!" );
	    free( split_str );
	    return 1;
	  }
	  
	n = strlen( colour_in );
    
	if( DBG )
	  fprintf( stdout,"  count %d  original colour %s RGB %d\n", n, colour_in, isRGB  );
	
	for( i=0 ; i< n-1; i++ )
	  {
	    if( i > 3 )
	      { 
		temp[i-4] = *colour_in;
	      }
	    colour_in++; 
	  }
	  
	temp[i-4] = '\0';
	
	if( DBG )
	  fprintf( stdout,"  count %d  modified color %s \n", (int)strlen(temp), temp  );
	
	sep_char =";";
	split_str_count = StringSplitWithSeperator( temp, sep_char, &split_str );
    
	if(  split_str_count != 3 ) 
	  {
	    cdoWarning( " Colour specified in Improper format!" );
	    free( split_str );
	    return 1;
	  }
    
	rgb_values[0] = atof( split_str[0] );
	rgb_values[1] = atof( split_str[1] );
	rgb_values[2] = atof( split_str[2] );
    
	if( rgb_values[0] + rgb_values[1] + rgb_values[2] > 3.0f  || 
	    rgb_values[0] + rgb_values[1] + rgb_values[2] < 0.0f 	   )
	  {
	    cdoWarning( " RGB Colour specified with Improper values!" );
	    free( split_str );
	    return 1;
	  }
	  
	free( split_str );  
      }
    else
      {
	if( strchr( colour_in,';') != NULL || strstr( colour_in,"RGB(") != NULL )
	  {
	    cdoWarning( "Found Colour with 'RGB(r;g;b)' format, set parameter RGB='true' !" );
	    free( split_str );
	    return 1;
	  }
	  
	StrToLowerCase( colour_in );
	for( i = 0 ; i < STD_COLOUR_COUNT; i++ )
	  {
	    if( !strcmp( STD_COLOUR_TABLE[i], colour_in ) )
	      {
		found = TRUE;
		return 0;
	      }
	  }
	  cdoWarning( "Specified Colour not in Standard colour list, resetting to blue(default colour)!" );
	  return 1;
      }
      
    if( DBG )  
      cdoWarning( "Colour %s verified!",ref );  
    return 0;
}


int ReadColourTable ( char *filepath )

{
    
    FILE *fp;
    int  i, num_colors;
    char **temp_table = NULL;
    char  orig_char = ';', rep_char = ',';
    
    fp = fopen( filepath,"r" );
    
    if( !fp )
      {
	fprintf( stdout, "File Not available!" );
	return 1;
      }
    
    fscanf( fp, "%d", &num_colors );
    
    if( DBG )
      fprintf( stderr, "Num Colours %d\n", num_colors );
    
    if( !num_colors )
      {
	cdoWarning("No colours found in File, proceeding with Standard Colour table!\n");
	fclose(fp);
	return 1;
      }
    
    USR_COLOUR_COUNT = 0;
    USR_COLOUR_TABLE = ( char **) malloc( num_colors * sizeof( char* ));
    temp_table  = ( char **) malloc( num_colors * sizeof( char* ));
    
    for( i =0; i < num_colors; i++ )
      {
         temp_table[i] = ( char *) malloc(  256 * sizeof( char ));
	 fscanf( fp, "%s", temp_table[i] );
	 if( DBG )
	   fprintf( stdout, "%s\n", temp_table[i] );
      } 
    
    for( i = 0; i < num_colors; i++ )
      {
	  if( DBG )
	    fprintf( stdout, "%s \n", temp_table[i] );
	  
	  if( !checkcolour( temp_table[i] ) )
	    {
	      if( isRGB )
		StrReplaceChar( temp_table[i], orig_char, rep_char ); /* replace ';' in RGB format to ',' */

	      if( DBG )
		fprintf( stdout, "Before appending %s\n", temp_table[i] );
	      
	      USR_COLOUR_TABLE[ USR_COLOUR_COUNT ] = strdup( temp_table[i] );
	      
	      /* strcpy( USR_COLOUR_TABLE[ USR_COLOUR_COUNT ], temp_table[i] ); */
	      USR_COLOUR_COUNT++;
	      
	      if( DBG )
		fprintf( stdout, "After appending %s\n", temp_table[i] );
	    }
      }
    
    if( USR_COLOUR_COUNT < num_colors )
      {
	  cdoWarning( " Discarding improper format colours and continuing!\n" );
      }
      
    fclose(fp);   
    return 0;
}

int checkstyle( char *style_in )

{
    int i, found = FALSE;
    StrToUpperCase( style_in );
    for( i = 0 ; i < STYLE_COUNT; i++ )
      {
	if( DBG )
	  fprintf( stderr, "Input %s ref %s\n",style_in, STYLE_TABLE[i] );
	
	if( !strcmp( STYLE_TABLE[i], style_in ) )
	  {
	    found = TRUE;
	    STYLE = style_in;
	    return 0;
	  }
      }
      
    if( !found )
	 cdoWarning( " Style specified with Improper value!\n" );
    
    return 1; 
}


int checkdevice( char *device_in )

{
    int i, found = FALSE;
    StrToUpperCase( device_in );
    for( i = 0 ; i < DEVICE_COUNT; i++ )
      {
	if( DBG )
	  fprintf( stderr, "Input %s ref %s\n",device_in, DEVICE_TABLE[i] );
	
	if( !strcmp( DEVICE_TABLE[i], device_in ) )
	  {
	    found = TRUE;

	    DEVICE = device_in;
	    if( DBG )
	      fprintf( stderr,"DEVICE %s\n",DEVICE );

	    if( !strcmp( "GIF_ANIMATION" , device_in ) || !strcmp( "KML", device_in )  )
              {
	         ANIM_FLAG = 1;
		 STEP_FREQ = 1;
	      }
	    return 0;
	  }
      }
      
    if( !found )
	 cdoWarning( " Device specified with Improper value!\n" );
    
    return 1; 
}
