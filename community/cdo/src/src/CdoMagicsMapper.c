#if defined(HAVE_CONFIG_H)
#  include "config.h" /* HAVE_LIBMAGICS */
#endif

#include "CdoMagicsMapper.h"
#include "magics_template_parser.h"

#define PARAM_COUNT  sizeof( mapper ) / sizeof ( CdoMagicsMapper )

/* extern int SetMagicsParameterValue( char *param_name, char *param_type, char *param_value ) */

int Set_magics_param_CCOLS( char *user_name, char *param_value );
int Reset_magics_param_CCOLS( char *user_name );

int Set_magics_param_CLEVS( char *user_name, char *param_value );
int Reset_magics_param_CLEVS( char *user_name );

int Set_magics_param_CTABLE( char *user_name, char *param_value );
int Reset_magics_param_CTABLE( char *user_name );

/* Define an array of Mapper structures to sort. */

typedef struct 

{
	char *cdo_name;
	char *magics_name;
	char *magics_type;
	int  (*Set_magics_param)( char *user_name, char *param_value );  /* Function to Update the Corresponding Magics parameters */
	int  (*Reset_magics_param)(char *user_name ); /* Function to Reset the Corresponding Magics parameters  */

} CdoMagicsMapper;


CdoMagicsMapper mapper[] =

{
	{ 
		"clevs", 
		"contour_level_list", 
		"floatarray",
		&Set_magics_param_CLEVS,
		&Reset_magics_param_CLEVS 
	},
	{ 
		"ccols",
		"contour_param_2",
		"intarray",
		&Set_magics_param_CCOLS,
		&Reset_magics_param_CCOLS
	},
	{ 
		"color_table",
                "contour_param_3",
                "intarray",
                &Set_magics_param_CTABLE,
                &Reset_magics_param_CTABLE
        }
};

int Compare( CdoMagicsMapper  *Parameter_one , CdoMagicsMapper *Parameter_two );

void PrintResult ( const CdoMagicsMapper *c );


int Set_magics_param_CCOLS( char *user_name, char *param_value )

{
	if( user_name == NULL )
		return 1;
	printf("Setting the CCOLS magics params \n"); 
        
        SetMagicsParameterValue( "contour_shade_colour_method", "string", "list" );
        SetMagicsParameterValue( "contour_shade_colour_list","stringarray", param_value );
#if 0
#endif
	return 0;
}


int Reset_magics_param_CCOLS( char *user_name )

{
	printf("Re-Setting the CCOLS magics params \n"); 
	return 0;
}


int Set_magics_param_CLEVS( char *user_name, char *param_value )

{
	if( user_name == NULL )
		return 1;

        SetMagicsParameterValue( "contour_level_selection_type","string", "level_list" );
        SetMagicsParameterValue( "contour_level_list","floatarray", param_value );

	return 0;
}


int Reset_magics_param_CLEVS( char *user_name )

{
        SetMagicsParameterValue( "contour_level_selection_type","string", "count" );
	printf("Re-Setting the CLEVS magics params \n"); 
	return 0;
}

int Set_magics_param_CTABLE( char *user_name, char *param_value )

{
	if( user_name == NULL )
		return 1;
	printf("Setting the CTABLE magics params \n"); 
#if 0
#if defined(HAVE_LIBMAGICS)
        SetMagicsParameterValue( "contour_level_list", "floatarray", param_value );
#endif
#endif
	return 0;
}


int Reset_magics_param_CTABLE( char *user_name )

{
	printf("Re-Setting the CTABLE magics params \n"); 
	return 0;
}


/* This is the comparison function used for sorting and searching. */
     
int Compare( CdoMagicsMapper *p1, CdoMagicsMapper *p2 )

{
	return strcmp ( p1->cdo_name, p2->cdo_name );
}
     
     
/* Print information about a critter. */
     
void PrintResult ( const CdoMagicsMapper *c )

{
	printf ( "CDO Name:%s\t MAGICS Name:%s\t MAGICS Type:%s\n", c->cdo_name, c->magics_name, c->magics_type );
}
     
     
/* Do the lookup into the sorted array. */

/* int GetMagicsParameterInfo( const char *user_name, char **magics_name, char **magics_type ) */

int GetMagicsParameterInfo( const char *user_name, char *param_value )

{
       static int once = 1;
       int ret_flag = 0;
       CdoMagicsMapper target, *result;
       target.cdo_name = (char *) user_name; 

       if( once )
       {
       		qsort ( mapper, PARAM_COUNT, sizeof ( CdoMagicsMapper ), ( void * )Compare );
		once = 0;
       }

       result = bsearch ( &target, mapper, PARAM_COUNT, sizeof ( CdoMagicsMapper ),
                          ( void * )Compare );
       if ( result )
       {
		result->Set_magics_param( result->cdo_name, param_value );

		/*
			*magics_name = result->magics_name;
			*magics_type = result->magics_type;
		*/
       }
       else
       {
	 /* Call the Reset functions of all the features to Reset the magics params to default in the calling function */
         ret_flag = 1;
       }
       return ret_flag;
}
