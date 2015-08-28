#if defined(HAVE_CONFIG_H)
#  include "config.h" /* HAVE_LIBMAGICS */
#endif

#include "magics_template_parser.h"
#include "StringUtilities.h"

#include "magics_api.h"

#define DBG 0 


extern xmlNode *magics_node;


/* Recursive function that sets the Magics parameters from the XML structure */

int magics_template_parser( xmlNode *a_node ) 

{
    int param_set_flag;
    xmlNode *cur_node = NULL;
    xmlChar    *param_name,*param_type,*param_value,*value;

    if( a_node == NULL )
        return 0;

#if 0
    fprintf( stdout,"Parsing the magics Node \n");
#endif

    if( !strcmp( a_node->name, "magics" ) )
    {
 	value = xmlGetProp( a_node, "version" );

        if( value )
        {
	    	if( DBG )
			printf( "Version %s \n", value ); 

		if( atof( value ) > 3.0f ) 
		{
			return 1;
		}
        }

    }


    for ( cur_node = a_node->children; cur_node; cur_node = cur_node->next )
    {
	param_name = NULL;
	param_type = NULL;
	param_value = NULL;

        if ( cur_node->type == XML_ELEMENT_NODE )
        {
		
	    if( DBG )
            	printf( "Node Name: %s \n", cur_node->name );

#if 0
            fprintf( stdout,"Node Name: %s \n", cur_node->name );
#endif

	    if( cur_node->properties == NULL )
	    {
		if( cur_node->children == NULL )
		{
			printf( "NO ATTRIBUTES!!!\n" );
		}
	    }
	    else
	    {
		
		param_name = xmlGetProp( cur_node,"parameter");
		param_type = xmlGetProp(cur_node,"type");
		param_value = xmlGetProp(cur_node,"value");
#if 0
    		printf( "\t\tAttr name: %s Type: %s Value: %s \n", param_name,param_type,param_value);
#endif
		
    		param_set_flag = SetMagicsParameterValue( param_name, param_type, param_value );
		
		if( param_set_flag )
			printf(" Error in Setting the Parameter %s\n",param_name );
	    }
        }
    }
    return 0;
}

int SetMagicsParameterValue( char *param_name, char *param_type, char *param_value )

{
	int i, ret_flag = 0;
	int split_str_count = 0;
	char **split_str = NULL;
	char *sep_char = ",";
	char *search_char = ";";
	double *float_param_list = NULL;
	int *int_param_list = NULL;

	if( param_name == NULL )
	  {
		ret_flag = 1;
		return ret_flag;
	  }

	if( param_value == NULL )
		ret_flag = 2;


	
   	/*   MAGICS++ ENV RELATED PARAMETERS   */
	if( !strcmp( param_type,"environvar" ) )
	  {
		if( !strcmp( param_name,"quiet_option" ) )
		  {
			if( !strcmp( param_value, "off" ) || !strcmp( param_value, "OFF" ) )
			  {
#if 0
				printf( "Quiet Option %s \n", param_value ); 
#endif
				if( !unsetenv( "MAGPLUS_QUIET" ) )
				  {
				      if( DBG )
					  fprintf( stderr, "Quiet Option %s is un-set successfully!!! \n", param_value ); 
				  }
				else
					fprintf( stderr, "Quiet Option %s COULDN'T be UNSET!!!\n", param_value ); 
			  }

			if( !strcmp( param_value, "on" ) || !strcmp( param_value, "ON" ) )
			  {
#if 0
				printf( "Quiet Option %s \n", param_value ); 
#endif
				if( !setenv( "MAGPLUS_QUIET","1",1 ) )
				  {
					if( DBG )
					  fprintf( stderr, "Quiet Option %s is set successfully!!! \n", param_value ); 
				  }
				else
					fprintf( stderr, "Quiet Option %s COULDN'T be SET!!!\n", param_value ); 
			  }
		  }
	  }

    	/*   MAGICS++ FLOAT TYPE PARAMETERS   */
	else if( !strcmp( param_type,"float" ) )
	  {
		mag_setr( param_name, atof( param_value ) );		
	  }

    	/*   MAGICS++ FLOAT ARRAY  TYPE    PARAMETERS   */
	else if( !strcmp( param_type,"floatarray" ) )
	  {

#if 0
	        fprintf(stderr, "param_name : %s\tparam_value: %s\n", param_name, param_value);
#endif
		if( strchr( param_value,';') )
		    sep_char = ";";
		split_str_count = StringSplitWithSeperator( param_value, sep_char, &split_str );
		if( split_str_count )
		  {
		    float_param_list = (double*) malloc(sizeof(double) * split_str_count );
			for( i = 0; i < split_str_count; i++ )
			  {
#if 0
			        fprintf(stderr, "%d %d %s\n", i, split_str_count, split_str[i]);
#endif
				float_param_list[i] = atof( split_str[i] );			
			  }
			mag_set1r( param_name, float_param_list, split_str_count );		
			free( float_param_list );
			free( split_str );
		  }

	  }

    	/*   MAGICS++ INT TYPE    PARAMETERS   */
	else if( !strcmp( param_type,"int" ) )
	  {
		mag_seti( param_name, atoi( param_value ) );		
	  }

    	/*   MAGICS++ INT ARRAY  TYPE    PARAMETERS   */
	else if( !strcmp( param_type,"intarray" ) )
	  {
	        if( strchr( param_value,';') )
		    sep_char = ";";
		split_str_count = StringSplitWithSeperator( param_value, sep_char, &split_str );
		if( split_str_count )
		  {
		    int_param_list = (double*) malloc(sizeof( int ) * split_str_count );
			for( i = 0; i < split_str_count; i++ )
			{
				int_param_list[i] = atoi( split_str[i] );			
			}
			mag_set1i( param_name, int_param_list, split_str_count );		
			free( int_param_list );
			free( split_str );
		  }
	  }

    	/*   MAGICS++ STRING TYPE    PARAMETERS   */
	else if( !strcmp( param_type,"string" ) )
	  {
		mag_setc( param_name, param_value );		
	  }

    	/*   MAGICS++ STRINGARRAY  TYPE    PARAMETERS   */
	else if( !strcmp( param_type,"stringarray" ) )
	  {
		if( DBG )
		  fprintf(stderr, "Input strarr is %s  Sep char is %s Search char is %s\n",param_value , sep_char, search_char );
		if( strstr( param_value,";") )
		{
		  sep_char = ";";
		}
		
		if( DBG )
		  fprintf( stderr, "Input strarr is %s  Sep char is %s\n",param_value , sep_char );
		split_str_count = StringSplitWithSeperator( param_value, sep_char, &split_str );
		
		if( DBG )
		  fprintf( stderr, "Input strarr is %s split str count is %d Sep char is %s\n",param_value, split_str_count, sep_char );
		
		mag_set1c( param_name, (const char**)split_str, split_str_count );		
		free( split_str );
	  }
	else 
	  {
		ret_flag = 3;
		fprintf(stderr, "Unknown Parameter Type\n" );
	  }

	return ret_flag;
}
