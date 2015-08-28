#include "template_parser.h"
#include "magics_template_parser.h"
#include "results_template_parser.h"

#define DBG_MSG 0 


/* extern int GetMagicsParameterInfo( const char *user_name, char **magics_name, char **magics_type ); */

extern int GetMagicsParameterInfo(  char *user_name, xmlChar *param_value );

extern xmlNode *results_node;


/* Recursive function that sets the results parameters from the XML structure */

int results_template_parser( xmlNode * a_node, const char *varname ) 

{
    xmlNode *cur_node = NULL;
    xmlAttrPtr attr = NULL;
    xmlChar    *param_name,*param_value,*value;
    char       *param_type;

	
    if( a_node == NULL )
	return 1;

    if( !strcmp( a_node->name, "results" ) )
    {
 	value = xmlGetProp( a_node, "version" );

	if( value )
	{
    		if( DBG_MSG )
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
		
	    if( DBG_MSG )
            	printf( "Node Name: %s \n", cur_node->name );

	    if( cur_node->properties == NULL )
	    {
		if( cur_node->children == NULL )
		{
			printf( "NO ATTRIBUTES!!!\n" );
		}
	    }
	    else
	    {
		
		/* 	Loop Over the attributes and get the corresponding
                   	Magics Parameter name and type, set the value 
		*/

#if 0
	      printf( "Finding varname = %s  result_name = %s\n", varname, xmlGetProp( cur_node,"name") );
#endif

	      if ( strcmp( varname, xmlGetProp( cur_node,"name" ) ) == 0 )
	      {
#if 0
	          printf( "Found varname = %s  result_name = %s\n", varname, xmlGetProp( cur_node,"name") );
#endif

		  for( attr = cur_node->properties; attr; attr = attr->next )
		  {	
		      if( attr != NULL )
		      {

			  param_value = xmlNodeGetContent( attr->children );

			  /* if( !GetMagicsParameterInfo( attr->name, &magics_param_name, &param_type ) ) */
			  if( !GetMagicsParameterInfo( (char *) attr->name, param_value ) )
			  {

#if 0
			      printf("Done corresponding Magics Parameter found!\n");
			      printf("Setting corresponding Magics Parameter %s and type %s!\n",magics_param_name, param_type );
	  		      fprintf(stderr, "param_value: %s\n", param_value);
			      SetMagicsParameterValue( magics_param_name, param_type, param_value );
#endif

			  }
			  else
			  {
#if 0
			      printf("No corresponding Magics Parameter found!\n");
#endif
			  }
		      }	
		  }	

		  break;
	      }
	      else
	      {
	   	  fprintf(stderr,"Var Name not matching resetting Magics Params!\n");
		  /* Call the Reset functions of all the features to Reset the magics params to default */
	      }
	    }
        }
    }
    return 0;
}
