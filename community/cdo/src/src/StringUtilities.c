#include "StringUtilities.h"

#define DBG 0

int StringSplitWithSeperator(  char *source_string, char *seperator, char*** ptr_split_string )

{
	char *duplicate_src = NULL, **temp_list = NULL , *temp_str = NULL, *saveptr;
	int n = 0, i;
	int str_len = 0, sep_count = 0;	

	str_len = strlen( source_string );

	if( !str_len )
	  return 0;

	if( DBG )
	  fprintf(stderr, "StringSplitWithSeperator Input str %s , seperator %s \n", source_string, seperator );
	
	duplicate_src = strdup(source_string);
        
	for( i = 0; i < str_len; i++ )
	  {
		 if( duplicate_src[i] == *seperator )  
		   sep_count++;
	  }	
 
	temp_list  = (char**) malloc( sizeof( char* ) * (sep_count+1));
	
	if( DBG )
	  fprintf(stderr, "Input str %s , seperator %s  sep count %d\n", duplicate_src, seperator, sep_count );
	
	while( ( temp_str = strtok_r( duplicate_src, seperator, &saveptr ) ) )
          {
	    temp_list[n] = temp_str;
            n++;
            duplicate_src = saveptr;
          }
          
        if( DBG )
	  {
	    for( i = 0; i <  n; i++ )
	      fprintf(stderr, "str  %s \n", temp_list[i] );  
	  }
	  
	*ptr_split_string = temp_list;

	return n;
}



int IsNumeric (const char *s)
{
    char *ptr;
    if (s == NULL || *s == '\0' || isspace(*s))
      return 0;
    
    strtod (s, &ptr);
    return *ptr == '\0';
}


void StrToUpperCase ( char *sPtr )
{
    while ( *sPtr != '\0' )
    {
      *sPtr = toupper ( ( unsigned char ) *sPtr );
      ++sPtr;
    }
}


void StrToLowerCase ( char *sPtr )
{
    while ( *sPtr != '\0' )
    {
      *sPtr = tolower ( ( unsigned char ) *sPtr );
      ++sPtr;
    }
}

/* To replace a single char with another single char in a given string */

void StrReplaceChar( char *str_in, char orig_char, char rep_char )

{
  
  char *ref = NULL;
  
  ref = str_in; 
  
  if( strchr( str_in, orig_char) == NULL )
    return;
  
  if( DBG )  
    fprintf( stderr,"Before %s\n", ref );        
  
  while( *str_in != '\0' )
    {
      if( *str_in == orig_char )
	  *str_in = rep_char;
      str_in++;
    } 
  if( DBG )  
    fprintf( stderr,"After %s\n", ref );      
  
  return  ;     
}

