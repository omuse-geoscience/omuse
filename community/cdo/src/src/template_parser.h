#ifndef TEMPLATE_PARSER_HH
#define TEMPLATE_PARSER_HH

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <locale.h>

#if defined(HAVE_LIBXML2)
#include<libxml/parser.h>
#include<libxml/tree.h>


int template_parser( char *Filename, const char *varname );
int init_XMLtemplate_parser( char *Filename );
int updatemagics_and_results_nodes( );
int quit_XMLtemplate_parser( );


#endif

#endif
