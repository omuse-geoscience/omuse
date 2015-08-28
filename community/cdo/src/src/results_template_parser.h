#ifndef RESULTS_TEMPLATE_PARSER_HH
#define RESULTS_TEMPLATE_PARSER_HH

#include<stdio.h>
#include<string.h>
#include<stdlib.h>

#include<libxml/parser.h>
#include<libxml/tree.h>

 int results_template_parser( xmlNode * a_node, const char *varname ); 

#endif
