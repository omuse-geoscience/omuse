#include <stdio.h>

#ifndef fileno
int fileno(FILE *stream);
#endif

#ifndef strdupx
#ifndef strdup
char *strdup(const char *s);
#endif
#define strdupx  strdup
/*
#define strdupx(s)			          \
({					      	  \
   const char *__old = (s);			  \
   size_t __len = strlen(__old) + 1;		  \
   char *__new = malloc(__len);	  \
   (char *) memcpy(__new, __old, __len);	  \
})
*/
#endif


typedef enum { typeCon, typeVar, typeFun, typeOpr } nodeEnum;

/* constants */
typedef struct {
  double value;               /* value of constant */
} conNodeType;

/* variables */
typedef struct {
  char *nm;                   /* variable name */
} varNodeType;

/* functions */
typedef struct {
  char *name;                 /* function name */
  struct nodeTypeTag *op;     /* operand       */
} funNodeType;

/* operators */
typedef struct {
  int oper;                   /* operator              */
  int nops;                   /* number of operands    */
  struct nodeTypeTag *op[1];  /* operands (expandable) */
} oprNodeType;

typedef struct nodeTypeTag {
  int tmpvar;
  int gridID, zaxisID;
  int nmiss;
  double missval;
  double *data;
  nodeEnum type;              /* type of node */

  /* union must be last entry in nodeType */
  /* because operNodeType may dynamically increase */
  union {
    conNodeType con;          /* constants   */
    varNodeType var;          /* variables   */
    funNodeType fun;          /* functions   */
    oprNodeType opr;          /* operators   */
  } u;
} nodeType;

#define MAX_VARS 1024

typedef struct{ /* prs_sct */
  int    vlistID1, vlistID2;
  int    nvars1, nvars2;
  int    nmiss[MAX_VARS];
  int    varID[MAX_VARS];
  int    var_needed[MAX_VARS];
  char   *var[MAX_VARS];
  int    init;
  int    debug;
  int    gridID2;
  int    zaxisID2;
  int    tsteptype2;
  double missval2;
  double **vardata1, **vardata2;
} parse_parm_t;


typedef union{
    double cvalue;              /* constant value */
    char *varnm;                /* variable name  */
    char *fname;                /* function name  */
    nodeType *nPtr;             /* node pointer   */
} stype_t;


#define YYSTYPE        stype_t
#define YY_EXTRA_TYPE  parse_parm_t *

#define YY_DECL int yylex(YYSTYPE *yylval_param, parse_parm_t *parse_arg, void *yyscanner)
YY_DECL;

int  yyparse(parse_parm_t *parse_arg, void*);
void yyerror(void *parse_arg, void *scanner, char *errstr);

int  yylex_init(void **);
int  yylex_destroy(void *);
void yyset_extra(YY_EXTRA_TYPE, void *);
