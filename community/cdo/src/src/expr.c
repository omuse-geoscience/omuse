#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <ctype.h>

#include <cdi.h>
#include "cdo.h"
#include "cdo_int.h"
#include "field.h"
#include "expr.h"
#include "expr_yacc.h"

#define    COMPLT(x,y)  ((x) < (y) ? 1 : 0)
#define    COMPGT(x,y)  ((x) > (y) ? 1 : 0)
#define    COMPLE(x,y)  ((x) <= (y) ? 1 : 0)
#define    COMPGE(x,y)  ((x) >= (y) ? 1 : 0)
#define    COMPNE(x,y)  (IS_NOT_EQUAL(x,y) ? 1 : 0)
#define    COMPEQ(x,y)  (IS_EQUAL(x,y) ? 1 : 0)
#define   COMPLEG(x,y)  ((x) < (y) ? -1 : ((x) > (y) ? 1 : 0))
#define   COMPAND(x,y)  (IS_NOT_EQUAL(x,0) && IS_NOT_EQUAL(y,0) ? 1 : 0)
#define    COMPOR(x,y)  (IS_NOT_EQUAL(x,0) || IS_NOT_EQUAL(y,0) ? 1 : 0)
#define  MVCOMPLT(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPLT(x,y))
#define  MVCOMPGT(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPGT(x,y))
#define  MVCOMPLE(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPLE(x,y))
#define  MVCOMPGE(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPGE(x,y))
#define  MVCOMPNE(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPNE(x,y))
#define  MVCOMPEQ(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPEQ(x,y))
#define MVCOMPLEG(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPLEG(x,y))
#define MVCOMPAND(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPAND(x,y))
#define  MVCOMPOR(x,y)  (DBL_IS_EQUAL((x),missval1) ? missval1 : COMPOR(x,y))

static double f_int(double x)  { return ((int)(x)); }
static double f_nint(double x) { return (round(x)); }
static double f_sqr(double x)  { return (x*x);      }

typedef struct {
  int type;
  char *name;                      /* function name            */
  double (*func)(double);          /* pointer to function      */
}
func_t;

static func_t fun_sym_tbl[] =
{
  /* scalar functions */
  {0, "abs",   fabs},
  {0, "floor", floor},
  {0, "ceil",  ceil},
  {0, "int",   f_int},
  {0, "nint",  f_nint},
  {0, "sqr",   f_sqr},
  {0, "sqrt",  sqrt},
  {0, "exp",   exp},
  {0, "erf",   erf},
  {0, "log",   log},
  {0, "log10", log10},
  {0, "sin",   sin},
  {0, "cos",   cos},
  {0, "tan",   tan},
  {0, "sinh",  sinh},
  {0, "cosh",  cosh},
  {0, "tanh",  tanh},
  {0, "asin",  asin},
  {0, "acos",  acos},
  {0, "atan",  atan},
  {0, "asinh", asinh},
  {0, "acosh", acosh},
  {0, "atanh", atanh},
  {0, "gamma", tgamma},

  /* array functions
  {1, "min",   min},
  {1, "max",   max},
  {1, "sum",   sum},
  {1, "avg",   avg},
  {1, "mean",  mean},
  {1, "std",   std},
  {1, "var",   var},
  */
};

static int NumFunc = sizeof(fun_sym_tbl) / sizeof(fun_sym_tbl[0]);

static
nodeType *expr_con_con(int oper, nodeType *p1, nodeType *p2)
{
  nodeType *p = (nodeType*) malloc(sizeof(nodeType));

  p->type = typeCon;

  double cval1 = p1->u.con.value;
  double cval2 = p2->u.con.value;

  switch ( oper )
    {
    case '+':  cval1 = cval1 + cval2; break;
    case '-':  cval1 = cval1 - cval2; break;
    case '*':  cval1 = cval1 * cval2; break;
    case '/':  cval1 = cval1 / cval2; break;
    case '^':  cval1 = pow(cval1, cval2); break;
    default:   cdoAbort("%s: operator %c unsupported!", __func__, oper); break;
    }

  p->u.con.value = cval1;

  return (p);
}

static
nodeType *expr_con_var(int oper, nodeType *p1, nodeType *p2)
{
  int gridID   = p2->gridID;
  int zaxisID  = p2->zaxisID;
  int nmiss    = p2->nmiss;
  double missval1 = p2->missval;
  double missval2 = p2->missval;

  int ngp  = gridInqSize(gridID);
  int nlev = zaxisInqSize(zaxisID);
  long n   = ngp*nlev;
  long i;

  nodeType *p = (nodeType*) malloc(sizeof(nodeType));

  p->type     = typeVar;
  p->tmpvar   = 1;
  p->u.var.nm = strdupx("tmp");
  p->gridID   = gridID;
  p->zaxisID  = zaxisID;
  p->missval  = missval1;

  p->data = (double*) malloc(n*sizeof(double));
  double *restrict odat = p->data;
  const double *restrict idat = p2->data;
  double cval = p1->u.con.value;

  switch ( oper )
    {
    case '+':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = ADD(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] = cval + idat[i];
      break;
    case '-':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = SUB(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] = cval - idat[i];
      break;
    case '*':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MUL(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] = cval * idat[i];
      break;
    case '/':
      for ( i=0; i<n; ++i ) odat[i] = DIV(cval, idat[i]);
      break;
    case '^':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = POW(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] = pow(cval, idat[i]);
      break;
    case '<':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPLT(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPLT(cval, idat[i]);
      break;
    case '>':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPGT(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPGT(cval, idat[i]);
      break;
    case LE:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPLE(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPLE(cval, idat[i]);
      break;
    case GE:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPGE(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPGE(cval, idat[i]);
      break;
    case NE:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPNE(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPNE(cval, idat[i]);
      break;
    case EQ:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPEQ(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPEQ(cval, idat[i]);
      break;
    case LEG:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPLEG(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPLEG(cval, idat[i]);
      break;
    case AND:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPAND(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPAND(cval, idat[i]);
    case OR:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPOR(cval, idat[i]);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPOR(cval, idat[i]);
      break;
    default:
      cdoAbort("%s: operator %c unsupported!", __func__, oper);
      break;
    }

  nmiss = 0;
  for ( i = 0; i < n; i++ )
    if ( DBL_IS_EQUAL(p->data[i], missval1) ) nmiss++;

  p->nmiss = nmiss;

  if ( p2->tmpvar ) free(p2->data);

  return (p);
}

static
nodeType *expr_var_con(int oper, nodeType *p1, nodeType *p2)
{
  int gridID   = p1->gridID;
  int zaxisID  = p1->zaxisID;
  int nmiss    = p1->nmiss;
  double missval1 = p1->missval;
  double missval2 = p1->missval;

  int ngp  = gridInqSize(gridID);
  int nlev = zaxisInqSize(zaxisID);
  long n   = ngp*nlev;
  long i;

  nodeType *p = (nodeType*) malloc(sizeof(nodeType));

  p->type     = typeVar;
  p->tmpvar   = 1;
  p->u.var.nm = strdupx("tmp");
  p->gridID   = gridID;
  p->zaxisID  = zaxisID;
  p->missval  = missval1;

  p->data = (double*) malloc(n*sizeof(double));
  double *restrict odat = p->data;
  const double *restrict idat = p1->data;
  double cval = p2->u.con.value;

  switch ( oper )
    {
    case '+':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = ADD(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] = idat[i] + cval;
      break;
    case '-':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = SUB(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] = idat[i] - cval;
      break;
    case '*':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MUL(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] = idat[i] * cval;
      break;
    case '/':
      if ( nmiss || IS_EQUAL(cval, 0) ) for ( i=0; i<n; ++i ) odat[i] = DIV(idat[i], cval);
      else                              for ( i=0; i<n; ++i ) odat[i] = idat[i] / cval;
      break;
    case '^':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = POW(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] = pow(idat[i], cval);
      break;
    case '<':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPLT(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPLT(idat[i], cval);
      break;
    case '>':
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPGT(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPGT(idat[i], cval);
      break;
    case LE:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPLE(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPLE(idat[i], cval);
      break;
    case GE:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPGE(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPGE(idat[i], cval);
      break;
    case NE:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPNE(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPNE(idat[i], cval);
      break;
    case EQ:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPEQ(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPEQ(idat[i], cval);
      break;
    case LEG:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPLEG(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPLEG(idat[i], cval);
      break;
    case AND:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPAND(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPAND(idat[i], cval);
      break;
    case OR:
      if ( nmiss ) for ( i=0; i<n; ++i ) odat[i] = MVCOMPOR(idat[i], cval);
      else         for ( i=0; i<n; ++i ) odat[i] =   COMPOR(idat[i], cval);
      break;
    default:
      cdoAbort("%s: operator %c unsupported!", __func__, oper);
      break;
    }

  nmiss = 0;
  for ( i = 0; i < n; i++ )
    if ( DBL_IS_EQUAL(p->data[i], missval1) ) nmiss++;

  p->nmiss = nmiss;

  if ( p1->tmpvar ) free(p1->data);

  return (p);
}

static
nodeType *expr_var_var(int oper, nodeType *p1, nodeType *p2)
{
  long i;
  long nlev, k;
  long loff, loff1, loff2;
  int nmiss;

  int nmiss1 = p1->nmiss;
  int nmiss2 = p2->nmiss;
  double missval1 = p1->missval;
  double missval2 = p2->missval;

  long ngp1 = gridInqSize(p1->gridID);
  long ngp2 = gridInqSize(p2->gridID);

  if ( ngp1 != ngp2 ) cdoAbort("Number of grid points differ. ngp1 = %ld, ngp2 = %ld", ngp1, ngp2);

  long ngp = ngp1;

  long nlev1 = zaxisInqSize(p1->zaxisID);
  long nlev2 = zaxisInqSize(p2->zaxisID);

  nodeType *p = (nodeType*) malloc(sizeof(nodeType));

  p->type     = typeVar;
  p->tmpvar   = 1;
  p->u.var.nm = strdupx("tmp");

  if ( nlev1 > nlev2 )
    {
      nlev = nlev1;
      p->gridID  = p1->gridID;
      p->zaxisID = p1->zaxisID;
      p->missval = p1->missval;
      if ( nlev2 != 1 ) cdoAbort("nlev2 = %d must be 1!", nlev2);
    }
  else if ( nlev2 > nlev1 )
    {
      nlev = nlev2;
      p->gridID  = p2->gridID;
      p->zaxisID = p2->zaxisID;
      p->missval = p2->missval;
      if ( nlev1 != 1 ) cdoAbort("nlev1 = %d must be 1!", nlev1);
    }
  else
    {
      nlev = nlev1;
      p->gridID  = p1->gridID;
      p->zaxisID = p1->zaxisID;
      p->missval = p1->missval;
    }

  p->data = (double*) malloc(ngp*nlev*sizeof(double));

  for ( k = 0; k < nlev; k++ )
    {
      loff = k*ngp;

      if ( nlev1 == 1 ) loff1 = 0;
      else              loff1 = k*ngp;

      if ( nlev2 == 1 ) loff2 = 0;
      else              loff2 = k*ngp;

      const double *restrict idat1 = p1->data+loff1;
      const double *restrict idat2 = p2->data+loff2;
      double *restrict odat = p->data+loff;
      int nmiss = nmiss1 > 0 || nmiss2 > 0;

      switch ( oper )
	{
	case '+':
	  if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = ADD(idat1[i], idat2[i]);
	  else         for ( i=0; i<ngp; ++i ) odat[i] = idat1[i] + idat2[i];
	  break;
	case '-':
	  if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = SUB(idat1[i], idat2[i]);
	  else         for ( i=0; i<ngp; ++i ) odat[i] = idat1[i] - idat2[i];
	  break;
	case '*':
	  if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MUL(idat1[i], idat2[i]);
	  else         for ( i=0; i<ngp; ++i ) odat[i] = idat1[i] * idat2[i];
	  break;
	case '/':
	  if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = DIV(idat1[i], idat2[i]);
	  else
	    {
	      for ( i = 0; i < ngp; ++i )
		{
		  if ( IS_EQUAL(idat2[i], 0.) ) odat[i] = missval1;
		  else                          odat[i] = idat1[i] / idat2[i];
		}
	    }
	  break;
	case '^':
	  if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = POW(idat1[i], idat2[i]);
	  else         for ( i=0; i<ngp; ++i ) odat[i] = pow(idat1[i], idat2[i]);
	  break;
	case '<':
	  if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPLT(idat1[i], idat2[i]);
	  else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPLT(idat1[i], idat2[i]);
	  break;
	case '>':
	  if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPGT(idat1[i], idat2[i]);
	  else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPGT(idat1[i], idat2[i]);
	  break;
	case LE:
	  if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPLE(idat1[i], idat2[i]);
	  else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPLE(idat1[i], idat2[i]);
	  break;
	case GE:
	  if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPGE(idat1[i], idat2[i]);
	  else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPGE(idat1[i], idat2[i]);
	  break;
	case NE:
	  if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPNE(idat1[i], idat2[i]);
	  else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPNE(idat1[i], idat2[i]);
	  break;
	case EQ:
	  if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPEQ(idat1[i], idat2[i]);
	  else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPEQ(idat1[i], idat2[i]);
	  break;
	case LEG:
	  if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPLEG(idat1[i], idat2[i]);
	  else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPLEG(idat1[i], idat2[i]);
	  break;
	case AND:
	  if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPAND(idat1[i], idat2[i]);
	  else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPAND(idat1[i], idat2[i]);
	  break;
	case OR:
	  if ( nmiss ) for ( i=0; i<ngp; ++i ) odat[i] = MVCOMPOR(idat1[i], idat2[i]);
	  else         for ( i=0; i<ngp; ++i ) odat[i] =   COMPOR(idat1[i], idat2[i]);
	  break;
	default:
	  cdoAbort("%s: operator %d (%c) unsupported!", __func__, (int)oper, oper);
          break;
	}
    }

  nmiss = 0;
  for ( i = 0; i < ngp*nlev; i++ )
    if ( DBL_IS_EQUAL(p->data[i], missval1) ) nmiss++;

  p->nmiss = nmiss;

  if ( p1->tmpvar ) free(p1->data);
  if ( p2->tmpvar ) free(p2->data);

  return (p);
}

static
void ex_copy(nodeType *p2, nodeType *p1)
{
  long i;

  if ( cdoVerbose ) printf("\tcopy %s\n", p1->u.var.nm);

  long ngp1 = gridInqSize(p1->gridID);
  long ngp2 = gridInqSize(p2->gridID);

  if ( ngp1 != ngp2 )
    cdoAbort("Number of grid points differ. ngp1 = %d, ngp2 = %d", ngp1, ngp2);

  long ngp = ngp2;
  long nlev = zaxisInqSize(p2->zaxisID);

  for ( i = 0; i < ngp*nlev; i++ ) p2->data[i] = p1->data[i];

  p2->missval = p1->missval;
  p2->nmiss   = p1->nmiss;
}

static
nodeType *expr(int oper, nodeType *p1, nodeType *p2)
{
  nodeType *p = NULL;

  if ( p1->type == typeVar && p2->type == typeVar )
    {
      p = expr_var_var(oper, p1, p2);
      if ( cdoVerbose )
	printf("\t%s %c %s\n", p1->u.var.nm, oper, p2->u.var.nm);
    }
  else if ( p1->type == typeCon && p2->type == typeCon )
    {
      p = expr_con_con(oper, p1, p2);
      if ( cdoVerbose )
	printf("\t%g %c %g\n", p1->u.con.value, oper, p2->u.con.value);
    }
  else if ( p1->type == typeVar && p2->type == typeCon )
    {
      p = expr_var_con(oper, p1, p2);
      if ( cdoVerbose )
	printf("\t%s %c %g\n", p1->u.var.nm, oper, p2->u.con.value);
    }
  else if ( p1->type == typeCon && p2->type == typeVar )
    {
      p = expr_con_var(oper, p1, p2);
      if ( cdoVerbose )
	printf("\t%g %c %s\n", p1->u.con.value, oper, p2->u.var.nm);
    }
  else
    cdoAbort("Internal problem!");

  return (p);
}

static
nodeType *ex_fun_con(char *fun, nodeType *p1)
{
  int i;
  int funcID = -1;

  nodeType *p = (nodeType*) malloc(sizeof(nodeType));

  p->type = typeCon;

  for ( i = 0; i < NumFunc; i++)
    if ( fun_sym_tbl[i].type == 0 )
      if ( strcmp(fun, fun_sym_tbl[i].name) == 0 )
	{ 
	  funcID = i;
	  break;
	}

  if ( funcID == -1 )
    cdoAbort("Function >%s< not available!", fun);

  p->u.con.value = fun_sym_tbl[funcID].func(p1->u.con.value);

  return (p);
}

static
nodeType *ex_fun_var(char *fun, nodeType *p1)
{
  long i;
  int funcID = -1;
  int gridID  = p1->gridID;
  int zaxisID = p1->zaxisID;
  int nmiss   = p1->nmiss;
  double missval = p1->missval;

  long ngp  = gridInqSize(gridID);
  long nlev = zaxisInqSize(zaxisID);

  nodeType *p = (nodeType*) malloc(sizeof(nodeType));

  p->type     = typeVar;
  p->tmpvar   = 1;
  p->u.var.nm = strdupx("tmp");
  p->gridID   = gridID;
  p->zaxisID  = zaxisID;
  p->missval  = missval;

  p->data = (double*) malloc(ngp*nlev*sizeof(double));

  for ( i = 0; i < NumFunc; i++)
    if ( strcmp(fun, fun_sym_tbl[i].name) == 0 )
      { 
	funcID = i;
	break;
      }

  if ( funcID == -1 )
    cdoAbort("Function >%s< not available!", fun);

  if ( nmiss > 0 )
    {
      for ( i = 0; i < ngp*nlev; i++ )
	{
	  errno = -1;
	  p->data[i] = DBL_IS_EQUAL(p1->data[i], missval) ? missval : fun_sym_tbl[funcID].func(p1->data[i]);
	  if ( errno == EDOM || errno == ERANGE ) p->data[i] = missval;
	  else if ( isnan(p->data[i]) )  p->data[i] = missval;
	}
    }
  else
    {
      for ( i = 0; i < ngp*nlev; i++ )
	{
	  errno = -1;
	  p->data[i] = fun_sym_tbl[funcID].func(p1->data[i]);
	  if ( errno == EDOM || errno == ERANGE ) p->data[i] = missval;
	  else if ( isnan(p->data[i]) )  p->data[i] = missval;
	}
    }

  nmiss = 0;
  for ( i = 0; i < ngp*nlev; i++ )
    if ( DBL_IS_EQUAL(p->data[i], missval) ) nmiss++;

  p->nmiss = nmiss;

  if ( p1->tmpvar ) free(p1->data);

  return (p);
}

static
nodeType *ex_fun(char *fun, nodeType *p1)
{
  nodeType *p = NULL;

  if ( p1->type == typeVar )
    {
      p = ex_fun_var(fun, p1);
      if ( cdoVerbose ) printf("\t%s (%s)\n", fun, p1->u.var.nm);
    }
  else if ( p1->type == typeCon )
    {
      p = ex_fun_con(fun, p1);
      if ( cdoVerbose ) printf("\t%s (%g)\n", fun, p1->u.con.value);
    }
  else
    cdoAbort("Internal problem!");

  return (p);
}

static
nodeType *ex_uminus_var(nodeType *p1)
{
  int gridID  = p1->gridID;
  int zaxisID = p1->zaxisID;
  int nmiss   = p1->nmiss;
  double missval = p1->missval;

  long ngp  = gridInqSize(gridID);
  long nlev = zaxisInqSize(zaxisID);

  nodeType *p = (nodeType*) malloc(sizeof(nodeType));

  p->type     = typeVar;
  p->tmpvar   = 1;
  p->u.var.nm = strdupx("tmp");
  p->gridID   = gridID;
  p->zaxisID  = zaxisID;
  p->missval  = missval;

  p->data = (double*) malloc(ngp*nlev*sizeof(double));

  if ( nmiss > 0 )
    {
      for ( long i = 0; i < ngp*nlev; i++ )
	p->data[i] = DBL_IS_EQUAL(p1->data[i], missval) ? missval : -(p1->data[i]);
    }
  else
    {
      for ( long i = 0; i < ngp*nlev; i++ )
	p->data[i] = -(p1->data[i]);
    }

  p->nmiss = nmiss;
  
  return (p);
}

static
nodeType *ex_uminus_con(nodeType *p1)
{
  nodeType *p = (nodeType*) malloc(sizeof(nodeType));

  p->type = typeCon;

  p->u.con.value = -(p1->u.con.value);

  return (p);
}

static
nodeType *ex_uminus(nodeType *p1)
{
  nodeType *p = NULL;

  if ( p1->type == typeVar )
    {
      p = ex_uminus_var(p1);
      if ( cdoVerbose ) printf("\t- (%s)\n", p1->u.var.nm);
    }
  else if ( p1->type == typeCon )
    {
      p = ex_uminus_con(p1);
      if ( cdoVerbose ) printf("\t- (%g)\n", p1->u.con.value);
    }
  else
    cdoAbort("Internal problem!");

  return (p);
}

static
nodeType *ex_ifelse(nodeType *p1, nodeType *p2, nodeType *p3)
{
  if ( cdoVerbose ) printf("\t %s ? %s : %s\n", p1->u.var.nm, p2->u.var.nm, p3->u.var.nm);

  if ( p1->type == typeCon ) cdoAbort("expr?expr:expr: First expression is a constant but must be a variable!");

  int nmiss1 = p1->nmiss;
  long ngp1 = gridInqSize(p1->gridID);
  long nlev1 = zaxisInqSize(p1->zaxisID);
  double missval1 = p1->missval;
  double *pdata1 = p1->data;

  long ngp = ngp1;
  long nlev = nlev1;
  nodeType *px = p1;

  double missval2 = missval1;
  double *pdata2;
  long ngp2 = 1;
  long nlev2 = 1;
  
  if ( p2->type == typeCon )
    {
      pdata2 = &p2->u.con.value;
    }
  else
    {
      ngp2 = gridInqSize(p2->gridID);
      nlev2 = zaxisInqSize(p2->zaxisID);
      missval2 = p2->missval;
      pdata2 = p2->data;
      if ( ngp2 > 1 && ngp2 != ngp1 )
	cdoAbort("expr?expr:expr: Number of grid points differ. ngp1 = %ld, ngp2 = %ld", ngp1, ngp2);
      if ( nlev2 > 1 && nlev2 != nlev )
	{
	  if ( nlev == 1 )
	    {
	      nlev = nlev2;
	      px = p2;
	    }
	  else
	    cdoAbort("expr?expr:expr: Number of levels differ. nlev = %ld, nlev2 = %ld", nlev, nlev2);
	}
    }

  double missval3 = missval1;
  double *pdata3;
  long ngp3 = 1;
  long nlev3 = 1;
  
  if ( p3->type == typeCon )
    {
      pdata3 = &p3->u.con.value;
    }
  else
    {
      ngp3 = gridInqSize(p3->gridID);
      nlev3 = zaxisInqSize(p3->zaxisID);
      missval3 = p3->missval;
      pdata3 = p3->data;
      if ( ngp3 > 1 && ngp3 != ngp1 )
	cdoAbort("expr?expr:expr: Number of grid points differ. ngp1 = %ld, ngp3 = %ld", ngp1, ngp3);
      if ( nlev3 > 1 && nlev3 != nlev )
	{
	  if ( nlev == 1 )
	    {
	      nlev = nlev3;
	      px = p3;
	    }
	  else
	    cdoAbort("expr?expr:expr: Number of levels differ. nlev = %ld, nlev3 = %ld", nlev, nlev3);
	}
    }

  nodeType *p = (nodeType*) malloc(sizeof(nodeType));

  p->type     = typeVar;
  p->tmpvar   = 1;
  p->u.var.nm = strdupx("tmp");

  p->gridID  = px->gridID;
  p->zaxisID = px->zaxisID;
  p->missval = px->missval;

  p->data = (double*) malloc(ngp*nlev*sizeof(double));

  long loff, loff1, loff2, loff3;

  for ( long k = 0; k < nlev; ++k )
    {
      loff = k*ngp;

      if ( nlev1 == 1 ) loff1 = 0;
      else              loff1 = k*ngp;

      if ( nlev2 == 1 ) loff2 = 0;
      else              loff2 = k*ngp;

      if ( nlev3 == 1 ) loff3 = 0;
      else              loff3 = k*ngp;

      const double *restrict idat1 = pdata1+loff1;
      const double *restrict idat2 = pdata2+loff2;
      const double *restrict idat3 = pdata3+loff3;
      double *restrict odat = p->data+loff;

      double ival2 = idat2[0];
      double ival3 = idat3[0];
      for ( long i = 0; i < ngp; ++i ) 
	{
	  if ( ngp2 > 1 ) ival2 = idat2[i];
	  if ( ngp3 > 1 ) ival3 = idat3[i];

	  if ( nmiss1 && DBL_IS_EQUAL(idat1[i], missval1) )
	    odat[i] = missval1;
	  else if ( IS_NOT_EQUAL(idat1[i], 0) )
	    odat[i] = DBL_IS_EQUAL(ival2, missval2) ? missval1 : ival2;
	  else
	    odat[i] = DBL_IS_EQUAL(ival3, missval3) ? missval1 : ival3;
	}
    }

  return (p);
}


int exNode(nodeType *p, parse_parm_t *parse_arg)
{
  if ( ! p ) return(0);

  /* node is leaf */
  if ( p->type == typeCon || p->type == typeVar || p->u.opr.nops == 0 )
    {
      return (0);
    }

  /* node has children */
  for ( int k = 0; k < p->u.opr.nops; k++ )
    {
      exNode(p->u.opr.op[k], parse_arg);
    }

  return (0);
}


nodeType *expr_run(nodeType *p, parse_parm_t *parse_arg)
{
  int gridID1 = -1, zaxisID1 = -1, tsteptype1 = -1;
  double missval = 0;
  char varname[256];
  int varID, nvars;
  nodeType *rnode = NULL;

  if ( ! p ) return (rnode);

  /*  if ( ! parse_arg->init ) { exNode(p, parse_arg); return (0); } */

  switch ( p->type )
    {
    case typeCon:       
      if ( parse_arg->init )
	{
	  if ( parse_arg->debug )
	    printf("\tpush const \t%g\n", p->u.con.value);
	}
      else
	{
	  rnode = p;
	}

      break;
    case typeVar:
      /*    if ( parse_arg->init ) */
	{
	  if ( parse_arg->debug )
	    printf("\tpush var \t%s\n", p->u.var.nm);

	  nvars = vlistNvars(parse_arg->vlistID1);
	  for ( varID = 0; varID < nvars; varID++ )
	    {
	      vlistInqVarName(parse_arg->vlistID1, varID, varname);
	      if ( strcmp(varname, p->u.var.nm) == 0 ) break;
	    }

	  if ( varID == nvars )
	    {
	      cdoAbort("Variable >%s< not found!", p->u.var.nm);
	    }
	  else
	    {
	      int nlev1, nlev2 = 0;
	      if ( varID >= MAX_VARS ) cdoAbort("Too many parameter (limit=%d)!", MAX_VARS);

	      if ( parse_arg->var_needed[varID] == 0 )
		{

		  parse_arg->var[varID] = strdupx(p->u.var.nm);
		  parse_arg->varID[varID] = varID;
		  parse_arg->var_needed[varID] = 1;
		}

	      gridID1    = vlistInqVarGrid(parse_arg->vlistID1, varID);
	      zaxisID1   = vlistInqVarZaxis(parse_arg->vlistID1, varID);
	      tsteptype1 = vlistInqVarTsteptype(parse_arg->vlistID1, varID);
	      missval    = vlistInqVarMissval(parse_arg->vlistID1, varID);
	      nlev1 = zaxisInqSize(zaxisID1);

	      parse_arg->missval2 = missval;

	      if ( parse_arg->gridID2 == -1 )
		parse_arg->gridID2 = gridID1;

	      if ( parse_arg->zaxisID2 != -1 ) nlev2 = zaxisInqSize(parse_arg->zaxisID2);

	      if ( parse_arg->zaxisID2 == -1 || (nlev1 > 1 && nlev2 == 1) )
		parse_arg->zaxisID2 = zaxisID1;

	      if ( parse_arg->tsteptype2 == -1 || parse_arg->tsteptype2 == TSTEP_CONSTANT )
		parse_arg->tsteptype2 = tsteptype1;
	    }
	}
	/* else */
	{ 
	  if ( parse_arg->debug )
	    printf("var: %s %d %d %d\n", p->u.var.nm, varID, gridID1, zaxisID1);
	  p->gridID  = gridID1;
	  p->zaxisID = zaxisID1;
	  p->missval = missval;
          p->nmiss   = 0;
	  if ( ! parse_arg->init )
	    {
	      p->data  = parse_arg->vardata1[varID];
	      p->nmiss = parse_arg->nmiss[varID];
	    }
	  p->tmpvar  = 0;
	  rnode = p;
	}

      break;
    case typeFun:
      if ( parse_arg->init )
	{
	  expr_run(p->u.fun.op, parse_arg);

	  if ( parse_arg->debug )
	    printf("\tcall \t%s\n", p->u.fun.name);
	}
      else
	{
	  rnode = ex_fun(p->u.fun.name, expr_run(p->u.fun.op, parse_arg));
	}
      break;
    case typeOpr:
      switch( p->u.opr.oper )
	{
        case '=':
	  parse_arg->gridID2    = -1;
	  parse_arg->zaxisID2   = -1;
          parse_arg->tsteptype2 = -1;

	  rnode = expr_run(p->u.opr.op[1], parse_arg);

	  if ( parse_arg->init )
	    {
	      if ( parse_arg->debug )
		printf("\tpop  var \t%s\n", p->u.opr.op[0]->u.var.nm);
	      /*
	      if ( p->u.opr.op[1]->type != typeVar )
		cdoAbort("Operand not variable!");
	      */
	      if ( parse_arg->gridID2 == -1 || parse_arg->zaxisID2 == -1 || parse_arg->tsteptype2 == -1 )
		cdoAbort("Operand not variable!");

	      varID = vlistDefVar(parse_arg->vlistID2, parse_arg->gridID2, parse_arg->zaxisID2, parse_arg->tsteptype2);
	      const char *varname = p->u.opr.op[0]->u.var.nm;
	      vlistDefVarName(parse_arg->vlistID2, varID, varname);
	      vlistDefVarMissval(parse_arg->vlistID2, varID, parse_arg->missval2);
	      if ( memcmp(varname, "var", 3) == 0 )
		{
		  if ( strlen(varname) > 3 && isdigit(varname[3]) )
		    {
		      int code = atoi(varname+3);
		      vlistDefVarCode(parse_arg->vlistID2, varID, code);
		    }
		}
	    }
	  else
	    {
	      if ( parse_arg->debug )
		printf("\tpop var\t%s\t%s\n", p->u.opr.op[0]->u.var.nm, rnode->u.var.nm);

	      nvars = vlistNvars(parse_arg->vlistID2);
	      for ( varID = nvars-1; varID >= 0; varID-- )
		{
		  vlistInqVarName(parse_arg->vlistID2, varID, varname);
		  if ( strcmp(varname, p->u.opr.op[0]->u.var.nm) == 0 ) break;
		}

	      if ( varID < 0 )
		{
		  cdoAbort("Variable >%s< not found!", p->u.opr.op[0]->u.var.nm);
		}
	      else
		{
		  parse_arg->gridID2  = vlistInqVarGrid(parse_arg->vlistID2, varID);
		  parse_arg->zaxisID2 = vlistInqVarZaxis(parse_arg->vlistID2, varID);
		  parse_arg->tsteptype2 = vlistInqVarTsteptype(parse_arg->vlistID2, varID);
		  missval  = vlistInqVarMissval(parse_arg->vlistID2, varID);
	      
		  p->gridID  = parse_arg->gridID2;
		  p->zaxisID = parse_arg->zaxisID2;
		  p->missval = missval;
		  p->data    = parse_arg->vardata2[varID];
		  p->tmpvar  = 0;

		  ex_copy(p, rnode);

		  if ( rnode->tmpvar ) free(rnode->data);
		}
	    }

	  break;
        case UMINUS:    
	  if ( parse_arg->init )
	    {
	      expr_run(p->u.opr.op[0], parse_arg);

	      if ( parse_arg->debug )
		printf("\tneg\n");
	    }
	  else
	    {
	      rnode = ex_uminus(expr_run(p->u.opr.op[0], parse_arg));
	    }

	  break;
        case '?':    
	  if ( parse_arg->init )
	    {
	      expr_run(p->u.opr.op[0], parse_arg);
	      expr_run(p->u.opr.op[1], parse_arg);
	      expr_run(p->u.opr.op[2], parse_arg);

	      if ( parse_arg->debug )
		printf("\t?:\n");
	    }
	  else
	    {
	      rnode = ex_ifelse(expr_run(p->u.opr.op[0], parse_arg),
			        expr_run(p->u.opr.op[1], parse_arg),
			        expr_run(p->u.opr.op[2], parse_arg));
	    }

	  break;
        default:
	  if ( parse_arg->init )
	    {
	      expr_run(p->u.opr.op[0], parse_arg);
	      expr_run(p->u.opr.op[1], parse_arg);
	      if ( parse_arg->debug )
		switch( p->u.opr.oper )
		  {
		  case '+':  printf("\tadd\n"); break;
		  case '-':  printf("\tsub\n"); break;
		  case '*':  printf("\tmul\n"); break;
		  case '/':  printf("\tdiv\n"); break;
		  case '<':  printf("\tcompLT\n"); break;
		  case '>':  printf("\tcompGT\n"); break;
		  case LE:   printf("\tcompLE\n"); break;
		  case GE:   printf("\tcompGE\n"); break;
		  case NE:   printf("\tcompNE\n"); break;
		  case EQ:   printf("\tcompEQ\n"); break;
		  case LEG:  printf("\tcompLEG\n"); break;
		  case AND:  printf("\tcompAND\n"); break;
		  case OR:   printf("\tcompOR\n"); break;
		  }
	    }
	  else
	    {
	      rnode = expr(p->u.opr.oper, expr_run(p->u.opr.op[0], parse_arg),
			                  expr_run(p->u.opr.op[1], parse_arg));
	    }
          break;
        }
      break;
    }

  return (rnode);
}
