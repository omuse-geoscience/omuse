/* This source code is copied from PINGO version 1.5 */

/* ********************************** */
/* HEADER FOR PARALLEL EIGEN SOLUTION */
/*  -->SEE END OF ROUTINE             */
/* ********************************** */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "cdo_int.h"

#define  FNORM_PRECISION  1e-12
#define  MAX_JACOBI_ITER  12

int jacobi_1side(double **M, double *A, long n);
void annihilate_1side(double **M, long i, long j, long n);

int n_finished;

// global variables to handle environment settings 
double fnorm_precision;
int max_jacobi_iter;

/* **************************************** */
/* ENDOF HEADER FOR PARALLEL EIGEN SOLUTION */
/*  -->SEE END OF ROUTINE                   */
/* **************************************** */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>

#include "cdo.h"
#include "cdo_int.h"
#include "statistic.h"

#ifndef  M_2_SQRTPI
# define M_2_SQRTPI	1.12837916709551257390	/* 2/sqrt(pi) */
#endif
#ifndef  M_SQRT2
# define M_SQRT2	1.41421356237309504880	/* sqrt(2) */
#endif


void make_symmetric_matrix_triangular (double **a, int n, double *d, double *e);
double pythagoras (double a, double b);
void eigen_solution_of_triangular_matrix (double *d, double *e, int n, double **a, const char *prompt);
int lu_decomposition (double **a, int n, int *index, int *sign);
void lu_backsubstitution (double **a, int n, int *index, double *b);

// moved to statistic.h to make heap_sort accessible for other modules
// void heap_sort (double *eig_val, double **a, int n);

static double gamma_help_1 (double a, double x, const char *prompt);
static double gamma_help_2 (double a, double x, const char *prompt);
static double beta_help (double a, double b, double x, const char *prompt);


void eigen_solution_of_symmetric_matrix (double **a, double *eig_val,
					 int n, const char *prompt)
/* After return the rows (!!!) of a are the eigenvectors */
{
  double *e;
  int i, j;
  double temp;
  
  e = (double*) malloc(n * sizeof(double));
  
  make_symmetric_matrix_triangular (a, n, eig_val, e);
  
  eigen_solution_of_triangular_matrix (eig_val, e, n, a, prompt);
  
  free (e);

  for (i = 0; i < n; i++)
    for (j = i + 1; j < n; j++)
      {
        temp = a[i][j];
        a[i][j] = a[j][i];
        a[j][i] = temp;
      }
  
  heap_sort (eig_val, a, n);
}


void heap_sort (double *eig_val, double **a, int n)
{
  int i, j, j_next, k1, k2;
  double temp;
  double *temp_p;
  
  /* First part of heap sort: */
  for (i = n / 2 - 1; i >= 0; i--)
    {
      for (j = i; 2 * j + 1 < n; j = j_next)
        {
          k1 = 2 * j + 1;
          k2 = 2 * j + 2;
          j_next = j;
          if (eig_val[k1] < eig_val[j_next])
            j_next = k1;
          if (k2 < n && eig_val[k2] < eig_val[j_next])
            j_next = k2;
          if (j == j_next)
            break;
          temp = eig_val[j_next];
          eig_val[j_next] = eig_val[j];
          eig_val[j] = temp;
          temp_p = a[j_next];
          a[j_next] = a[j];
          a[j] = temp_p;
        }
    }
  /* Second part of head sort: */
  for (i = n - 1; i > 0; i--)
    {
      temp = eig_val[i];
      eig_val[i] = eig_val[0];
      eig_val[0] = temp;
      temp_p = a[i];
      a[i] = a[0];
      a[0] = temp_p;
      for (j = 0; 2 * j + 1 < i; j = j_next)
        {
          k1 = 2 * j + 1;
          k2 = 2 * j + 2;
          j_next = j;
          if (eig_val[k1] < eig_val[j_next])
            j_next = k1;
          if (k2 < i && eig_val[k2] < eig_val[j_next])
            j_next = k2;
          if (j == j_next)
            break;
          temp = eig_val[j_next];
          eig_val[j_next] = eig_val[j];
          eig_val[j] = temp;
          temp_p = a[j_next];
          a[j_next] = a[j];
          a[j] = temp_p;
        }
    }
}


void make_symmetric_matrix_triangular(double **a, int n, double *d, double *e)
{
  int i, j, k;
  double f, g, h, hh, scale;
 
  for ( i = n - 1; i >= 1; i-- )
    {
      h = scale = 0;
      if (i > 1)
        {
          for (k = 0; k < i; k++)
            scale += fabs (a[i][k]);
          if ( DBL_IS_EQUAL(scale, 0.) )
            e[i] = a[i][i - 1];
          else
            {
              for (k = 0; k < i; k++)
                {
                  a[i][k] /= scale;
                  h += a[i][k] * a[i][k];
                }
              f = a[i][i - 1];
              g = f >= 0 ? -sqrt (h) : sqrt (h);
              e[i] = scale * g;
              h -= f * g;
              a[i][i - 1] = f - g;
              f = 0;
              for (j = 0; j < i; j++)
                {
                  a[j][i] = a[i][j] / h;
                  g = 0;
                  for (k = 0; k <= j; k++)
                    g += a[j][k] * a[i][k];
                  for (k = j + 1; k < i; k++)
                    g += a[k][j] * a[i][k];
                  e[j] = g / h;
                  f += e[j] * a[i][j];
                }
              hh = f / (2 * h);
              for (j = 0; j < i; j++)
                {
                  f = a[i][j];
                  e[j] = g = e[j] - hh * f;
                  for (k = 0; k <= j; k++)
                    a[j][k] -= f * e[k] + g * a[i][k];
                }
            }
        }
      else
        e[i] = a[i][i - 1];

      d[i] = h;
    }
  
  d[0] = e[0] = 0;
  for ( i = 0; i < n; i++ )
    {
      if ( fabs(d[i]) > 0 )
        {
          for ( j = 0; j < i; j++ )
            {
              g = 0;
              for ( k = 0; k < i; k++ )  g += a[i][k] * a[k][j];
              for ( k = 0; k < i; k++ )  a[k][j] -= g * a[k][i];
            }
        }
      d[i] = a[i][i];
      a[i][i] = 1;
      for ( j = 0; j < i; j++ ) a[j][i] = a[i][j] = 0;
    }
}


double pythagoras (double a, double b)
{
  double abs_a, abs_b, sqr;
  abs_a = fabs (a);
  abs_b = fabs (b);
  if (abs_a > abs_b)
    {
      sqr = abs_b / abs_a;
      sqr *= sqr;
      return abs_a * sqrt (1 + sqr);
    }
  else if (abs_b > abs_a)
    {
      sqr = abs_a / abs_b;
      sqr *= sqr;
      return abs_b * sqrt (1 + sqr);
    }
  else
    return M_SQRT2 * abs_a;
}

#define MAX_ITER 1000

void eigen_solution_of_triangular_matrix (double *d, double *e, int n,
					  double **a, const char *prompt)
{
  int i, k, l, m, iter;
  double b, c, f, g, p, r, s;
  static const double eps = 1e-6;
  
  for (i = 1; i < n; i++)
    e[i - 1] = e[i];
  
  e[n - 1] = 0;
  for (l = 0; l < n; l++)
    {
      iter = 0;
      while (1)
        {
          for (m = l; m < n - 1; m++)
            if (fabs (e[m]) <= eps * (fabs (d[m]) + fabs (d[m + 1])))
              break;
          if (m == l)
            {
	      // printf("found solution after %i Iteration\n", iter++);
              break;
            }
          iter++;
          if (iter == MAX_ITER)
            {
              fprintf (stderr, "%s: ERROR! Too many iterations while"
                       " determining the eigensolution!\n", prompt);
              exit (1);
            }
          g = (d[l + 1] - d[l]) / (2 * e[l]);
          r = pythagoras (g, 1);
          g = d[m] - d[l] + e[l] / (g + (fabs(g) > 0 ? (g >= 0 ? fabs (r) : -fabs (r)) : r));
          s = c = 1;
          p = 0;
          for (i = m - 1; i >= l; i--)
            {
              f = s * e[i];
              b = c * e[i];
              e[i + 1] = r = pythagoras (f, g);
              
              if ( DBL_IS_EQUAL(r, 0.) )
                {
                  d[i + 1] -= p;
                  e[m] = 0;
                  break;
                }
              
              s = f / r;
              c = g / r;
              g = d[i + 1] - p;
              r = (d[i] - g) * s + 2 * c * b;
              p = s * r;
              d[i + 1] = g + p;
              g = c * r - b;
              for (k = 0; k < n; k++)
                {
                  f = a[k][i + 1];
                  a[k][i + 1] = s * a[k][i] + c * f;
                  a[k][i] = c * a[k][i] - s * f;
                }
            }
          
          if ( DBL_IS_EQUAL(r, 0.) && i >= l ) continue;
          
          d[l] -= p;
          e[l] = g;
          e[m] = 0;
        }
      /*
       if (user_asked)
       {
       lock ();
       fprintf (stderr,
		   "%s: Status: Computing eigen solution pass 3 of 3"
		   " cycle %ld of %ld.\n", prompt, (long) (l + 1), (long) n);
       fflush (stderr);
       unlock ();
       user_asked = FALSE;
       }
       */
    }
}


int solution_of_linear_equation (double **a, double *b, int n)
{
  int *index;
  int sign;
  int not_singular;
  
  index = (int*) malloc(n * sizeof(int));
  
  not_singular = lu_decomposition (a, n, index, &sign);
  
  if (not_singular)
    lu_backsubstitution (a, n, index, b);
  
  free (index);
  
  return not_singular;
}


int inverse_of_matrix (double **a, double **b, int n)
{
  int *index;
  int sign;
  int i, j;
  int not_singular;
  double *col;
  
  index = (int*) malloc(n * sizeof(int));
  col = (double*) malloc(n * sizeof(double));
  
  not_singular = lu_decomposition (a, n, index, &sign);
  
  if (not_singular)
    {
      for (i = 0; i < n; i++)
        {
          for (j = 0; j < n; j++)
            col[j] = 0;
          col[i] = 1;
          lu_backsubstitution (a, n, index, col);
          for (j = 0; j < n; j++)
            b[j][i] = col[j];
        }
    }
  
  free (index);
  free (col);

  return not_singular;
}


int lu_decomposition (double **a, int n, int *index, int *sign)
{
  int i, imax = 0, j, k;
  double big, sum, temp;
  double *v;
  
  v = (double*) malloc(n * sizeof(double));
  *sign = 1;
  for (i = 0; i < n; i++)
    {
      big = 0;
      for (j = 0; j < n; j++)
        if ((temp = fabs (a[i][j])) > big)
          big = temp;
      
      if ( DBL_IS_EQUAL(big, 0.) ) return 0;
      
      v[i] = 1 / big;
    }
  for (j = 0; j < n; j++)
    {
      for (i = 0; i < j; i++)
        {
          sum = a[i][j];
          for (k = 0; k < i; k++)
            sum -= a[i][k] * a[k][j];
          a[i][j] = sum;
        }
      big = 0;
      for (i = j; i < n; i++)
        {
          sum = a[i][j];
          for (k = 0; k < j; k++)
            sum -= a[i][k] * a[k][j];
          a[i][j] = sum;
          if ((temp = v[i] * fabs (sum)) >= big)
            {
              big = temp;
              imax = i;
            }
        }
      if (j != imax)
        {
          for (k = 0; k < n; k++)
            {
              temp = a[imax][k];
              a[imax][k] = a[j][k];
              a[j][k] = temp;
            }
          *sign = -*sign;
          v[imax] = v[j];
        }
      index[j] = imax;
      
      if ( DBL_IS_EQUAL(a[j][j], 0.) ) return 0;
      
      if (j != n)
        {
          temp = 1 / a[j][j];
          for (i = j + 1; i < n; i++)
            a[i][j] *= temp;
        }
    }

  free (v);

  return 1;
}


void lu_backsubstitution (double **a, int n, int *index, double *b)
{
  int i, ii, ip, j;
  double sum;
  
  ii = 0;
  for (i = 0; i < n; i++)
    {
      ip = index[i];
      sum = b[ip];
      b[ip] = b[i];
      
      if (ii)
        for (j = ii; j < i; j++) sum -= a[i][j] * b[j];
      else if ( fabs(sum) > 0 )
        ii = i;
      
      b[i] = sum;
    }
  for (i = n - 1; i >= 0; i--)
    {
      sum = b[i];
      for (j = i + 1; j < n; j++)
        sum -= a[i][j] * b[j];
      b[i] = sum / a[i][i];
    }
}


void fft(double *real, double *imag, int n, int sign)
{
  /* n must be a power of 2 */
  /* sign should be 1 (FT) or -1 (reverse FT) */
  int i, j, j1, j2;
  int bit;
  int step;
  double temp_r, temp_i, norm;
  double w_r, w_i, ww_r, ww_i;
  
  /* Bit reversal part */
  for ( i = j = 0; i < n; i++ )	/* The bit pattern of i and j are reverse */
    {
      if ( i > j )
        {
          /* swap real part */
          temp_r = real[i];
          real[i] = real[j];
          real[j] = temp_r;

          /* swap imaginary part */
          temp_i = imag[i];
          imag[i] = imag[j];
          imag[j] = temp_i;
        }

      for ( bit = n >> 1; j & bit; bit >>= 1 ) j ^= bit;
      j |= bit;
    }

  /* Danielson-Lanczos Part */
  for (step = 1; step < n; step <<= 1)
    {
      w_r = cos (M_PI / step);
      w_i = sin (M_PI / step) * sign;
      ww_r = 1;
      ww_i = 0;
      for (i = 0; i < step; i++)
        {
          for (j1 = i, j2 = i + step; j2 < n; j1 += 2 * step, j2 += 2 * step)
            {
              temp_r = ww_r * real[j2] - ww_i * imag[j2];
              temp_i = ww_r * imag[j2] + ww_i * real[j2];
              real[j2] = real[j1] - temp_r;
              imag[j2] = imag[j1] - temp_i;
              real[j1] += temp_r;
              imag[j1] += temp_i;
            }
          temp_r = ww_r;
          ww_r = ww_r * w_r - ww_i * w_i;
          ww_i = temp_r * w_i + ww_i * w_r;
        }
    }

  norm = 1./sqrt(n);
  for ( i = 0; i < n; i++ )
    {
      real[i] *= norm;
      imag[i] *= norm;
    }
}


void ft(double *real, double *imag, int n, int sign)
{
  /* sign should be 1 (FT) or -1 (reverse FT) */
  int j, k;
  static double *work_r = 0, *work_i = 0;
  double sum_r, sum_i, norm;
  double w_r, w_i, ww_r, ww_i, temp_r;
  
  if (!work_r)
    {
      work_r = (double*) malloc(n * sizeof(double));
      /* free_at_exit (work_r); */
    }
  if (!work_i)
    {
      work_i = (double*) malloc(n * sizeof(double));
      /* free_at_exit (work_i); */
    }
  
  for (k = 0; k < n; k++)
    {
      w_r = cos (2 * M_PI * k / n);
      w_i = sin (2 * M_PI * k / n) * sign;
      ww_r = 1;
      ww_i = 0;
      sum_r = 0;
      sum_i = 0;
      for (j = 0; j < n; j++)
        {
          sum_r += real[j] * ww_r - imag[j] * ww_i;
          sum_i += real[j] * ww_i + imag[j] * ww_r;
          temp_r = ww_r;
          ww_r = ww_r * w_r - ww_i * w_i;
          ww_i = temp_r * w_i + ww_i * w_r;
        }
      work_r[k] = sum_r;
      work_i[k] = sum_i;
    }

  norm = 1./sqrt(n);
  for (k = 0; k < n; k++)
    {
      real[k] = work_r[k] * norm;
      imag[k] = work_i[k] * norm;
    }
}

/* reentrant version of ft */
void ft_r(double * restrict real, double * restrict imag, int n, int sign, double * restrict work_r, double * restrict work_i)
{
  /* sign should be 1 (FT) or -1 (reverse FT) */
  int j, k;
  double sum_r, sum_i, norm;
  double w_r, w_i, ww_r, ww_i, temp_r;
    
  for ( k = 0; k < n; k++ )
    {
      w_r = cos (2 * M_PI * k / n);
      w_i = sin (2 * M_PI * k / n) * sign;
      ww_r = 1;
      ww_i = 0;
      sum_r = 0;
      sum_i = 0;
      for ( j = 0; j < n; j++ )
        {
          sum_r += real[j] * ww_r - imag[j] * ww_i;
          sum_i += real[j] * ww_i + imag[j] * ww_r;
          temp_r = ww_r;
          ww_r = ww_r * w_r - ww_i * w_i;
          ww_i = temp_r * w_i + ww_i * w_r;
        }
      work_r[k] = sum_r;
      work_i[k] = sum_i;
    }

  norm = 1./sqrt(n);
  for ( k = 0; k < n; k++ )
    {
      real[k] = work_r[k] * norm;
      imag[k] = work_i[k] * norm;
    }
}


double lngamma (double x)
{
  double a, b, temp, ser;
  static double cof[6] =
  { 76.18009172947146, -86.50532032941677, 24.01409824083091,
    -1.231739572450155, 0.1208650973866179e-2, -0.5395239384953e-5
  };
  int j;
  
  a = b = x;
  temp = a + 5.5;
  temp -= (a + 0.5) * log (temp);
  ser = 1.000000000190015;
  for (j = 0; j <= 5; j++)
    ser += cof[j] / ++b;
  return -temp + log (2.5066282746310005 * ser / a);
}


double incomplete_gamma (double a, double x, const char *prompt)
{
  if (x < 0 || a <= 0)
    {
      fprintf (stderr,
               "%s: IMPLEMENTATION ERROR! (Invalid argument in function "
               "\"incomplete_gamma\")\n", prompt);
      exit (4);
    }
  if (x < (a + 1))
    return gamma_help_1 (a, x, prompt);
  else
    return 1 - gamma_help_2 (a, x, prompt);
}

static
double gamma_help_1 (double a, double x, const char *prompt)
{
  int i;
  double ap, del, gln, sum;
  static const double eps = 1e-20;
  
  gln = lngamma (a);
  ap = a;
  del = sum = 1 / a;
  for (i = 1; i <= 100; i++)
    {
      ap++;
      del *= x / ap;
      sum += del;
      if (fabs (del) < fabs (sum) * eps)
        return sum * exp (-x + a * log (x) - (gln));
    }
  fprintf (stderr,
           "%s: ERROR! Too many iterations in routine \"gamma_help_1\"!\n",
           prompt);
  exit (1);
  return 0;
}

static
double gamma_help_2 (double a, double x, const char *prompt)
{
  int i;
  double an, b, c, d, del, gln, h;
  double const very_small = 1000 * DBL_MIN;
  static const double eps = 1e-20;
  
  gln = lngamma (a);
  b = x + 1 - a;
  c = 1 / very_small;
  d = 1 / b;
  h = d;
  for (i = 1; i <= 100; i++)
    {
      an = -i * (i - a);
      b += 2;
      d = an * d + b;
      if (fabs (d) < very_small)
        d = very_small;
      c = b + an / c;
      if (fabs (c) < very_small)
        c = very_small;
      d = 1 / d;
      del = d * c;
      h *= del;
      if (fabs (del - 1) < eps)
        return exp (-x + a * log (x) - gln) * h;
    }
  fprintf (stderr,
           "%s: ERROR! Too many iterations in routine \"gamma_help_2\"!\n",
           prompt);
  exit (1);
  return -1;
}


double beta (double a, double b, const char *prompt)
{
  if (a <= 0 || b <= 0)
    {
      fprintf (stderr,
               "%s: IMPLEMENTATION ERROR! (Invalid argument in function "
               "\"beta\")\n", prompt);
      exit (4);
    }
  return exp (lngamma (a) + lngamma (b) - lngamma (a + b));
}


double incomplete_beta (double a, double b, double x, const char *prompt)
{
  double c;
  
  if (a <= 0 || b <= 0)
    {
      fprintf (stderr,
               "%s: IMPLEMENTATION ERROR! (Invalid argument in function "
               "\"incomplete_beta\")\n", prompt);
      exit (4);
    }
  
  if (x < 0 || x > 1)
    {
      fprintf (stderr, "%s: Value out of range (0-1)!\n", prompt);
      exit (4);
    }
  
  c = (DBL_IS_EQUAL(x, 0.) || DBL_IS_EQUAL(x, 1.)) ? 0 :
  exp (lngamma (a + b) - lngamma (a) - lngamma (b) + a * log (x) + b * log (1 - x));
  
  if (x < (a + 1) / (a + b + 2))
    return c * beta_help (a, b, x, prompt) / a;
  else
    return 1 - c * beta_help (b, a, 1 - x, prompt) / b;
}

static
double beta_help (double a, double b, double x, const char *prompt)
{
  int m, m2;
  double aa, c, d, del, h, qab, qam, qap;
  double const very_small = 1000 * DBL_MIN;
  static const double eps = 3e-07;
  
  qab = a + b;
  qap = a + 1;
  qam = a - 1;
  c = 1;
  d = 1 - qab * x / qap;
  if (fabs (d) < very_small)
    d = very_small;
  d = 1 / d;
  h = d;
  for (m = 1; m <= 100; m++)
    {
      m2 = 2 * m;
      aa = m * (b - m) * x / ((qam + m2) * (a + m2));
      d = 1 + aa * d;
      if (fabs (d) < very_small)
        d = very_small;
      c = 1 + aa / c;
      if (fabs (c) < very_small)
        c = very_small;
      d = 1 / d;
      h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      d = 1 + aa * d;
      if (fabs (d) < very_small)
        d = very_small;
      c = 1 + aa / c;
      if (fabs (c) < very_small)
        c = very_small;
      d = 1 / d;
      del = d * c;
      h *= del;
      if (fabs (del - 1.0) < eps)
        return h;
    }
  fprintf (stderr,
           "%s: ERROR! Too many iterations in routine \"beta_help\"!\n",
           prompt);
  exit (1);
  return -1;
}


double normal_density (double x)
{
  return M_2_SQRTPI / 2 / M_SQRT2 * exp (-x * x / 2);
}


double normal (double x, const char *prompt)
{
  return x > 0 ? 0.5 * (1 + incomplete_gamma (0.5, x * x / 2, prompt)) :
  x < 0 ? 0.5 * (1 - incomplete_gamma (0.5, x * x / 2, prompt)) : 0.5;
}


double normal_inv (double p, const char *prompt)
{
  static double last_p = 0.5, last_x = 0;
  double x, xx;
  static const double eps = 1e-10;
  
  if (p <= 0 || p >= 1)
    {
      fprintf (stderr,
               "%s: IMPLEMENTATION ERROR! (Invalid argument in function "
               "\"normal_inv\")\n", prompt);
      exit (4);
    }
  
  if ( DBL_IS_EQUAL(p, last_p) )
    return last_x;
  if ( DBL_IS_EQUAL(p, 0.5) )
    return 0;
  else if (p < 0.5)
    return -normal_inv (1 - p, prompt);
  else
    {
      x = 0;
      while (TRUE)
        {
          xx = x - (normal (x, prompt) - p) / normal_density (x);
          if (fabs (xx - x) < x * eps)
            break;
          x = xx;
        }
      last_p = p;
      last_x = x;
      return x;
    }
}


double student_t_density (double n, double x, const char *prompt)
{
  if (n <= 0)
    {
      fprintf (stderr,
               "%s: IMPLEMENTATION ERROR! (Invalid argument in function "
               "\"student_t_density\")\n", prompt);
      exit (4);
    }
  return exp (lngamma ((n + 1) / 2) - lngamma (n / 2)) / sqrt (n / 2)
  * pow ((1 + x * x / n), -(n + 1) / 2);
}


double student_t (double n, double x, const char *prompt)
{
  if (n <= 0)
    {
      fprintf (stderr,
               "%s: IMPLEMENTATION ERROR! (Invalid argument in function "
               "\"student_t\")\n", prompt);
      exit (4);
    }
  if (x > 0)
    return 1 - 0.5 * incomplete_beta (n / 2, 0.5, n / (n + x * x), prompt);
  else if (x < 0)
    return 0.5 * incomplete_beta (n / 2, 0.5, n / (n + x * x), prompt);
  else
    return 0.5;
}


double student_t_inv (double n, double p, const char *prompt)
{
  static double last_n = 1, last_p = 0.5, last_x = 0;
  double x, xx;
  static const double eps = 1e-10;
  
  if (n <= 0 || p <= 0 || p >= 1)
    {
      fprintf (stderr,
               "%s: IMPLEMENTATION ERROR! (Invalid argument in function "
               "\"student_t_inv\")\n", prompt);
      exit (4);
    }
  
  if ( DBL_IS_EQUAL(n, last_n) && DBL_IS_EQUAL(p, last_p) )
    return last_x;
  if ( DBL_IS_EQUAL(p, 0.5) )
    return 0;
  else if (p < 0.5)
    return -student_t_inv (n, 1 - p, prompt);
  else
    {
      x = 0;
      while (TRUE)
        {
          xx = x - (student_t (n, x, prompt) - p)
          / student_t_density (n, x, prompt);
          if (fabs (xx - x) < x * eps)
            break;
          x = xx;
        }
      last_n = n;
      last_p = p;
      last_x = x;
      return x;
    }
}


double chi_square_density (double n, double x, const char *prompt)
{
  if (n <= 0)
    {
      fprintf (stderr,
               "%s: IMPLEMENTATION ERROR! (Invalid argument in function "
               "\"chi_square_density\")\n", prompt);
      exit (4);
    }
  return x <= 0 ? 0 : pow (2, -n / 2) * 
  pow (x, n / 2 - 1) * exp (-x / 2 - lngamma (n / 2));
}


double chi_square (double n, double x, const char *prompt)
{
  if (n <= 0)
    {
      fprintf (stderr,
               "%s: IMPLEMENTATION ERROR! (Invalid argument in function "
               "\"chi_square\")\n", prompt);
      exit (4);
    }
  return x <= 0 ? 0 : incomplete_gamma (n / 2, x / 2, prompt);
}


double chi_square_inv (double n, double p, const char *prompt)
{
  static double last_n = -1, last_p = -1, last_x = -1;
  static double last_last_n = -1, last_last_p = -1, last_last_x = -1;
  double x, xx;
  static const double eps = 1e-10;
  
  if (n <= 0 || p <= 0 || p >= 1)
    {
      fprintf (stderr,
               "%s: IMPLEMENTATION ERROR! (Invalid argument in function "
               "\"chi_square_inv\")\n", prompt);
      exit (4);
    }
  
  if ( DBL_IS_EQUAL(n, last_n) && DBL_IS_EQUAL(p, last_p) )
    return last_x;
  
  if (DBL_IS_EQUAL(n, last_last_n) && DBL_IS_EQUAL(p, last_last_p) )
    return last_last_x;
  
  x = n;
  while (TRUE)
    {
      xx = x - (chi_square (n, x, prompt) - p)
      / chi_square_density (n, x, prompt);
      if (fabs (xx - x) < x * eps)
        break;
      if (xx < 0)
        x /= 2;
      else
        x = xx;
    }
  
  last_last_n = last_n;
  last_last_p = last_p;
  last_last_x = last_x;
  last_n = n;
  last_p = p;
  last_x = x;
  return x;
}


void chi_square_constants (double n, double p, double *c1, double *c2, const char *prompt)
{
  double delta_c1, delta_c2;
  static double last_n, last_p, last_c1, last_c2;
  double a11, a12, a21, a22, b1, b2;
  double det;
  static const double eps = 1e-10;
  
  if (n <= 0 || p <= 0 || p >= 1)
    {
      fprintf (stderr,
               "%s: IMPLEMENTATION ERROR! (Invalid argument in function "
               "\"chi_square_constants\")\n", prompt);
      exit (4);
    }
  
  if ( DBL_IS_EQUAL(n, last_n) && DBL_IS_EQUAL(p, last_p) )
    {
      *c1 = last_c1;
      *c2 = last_c2;
      return;
    }
  
  *c1 = n;
  *c2 = n;
  
  while (TRUE)
    {
      a11 = -chi_square_density (n, *c1, prompt);
      a12 = chi_square_density (n, *c2, prompt);
      a21 = -chi_square_density (n + 2, *c1, prompt);
      a22 = chi_square_density (n + 2, *c2, prompt);
      b1 = p + chi_square (n, *c1, prompt) - chi_square (n, *c2, prompt);
      b2 =
      p + chi_square (n + 2, *c1, prompt) - chi_square (n + 2, *c2, prompt);
      /* Solve ((a11,a12),(a21,a22))*(delta_c1,delta_c2)==(b1,b2) */
      det = a11 * a22 - a12 * a21;
      delta_c1 = (b1 * a22 - b2 * a12) / det;
      delta_c2 = (b2 * a11 - b1 * a21) / det;
      if (fabs (delta_c1) < *c1 * eps && fabs (delta_c2) < *c2 * eps)
        break;
      if (*c1 + delta_c1 >= n)
        *c1 = (n + *c1) / 2;
      else if (*c1 + delta_c1 <= 0)
        *c1 /= 2;
      else
        *c1 += delta_c1;
      if (*c2 + delta_c2 <= n)
        *c2 = (n + *c2) / 2;
      else
        *c2 += delta_c2;
    }
  
  last_n = n;
  last_p = p;
  last_c1 = *c1;
  last_c2 = *c2;
  return;
}


double beta_distr_density (double a, double b, double x, const char *prompt)
{
  if (a <= 0 || b <= 0)
    {
      fprintf (stderr,
               "%s: IMPLEMENTATION ERROR! (Invalid argument in function "
               "\"beta_distr_density\")\n", prompt);
      exit (4);
    }
  return x <= 0 ? 0 : x >= 1 ? 1 : pow (x, a - 1) 
  * pow (1 - x, b - 1) / beta (a, b, prompt);
}


double beta_distr (double a, double b, double x, const char *prompt)
{
  return incomplete_beta (a, b, x, prompt);
}


double beta_distr_inv (double a, double b, double p, const char *prompt)
{
  static double last_a = -1, last_b, last_p = -1, last_x = -1;
  static double last_last_a = -1, last_last_b = -1, last_last_p = -1,
  last_last_x = -1;
  double xx, x;
  static const double eps = 1e-10;
  
  if (a <= 0 || b <= 0 || p <= 0 || p >= 1)
    {
      fprintf (stderr,
               "%s: IMPLEMENTATION ERROR! (Invalid argument in function "
               "\"beta_distr_inv\")\n", prompt);
      exit (4);
    }
  
  if ( DBL_IS_EQUAL(a, last_a) && DBL_IS_EQUAL(b, last_b) && DBL_IS_EQUAL(p, last_p) )
    return last_x;
  
  if ( DBL_IS_EQUAL(a, last_last_a) && DBL_IS_EQUAL(b, last_last_b) && DBL_IS_EQUAL(p, last_last_p) )
    return last_last_x;
  
  x = a / (a + b);
  while (TRUE)
    {
      xx = x - (beta_distr (a, b, x, prompt) - p)
      / beta_distr_density (a, b, x, prompt);
      if (fabs (xx - x) < x * eps)
        break;
      if (xx <= 0)
        x /= 2;
      else if (xx >= 1)
        x = (1 + x) / 2;
      else
        x = xx;
    }
#if 0
  for (x_l = 0, x_r = 1; fabs (x_l - x_r) > eps;
       x = (x_l+x_r) / 2, beta_distr (a, b, x, prompt) < p ? (x_l=x):(x_r=x));
#endif
  
  last_last_a = last_a;
  last_last_b = last_b;
  last_last_p = last_p;
  last_last_x = last_x;
  last_a = a;
  last_b = b;
  last_p = p;
  last_x = x;
  return x;
}


void beta_distr_constants(double a, double b, double p, double *c1, double *c2,
                          const char *prompt)
{
  double delta_c1, delta_c2;
  static double last_a, last_b, last_p, last_c1, last_c2;
  double a11, a12, a21, a22, b1, b2;
  double det;
  static const double eps = 1e-10;
  
  if (a <= 0 || b <= 0 || p <= 0 || p >= 1)
    {
      fprintf (stderr,
               "%s: IMPLEMENTATION ERROR! (Invalid argument in function "
               "\"beta_distr_constants\")\n", prompt);
      exit (4);
    }
  
  if ( DBL_IS_EQUAL(a, last_a) && DBL_IS_EQUAL(b, last_b) && DBL_IS_EQUAL(p, last_p) )
    {
      *c1 = last_c1;
      *c2 = last_c2;
      return;
    }
  
#if 0
  *c1 = a / (a + b);
  *c2 = a / (a + b);
#endif
  *c1 = beta_distr_inv (a, b, p / 2, prompt);
  *c2 = beta_distr_inv (a, b, 1 - p / 2, prompt);
  
  while (TRUE)
    {
      a11 = -beta_distr_density (a, b, *c1, prompt);
      a12 =  beta_distr_density (a, b, *c2, prompt);
      a21 = -beta_distr_density (a + 1, b, *c1, prompt);
      a22 =  beta_distr_density (a + 1, b, *c2, prompt);
      b1 =
      p + beta_distr (a, b, *c1, prompt) - beta_distr (a, b, *c2, prompt);
      b2 =
      p + beta_distr (a + 1, b, *c1, prompt) - beta_distr (a + 1, b, *c2,
                                                           prompt);
      /* Solve ((a11,a12),(a21,a22))*(delta_c1,delta_c2)==(b1,b2) */
      det = a11 * a22 - a12 * a21;
      delta_c1 = (b1 * a22 - b2 * a12) / det;
      delta_c2 = (b2 * a11 - b1 * a21) / det;
      if (fabs (delta_c1) < *c1 * eps && fabs (delta_c2) < *c2 * eps)
        break;
      if (*c1 + delta_c1 >= a / (a + b))
        *c1 = (a / (a + b) + *c1) / 2;
      else if (*c1 + delta_c1 <= 0)
        *c1 /= 2;
      else
        *c1 += delta_c1;
      if (*c2 + delta_c2 >= 1)
        *c2 = (1 + *c2) / 2;
      else if (*c2 + delta_c2 <= a / (a + b))
        *c2 = (a / (a + b) + *c2) / 2;
      else
        *c2 += delta_c2;
    }
  
  last_a = a;
  last_b = b;
  last_p = p;
  last_c1 = *c1;
  last_c2 = *c2;
  return;
}


double fisher(double m, double n, double x, const char *prompt)
{
  if (m <= 0 || n <= 0)
    {
      fprintf (stderr,
               "%s: IMPLEMENTATION ERROR! (Invalid argument in function "
               "\"fisher\")\n", prompt);
      exit (4);
    }
  return incomplete_beta (m / 2, n / 2, n / (n + m * x), prompt);
}

/* ******************************************************************************** */
/*                                                                                  */
/*   P A R A L L E L   S O L U T I O N   O F   T H E   E I G E N   P R O B L E M    */
/*                     WITH ONE SIDED JACOBI ALGORITHM                              */
/*                                                                                  */
/* ******************************************************************************** */


void parallel_eigen_solution_of_symmetric_matrix(double **M, double *A, int n, const char func[])
{
  UNUSED(func);

  char *envstr;
  /* Get Environment variables if set */
  envstr = getenv("MAX_JACOBI_ITER");
  if ( envstr ) 
    max_jacobi_iter = atoi(envstr);
  else
    max_jacobi_iter = MAX_JACOBI_ITER;
  if ( cdoVerbose )
    cdoPrint("Using MAX_JACOBI_ITER %i from %s",
	     max_jacobi_iter, envstr?"Environment":"default");

   
  envstr = getenv("FNORM_PRECISION");
  if ( envstr ) 
    fnorm_precision = strtod(envstr,NULL);
  else
    fnorm_precision = FNORM_PRECISION;

  if ( cdoVerbose ) 
    cdoPrint("Using FNORM_PRECISION %g from %s",
	     fnorm_precision,envstr?"Environment":"default");

  // eigen_solution_of_symmetric_matrix(M, A, n, func);
  jacobi_1side(M, A, n);

  return;
}

/* ******************************************************************************** */
/* This routine rotates columns/rows i and j of a symmetric Matrix M in a fashion,  */
/* thus that the dot product of columns i and j 0 afterwards                        */
/*                                                                                  */
/* As this is done by a right-multiplication with a rotation matrix, which only     */
/* changes columns i and j, this can be carried out for n/2 pairs of columns at     */
/* the same time.                                                                   */
/* ******************************************************************************** */
void annihilate_1side(double **M, long i, long j, long n)
{

  double tk, ck, sk, alpha=0, beta=0, gamma=0, zeta=0;
  //  int first_annihilation = 0;
  long r;

  i--; j--;

  double *restrict mi = (double*) malloc(n*sizeof(double));
  double *restrict mj = (double*) malloc(n*sizeof(double));

  if ( ! mj || ! mi) 
    fprintf(stderr, "ERROR: allocation error - cannot allocate memory\n"
	            "ERROR: check stacksize and physically available memory\n");

  if ( j < i ) { int tmp = i; i = j; j = tmp; }
  
  double *restrict Mi = M[i];
  double *restrict Mj = M[j];

#if defined(HAVE_OPENMP4)
#pragma omp simd
#endif
  for ( r = 0; r < n; r++ )
    {
      alpha += Mj[r]*Mj[r];
      beta  += Mi[r]*Mi[r];
      gamma += Mi[r]*Mj[r];
    }

  // 2011-08-15 Cedrick Ansorge: bug fix
  //  tmp = fabs(gamma/sqrt(alpha/beta));
  double tmp = fabs(gamma/sqrt(alpha*beta));

  if ( tmp < fnorm_precision ) {
#if defined(_OPENMP)
#pragma omp critical 
#endif
    {
      n_finished++;
    }
    free(mi);
    free(mj);
    return;
  }
  
  zeta = (beta-alpha)/(2.*gamma);  // tan(2*theta)
  tk = 1./(fabs(zeta)+sqrt(1.+zeta*zeta)); 
  tk = zeta>0? tk : -tk;           // = cot(2*theta)
  ck = 1./sqrt(1.+tk*tk);          // = cos(theta)
  sk = ck*tk;                      // = sin(theta)
  
  // calculate a_i,j - tilde
  for ( r = 0; r < n; r++ )
    {
      mi[r] =  ck*Mi[r]  + sk*Mj[r];
      mj[r] = -sk*Mi[r]  + ck*Mj[r];
    }

  for ( r = 0; r < n; r++ ) Mi[r] = mi[r];
  for ( r = 0; r < n; r++ ) Mj[r] = mj[r];

  free(mi);
  free(mj);

  return;
}


int jacobi_1side(double **M, double *A, long n)
{  
  long i,j,k,m,r;
  long idx;
  long i_ann,j_ann;
  int n_iter = 0;
  int count=0;
  int **annihilations = NULL, *annihilations_buff = NULL;

  if ( n > 0 )
    {
      annihilations_buff = (int*) malloc(n*n*2*sizeof(int));
      annihilations = (int**) malloc((n*n)*sizeof(int*));
    }

  for(i=0;i<n*n;i++)
    annihilations[i] = & annihilations_buff[2*i];

  for( k=1; k<n+1; k++ ) {
    if ( k < n ) {
      for ( i=1;i<=(int)ceil(1./2.*(n-k));i++ ) {
	j = n-k+2-i;
        annihilations[count][0] = i;
        annihilations[count][1] = j;
	count++;
      }
      if ( k > 2 ) {
        for ( i=n-k+2;i<=n-(int)floor(1./2.*k);i++ ) {
          j = 2*n-k+2-i;
          annihilations[count][0] = i;
          annihilations[count][1] = j;
          count++;
        }
      }
    }
    else if ( k == n ) {
      for(i=2;i<=(int)ceil(1./2.*n);i++) {
        j = n+2-i;
        annihilations[count][0] = i;
        annihilations[count][1] = j;
        count++;
      }
    }
  }

  //  fprintf(stderr, "%i annihilations per sweep\n",count);

  n_finished = 0;

  //  override global openmp settings works
  //  omp_set_num_threads(2);

  while ( n_iter < max_jacobi_iter && n_finished < count ) {
    n_finished = 0;
    if ( n%2 == 1 ) {
      for ( m = 0; m < n; m++ ) {
#if defined(_OPENMP)
#pragma omp parallel for private(i,idx,i_ann,j_ann) shared(M,annihilations,n) reduction(+:n_finished)
#endif
        for ( i = 0; i < n/2; i++) {
          idx = m*(n/2)+i;
	  i_ann = annihilations[idx][0];
	  j_ann = annihilations[idx][1];
          if ( i_ann != j_ann && i_ann && j_ann )
	    annihilate_1side(M, i_ann, j_ann, n);
	}
      }
    }
    else { // n%2 == 0                                                                               
      for( m = 0; m < n; m++) {
#if defined(_OPENMP)
#pragma omp parallel for private(i,idx,i_ann,j_ann) shared(M,annihilations,n) reduction(+:n_finished)
#endif
        for( i = 0; i < n/2-(m%2); i++) {
	  idx = m/2 * ( n/2 + n/2-1);
          if ( m % 2 ) idx += n/2;
	  i_ann = annihilations[idx+i][0];
	  j_ann = annihilations[idx+i][1];
          if ( i_ann && j_ann && i_ann != j_ann ) 
	    annihilate_1side(M, i_ann, j_ann, n);
        }
      }
    }
    n_iter++;
  }

  if ( cdoVerbose ) 
    cdoPrint("Finished one-sided jacobi scheme for eigenvalue computation after %i iterations", n_iter);

  //  fprintf(stderr,"finished after %i sweeps (n_finished %i)\n",n_iter,n_finished);

  if ( n_iter == max_jacobi_iter && n_finished < count)
    {
      fprintf(stderr, 
	      "statistics-module (Warning): Eigenvalue computation with one-sided jacobi scheme\n"
	      "                             did not converge properly. %i of %i pairs of columns did\n"
	      "                             not achieve requested orthogonality of %10.6g\n",
	      count-n_finished,count, fnorm_precision);

      if ( n_finished == 0 ) {
	//	Do not overwrite results in case of insufficient convergence
	cdoWarning("Setting Matrix and Eigenvalues to 0 before return");
	for ( i=0; i<n; i++ ) {
	  memset(M[i],0,n*sizeof(double));
	  memset(A,0,n*sizeof(double));
	}
	return -1;
      }
    }
  // calculate  eigen values as sqrt(||m_i||)
  for ( i=0; i<n; i++ )
    {
      A[i] = 0;
      for ( r=0; r<n; r++ )
	A[i] += M[i][r] * M[i][r];
      A[i] = sqrt(A[i]);
      for ( r=0; r<n; r++ )
	M[i][r] /= A[i];
    }

  heap_sort(A,M,n);
  
  if ( annihilations      ) free(annihilations);
  if ( annihilations_buff ) free(annihilations_buff);
  
  return n_iter;
}

