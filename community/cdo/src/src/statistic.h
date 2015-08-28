#ifndef _STATISTIC_H

#define _STATISTIC_H

void eigen_solution_of_symmetric_matrix (double **a, double *eig_val, int n, const char *prompt);
int solution_of_linear_equation (double **a, double *b, int n);
int inverse_of_matrix (double **a, double **b, int n);
void fft(double *real, double *imag, int n, int sign);
//void ft (double *real, double *imag, int n, int sign);
void ft_r(double * restrict real, double * restrict imag, int n, int sign, double * restrict work_r, double * restrict work_i);
double lngamma (double x);
double beta (double a, double b, const char *prompt);
double incomplete_gamma (double a, double x, const char *prompt);
double incomplete_beta (double a, double b, double x, const char *prompt);
double normal_density (double x);
double normal (double x, const char *prompt);
double normal_inv (double p, const char *prompt);
double student_t_density (double n, double x, const char *prompt);
double student_t (double n, double x, const char *prompt);
double student_t_inv (double n, double p, const char *prompt);
double chi_square_density (double n, double x, const char *prompt);
double chi_square (double n, double x, const char *prompt);
double chi_square_inv (double n, double p, const char *prompt);
void chi_square_constants (double n, double p, double *c1, double *c2, const char *prompt);
double beta_distr_density (double a, double b, double x, const char *prompt);
double beta_distr (double a, double b, double x, const char *prompt);
double beta_distr_inv (double a, double b, double p, const char *prompt);
void beta_distr_constants (double a, double b, double p, double *c1, double *c2, const char *prompt);
double fisher (double m, double n, double x, const char *prompt);

// make heap sort accessible for eigen value computation in EOFs.c
void heap_sort (double *eig_val, double **a, int n);

// make parallel eigen solution accessible for eigen value computation in EOF3d.c
void parallel_eigen_solution_of_symmetric_matrix(double **M, double *A, int n, const char func[]);


#endif
