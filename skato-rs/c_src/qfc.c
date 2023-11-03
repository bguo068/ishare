/**************************************************************
 *
 *	qfc function from CompQuadForm package
 *		date: 12/05/2011
 *
 ****************************************************************/

#define UseDouble 0 /* all floating point double */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
/*#include <malloc.h>*/
#include <setjmp.h>

#define TRUE 1
#define FALSE 0
typedef int BOOL;

#define pi 3.14159265358979
#define log28 .0866 /*  log(2.0) / 8.0  */

// extern "C" {

static double sigma_sq, lambda_max, lambda_min, mean, q_rho;
static double intl, ersm;
static int count, n_lambdas, limits;
static BOOL ndtsrt, fail;
static int *df_arr, *th;
static double *lambda_arr, *non_centrality_arr;
static jmp_buf env;

static double exp1(double x) /* to avoid underflows  */
{
  return x < -50.0 ? 0.0 : exp(x);
}

static void counter(void)
/*  count number of calls to errbd, truncation, cfe */
{
  extern int count, limits;
  count = count + 1;
  if (count > limits)
    longjmp(env, 1);
}

static double square(double x) { return x * x; }

static double cube(double x) { return x * x * x; }

static double log1(double x, BOOL first)
/* if (first) log(1 + x) ; else  log(1 + x) - x */
{
  if (fabs(x) > 0.1) {
    return (first ? log(1.0 + x) : (log(1.0 + x) - x));
  } else {
    double s, s1, term, y, k;
    y = x / (2.0 + x);
    term = 2.0 * cube(y);
    k = 3.0;
    s = (first ? 2.0 : -x) * y;
    y = square(y);
    for (s1 = s + term / k; s1 != s; s1 = s + term / k) {
      k = k + 2.0;
      term = term * y;
      s = s1;
    }
    return s;
  }
}

static void order(void)
/* find order of absolute values of lb */
{
  int j, k;
  double lj;
  extern double *lambda_arr;
  extern int *th;
  extern int n_lambdas;
  extern BOOL ndtsrt;
  for (j = 0; j < n_lambdas; j++) {
    lj = fabs(lambda_arr[j]);
    for (k = j - 1; k >= 0; k--) {
      if (lj > fabs(lambda_arr[th[k]]))
        th[k + 1] = th[k];
      else
        goto l1;
    }
    k = -1;
  l1:
    th[k + 1] = j;
  }
  ndtsrt = FALSE;
}

static double errbd(double u, double *cx)
/*  find bound on tail probability using mgf, cutoff
   point returned to *cx */
{
  double sum1, lj, ncj, x, y, xconst;
  int j, nj;
  extern double sigma_sq, *lambda_arr, *non_centrality_arr;
  extern int *df_arr;
  extern int n_lambdas;
  counter();
  xconst = u * sigma_sq;
  sum1 = u * xconst;
  u = 2.0 * u;
  for (j = n_lambdas - 1; j >= 0; j--) {
    nj = df_arr[j];
    lj = lambda_arr[j];
    ncj = non_centrality_arr[j];
    x = u * lj;
    y = 1.0 - x;
    xconst = xconst + lj * (ncj / y + nj) / y;
    sum1 = sum1 + ncj * square(x / y) + nj * (square(x) / y + log1(-x, FALSE));
  }
  *cx = xconst;
  return exp1(-0.5 * sum1);
}

static double ctff(double accx, double *upn)
/*  find ctff so that p(qf > ctff) < accx  if (upn > 0,
    p(qf < ctff) < accx otherwise */
{
  double u1, u2, u, rb, xconst, c1, c2;
  extern double lambda_min, lambda_max, mean;
  u2 = *upn;
  u1 = 0.0;
  c1 = mean;
  rb = 2.0 * ((u2 > 0.0) ? lambda_max : lambda_min);
  for (u = u2 / (1.0 + u2 * rb); errbd(u, &c2) > accx;
       u = u2 / (1.0 + u2 * rb)) {
    u1 = u2;
    c1 = c2;
    u2 = 2.0 * u2;
  }
  for (u = (c1 - mean) / (c2 - mean); u < 0.9; u = (c1 - mean) / (c2 - mean)) {
    u = (u1 + u2) / 2.0;
    if (errbd(u / (1.0 + u * rb), &xconst) > accx) {
      u1 = u;
      c1 = xconst;
    } else {
      u2 = u;
      c2 = xconst;
    }
  }
  *upn = u2;
  return c2;
}

static double truncation(double u, double tausq)
/* bound integration error due to truncation at u */
{
  double sum1, sum2, prod1, prod2, prod3, lj, ncj, x, y, err1, err2;
  int j, nj, s;
  extern double sigma_sq, *lambda_arr, *non_centrality_arr;
  extern int *df_arr;
  extern int n_lambdas;

  counter();
  sum1 = 0.0;
  prod2 = 0.0;
  prod3 = 0.0;
  s = 0;
  sum2 = (sigma_sq + tausq) * square(u);
  prod1 = 2.0 * sum2;
  u = 2.0 * u;
  for (j = 0; j < n_lambdas; j++) {
    lj = lambda_arr[j];
    ncj = non_centrality_arr[j];
    nj = df_arr[j];
    x = square(u * lj);
    sum1 = sum1 + ncj * x / (1.0 + x);
    if (x > 1.0) {
      prod2 = prod2 + nj * log(x);
      prod3 = prod3 + nj * log1(x, TRUE);
      s = s + nj;
    } else
      prod1 = prod1 + nj * log1(x, TRUE);
  }
  sum1 = 0.5 * sum1;
  prod2 = prod1 + prod2;
  prod3 = prod1 + prod3;
  x = exp1(-sum1 - 0.25 * prod2) / pi;
  y = exp1(-sum1 - 0.25 * prod3) / pi;
  err1 = (s == 0) ? 1.0 : x * 2.0 / s;
  err2 = (prod3 > 1.0) ? 2.5 * y : 1.0;
  if (err2 < err1)
    err1 = err2;
  x = 0.5 * sum2;
  err2 = (x <= y) ? 1.0 : y / x;
  return (err1 < err2) ? err1 : err2;
}

static void findu(double *utx, double accx)
/*  find u such that truncation(u) < accx and truncation(u / 1.2) > accx */
{
  double u, ut;
  int i;
  static double divis[] = {2.0, 1.4, 1.2, 1.1};
  ut = *utx;
  u = ut / 4.0;
  if (truncation(u, 0.0) > accx) {
    for (u = ut; truncation(u, 0.0) > accx; u = ut)
      ut = ut * 4.0;
  } else {
    ut = u;
    for (u = u / 4.0; truncation(u, 0.0) <= accx; u = u / 4.0)
      ut = u;
  }
  for (i = 0; i < 4; i++) {
    u = ut / divis[i];
    if (truncation(u, 0.0) <= accx)
      ut = u;
  }
  *utx = ut;
}

static void integrate(int nterm, double interv, double tausq, BOOL mainx)
/*  carry out integration with nterm terms, at stepsize
   interv.  if (! mainx) multiply integrand by
      1.0-exp(-0.5*tausq*u^2) */
{
  double inpi, u, sum1, sum2, sum3, x, y, z;
  int k, j, nj;
  extern double intl, ersm;
  extern double sigma_sq, q_rho;
  extern int *df_arr;
  extern double *lambda_arr, *non_centrality_arr;
  extern int n_lambdas;
  inpi = interv / pi;
  for (k = nterm; k >= 0; k--) {
    u = (k + 0.5) * interv;
    sum1 = -2.0 * u * q_rho;
    sum2 = fabs(sum1);
    sum3 = -0.5 * sigma_sq * square(u);
    for (j = n_lambdas - 1; j >= 0; j--) {
      nj = df_arr[j];
      x = 2.0 * lambda_arr[j] * u;
      y = square(x);
      sum3 = sum3 - 0.25 * nj * log1(y, TRUE);
      y = non_centrality_arr[j] * x / (1.0 + y);
      z = nj * atan(x) + y;
      sum1 = sum1 + z;
      sum2 = sum2 + fabs(z);
      sum3 = sum3 - 0.5 * x * y;
    }
    x = inpi * exp1(sum3) / u;
    if (!mainx)
      x = x * (1.0 - exp1(-0.5 * tausq * square(u)));
    sum1 = sin(0.5 * sum1) * x;
    sum2 = 0.5 * sum2 * x;
    intl = intl + sum1;
    ersm = ersm + sum2;
  }
}

static double cfe(double x)
/*  coef of tausq in error when convergence factor of
   exp1(-0.5*tausq*u^2) is used when df is evaluated at x */
{
  double axl, axl1, axl2, sxl, sum1, lj;
  int j, k, t;
  extern BOOL ndtsrt, fail;
  extern int *th, *df_arr;
  extern double *lambda_arr, *non_centrality_arr;
  extern int n_lambdas;
  counter();
  if (ndtsrt)
    order();
  axl = fabs(x);
  sxl = (x > 0.0) ? 1.0 : -1.0;
  sum1 = 0.0;
  for (j = n_lambdas - 1; j >= 0; j--) {
    t = th[j];
    if (lambda_arr[t] * sxl > 0.0) {
      lj = fabs(lambda_arr[t]);
      axl1 = axl - lj * (df_arr[t] + non_centrality_arr[t]);
      axl2 = lj / log28;
      if (axl1 > axl2)
        axl = axl1;
      else {
        if (axl > axl2)
          axl = axl2;
        sum1 = (axl - axl1) / lj;
        for (k = j - 1; k >= 0; k--)
          sum1 = sum1 + (df_arr[th[k]] + non_centrality_arr[th[k]]);
        goto l;
      }
    }
  }
l:
  if (sum1 > 100.0) {
    fail = TRUE;
    return 1.0;
  } else
    return pow(2.0, (sum1 / 4.0)) / (pi * square(axl));
}

void qfc_1(double *lb1, double *nc1, int *n1, int *r1, double *sigma,
           double *c1, int *lim1, double *acc, double *trace_arr,
           int *idx_fault, double *res)

/*  distribution function of a linear combination of non-central
   chi-squared random variables :

input:
   lb[j]            coefficient of j-th chi-squared variable
   nc[j]            non-centrality parameter
   n[j]             degrees of freedom
   j = 0, 2 ... r-1
   sigma            coefficient of standard normal variable
   c                point at which df is to be evaluated
   lim              maximum number of terms in integration
   acc              maximum error

output:
   ifault = 1       required accuracy NOT achieved
            2       round-off error possibly significant
            3       invalid parameters
            4       unable to locate integration parameters
            5       out of memory

   trace[0]         absolute sum
   trace[1]         total number of integration terms
   trace[2]         number of integrations
   trace[3]         integration interval in final integration
   trace[4]         truncation point in initial integration
   trace[5]         s.d. of initial convergence factor
   trace[6]         cycles to locate integration parameters     */

{
  int j, df_j, nt, ntm;
  double acc1, abs_lambda_max, xlim, xnt, xntm;
  double utx, tausq, sd, intv, intv1, x, up, un, d1, d2, lambda_j, noncentral_j;
  extern double sigma_sq, lambda_max, lambda_min, mean;
  extern double intl, ersm;
  extern int n_lambdas, limits;
  extern double q_rho;
  extern int *df_arr, *th;
  extern double *lambda_arr, *non_centrality_arr;
  double qfval = -1.0;
  static int rats[] = {1, 2, 4, 8};

  if (setjmp(env) != 0) {
    *idx_fault = 4;
    goto endofproc;
  }
  n_lambdas = r1[0];
  limits = lim1[0];
  q_rho = c1[0];
  df_arr = n1;
  lambda_arr = lb1;
  non_centrality_arr = nc1;
  for (j = 0; j < 7; j++)
    trace_arr[j] = 0.0;
  *idx_fault = 0;
  count = 0;
  intl = 0.0;
  ersm = 0.0;
  qfval = -1.0;
  acc1 = acc[0];
  ndtsrt = TRUE;
  fail = FALSE;
  xlim = (double)limits;
  th = (int *)malloc(n_lambdas * (sizeof(int)));
  if (!th) {
    *idx_fault = 5;
    goto endofproc;
  }

  /* find mean, sd, max and min of lb,
     check that parameter values are valid */
  sigma_sq = square(sigma[0]);
  sd = sigma_sq;
  lambda_max = 0.0;
  lambda_min = 0.0;
  mean = 0.0;
  for (j = 0; j < n_lambdas; j++) {
    df_j = df_arr[j];
    lambda_j = lambda_arr[j];
    noncentral_j = non_centrality_arr[j];
    if (df_j < 0 || noncentral_j < 0.0) {
      *idx_fault = 3;
      goto endofproc;
    }
    sd = sd + square(lambda_j) * (2 * df_j + 4.0 * noncentral_j);
    mean = mean + lambda_j * (df_j + noncentral_j);
    if (lambda_max < lambda_j)
      lambda_max = lambda_j;
    else if (lambda_min > lambda_j)
      lambda_min = lambda_j;
  }
  if (sd == 0.0) {
    qfval = (q_rho > 0.0) ? 1.0 : 0.0;
    goto endofproc;
  }
  if (lambda_min == 0.0 && lambda_max == 0.0 && sigma[0] == 0.0) {
    *idx_fault = 3;
    goto endofproc;
  }
  sd = sqrt(sd);
  abs_lambda_max = (lambda_max < -lambda_min) ? -lambda_min : lambda_max;

  /* starting values for findu, ctff */
  utx = 16.0 / sd;
  up = 4.5 / sd;
  un = -up;
  /* truncation point with no convergence factor */
  findu(&utx, .5 * acc1);
  /* does convergence factor help */
  if (q_rho != 0.0 && (abs_lambda_max > 0.07 * sd)) {
    tausq = .25 * acc1 / cfe(q_rho);
    if (fail)
      fail = FALSE;
    else if (truncation(utx, tausq) < .2 * acc1) {
      sigma_sq = sigma_sq + tausq;
      findu(&utx, .25 * acc1);
      trace_arr[5] = sqrt(tausq);
    }
  }
  trace_arr[4] = utx;
  acc1 = 0.5 * acc1;

  /* find RANGE of distribution, quit if outside this */
l1:
  d1 = ctff(acc1, &up) - q_rho;
  if (d1 < 0.0) {
    qfval = 1.0;
    goto endofproc;
  }
  d2 = q_rho - ctff(acc1, &un);
  if (d2 < 0.0) {
    qfval = 0.0;
    goto endofproc;
  }
  /* find integration interval */
  intv = 2.0 * pi / ((d1 > d2) ? d1 : d2);
  /* calculate number of terms required for main and
     auxillary integrations */
  xnt = utx / intv;
  xntm = 3.0 / sqrt(acc1);
  if (xnt > xntm * 1.5) {
    /* parameters for auxillary integration */
    if (xntm > xlim) {
      *idx_fault = 1;
      goto endofproc;
    }
    ntm = (int)floor(xntm + 0.5);
    intv1 = utx / ntm;
    x = 2.0 * pi / intv1;
    if (x <= fabs(q_rho))
      goto l2;
    /* calculate convergence factor */
    tausq = .33 * acc1 / (1.1 * (cfe(q_rho - x) + cfe(q_rho + x)));
    if (fail)
      goto l2;
    acc1 = .67 * acc1;
    /* auxillary integration */
    integrate(ntm, intv1, tausq, FALSE);
    xlim = xlim - xntm;
    sigma_sq = sigma_sq + tausq;
    trace_arr[2] = trace_arr[2] + 1;
    trace_arr[1] = trace_arr[1] + ntm + 1;
    /* find truncation point with new convergence factor */
    findu(&utx, .25 * acc1);
    acc1 = 0.75 * acc1;
    goto l1;
  }

  /* main integration */
l2:
  trace_arr[3] = intv;
  if (xnt > xlim) {
    *idx_fault = 1;
    goto endofproc;
  }
  nt = (int)floor(xnt + 0.5);
  integrate(nt, intv, 0.0, TRUE);
  trace_arr[2] = trace_arr[2] + 1;
  trace_arr[1] = trace_arr[1] + nt + 1;
  qfval = 0.5 - intl;
  trace_arr[0] = ersm;

  /* test whether round-off error could be significant
     allow for radix 8 or 16 machines */
  up = ersm;
  x = up + acc[0] / 10.0;
  for (j = 0; j < 4; j++) {
    if (rats[j] * x == rats[j] * up)
      *idx_fault = 2;
  }

endofproc:
  free((char *)th);
  trace_arr[6] = (double)count;
  res[0] = qfval;
  return;
}

//}
