/*
 * quadprog solver.
 * Copyright (C) shouran.ma@rwth-aachen.de
 *
 * This file is part of HECTR.
 *
 * HECTR is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * HECTR is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this program; if not, see <http://www.gnu.org/licenses/>.
 */

#include "hectr.h"
#include <stdio.h>
#include <float.h>

BEGIN_DECLS

/**
 * qp_lagrange - Lagrangian multiplier method to solve quadratic programming
 * with linear equality constraints.
 *
 * @w[n*1] the state
 * @lambda[m*1]: Lagrangian multiplier
 * @n: dimension of the state / the degree of freedom
 * @m: number of equality constraint
 * @H[n*n]: Hessian matrix, positive and symmetric.
 * @c[n*1] 
 * @A[m*n] Aw+b=0, coefficient matrix for linear equality constraint.
 * @b[m*1] Aw+b=0, coefficient matrix for linear equality constraint.
 */
static void qp_lagrange(double w[], double lambda[],
  const unsigned int n, const unsigned int m,
  const double H[], const double c[],
  const double A[], const double b[])
{
  /**
   * min 0.5 w'*H*w + c'*w
   *  w
   * s.t. A*w+b = 0
   *
   * L = 0.5 w'*H*w + c'*w + lambda*(A*w+b)
   *
   * [H, A.T] * [  w   ] = [-c]  =>  [  w   ] = [M11, M12] * [-c]
   * [A, 0  ]   [lambda]   [-b]      [lambda]   [M21, M22]   [-b]
   *
   * M22 = -(A * H.I * A.T).I
   * M21 = (A * H.I * A.T).I * A*H.I 
   *     = -M22*A*H.I
   * M12 = M21.T
   * M11 = H.I - H.I*A.T * (A * H.I * A.T).I * A*H.I
   *     = H.I + H.I*A.T * M22 * A*H.I
   *     = H.I - H.I*A.T * M21
   */

  int err = 0;

  double H_I[n*n];
  err |= dgeinv(H_I, n, H);

  if (!m) {
    /* w = M11*(-c) + M12*(-b) */
    memset(w, 0, n*sizeof(double));
    for (unsigned int i=0; i<n; i++)
      for (unsigned int j=0; j<n; j++)
        w[i] += H_I[i*n+j]*(-c[j]);
    /* lambda = M21*(-c) + M22*(-b) */
    memset(lambda, 0, m*sizeof(double));
    return ;
  }

  /* A * H.I [m*n] */
  double A_mul_H_I[m*n];
  memset(A_mul_H_I, 0, sizeof(A_mul_H_I));
  for (unsigned int i=0, all=m*n; i<all; i++)
    for (unsigned int j=0; j<n; j++)
      A_mul_H_I[i] += A[(i/n)*n+j]*H_I[j*n+(i%n)];

  /* H.I[n*n] * A.T[n*m] */
  double H_I_mul_A_T[n*m];
  memset(H_I_mul_A_T, 0, sizeof(H_I_mul_A_T));
  for (unsigned int i=0, all=n*m; i<all; i++)
    for (unsigned int j=0; j<n; j++)
      H_I_mul_A_T[i] += H_I[(i/m)*n+j]*A[(i%m)*n+j];

  /* M22 = -(A * H.I * A.T).I */
  double M22[m*m];
  memset(M22, 0, sizeof(M22));
  for (unsigned int i=0, all=m*m; i<all; i++)
    for (unsigned int j=0; j<n; j++)
      M22[i] += A[(i/m)*n+j] * H_I_mul_A_T[j*m+(i%m)];
  err |= dgeinv(M22, m, M22);
  for (unsigned int i=0, all=m*m; i<all; i++)
    M22[i] = -M22[i];

  /* M21[m*n] = -M22[m*m] * (A*H.I)[m*n] */
  double M21[m*n];
  memset(M21, 0, sizeof(M21));
  for (unsigned int i=0, all=m*n; i<all; i++)
    for (unsigned int j=0; j<m; j++)
      M21[i] += -M22[(i/n)*m+j] * A_mul_H_I[j*n+(i%n)];

  double M12[n*m];
  memset(M12, 0, sizeof(M12));
  for (unsigned int i=0, all=n*m; i<all; i++)
    M12[(i/m)*m+(i%m)] = M21[(i%m)*n+(i/m)];

  /* M11 = H.I - M12 * A * H.I */
  double M11[n*n];
  memset(M11, 0, sizeof(M11));
  for (unsigned int i=0, all=n*n; i<all; i++) {
    for (unsigned int j=0; j<m; j++)
      M11[i] += -H_I_mul_A_T[(i/n)*m+j] * M21[j*n+(i%n)];
    M11[i] += H_I[i];
  }

  double KKTinv[(n+m)*(n+m)];
  memset(KKTinv, 0, sizeof(KKTinv));
  for (unsigned int i=0; i<n+m; i++) {
    for (unsigned int j=0; j<n+m; j++) {
      if ((i<n)&&(j<n))
        KKTinv[i*(n+m)+j] = M11[i*n+j];
      if ((i<n)&&(n<=j))
        KKTinv[i*(n+m)+j] = M12[i*m+j-n];
      if ((n<=i)&&(j<n))
        KKTinv[i*(n+m)+j] = M21[(i-n)*n+j];
      if ((n<=i)&&(n<=j))
        KKTinv[i*(n+m)+j] = M22[(i-n)*m+j-n];
    }
  }

  double vec[n+m];
  for (unsigned int i=0; i<n+m; i++)
    vec[i] = (i<n)? -c[i] : -b[i-n];

  double ret[n+m];
  memset(ret, 0, sizeof(ret));
  for (unsigned int i=0; i<n+m; i++)
    for (unsigned int j=0; j<n+m; j++)
      ret[i] += KKTinv[i*(n+m)+j]*vec[j];
  memcpy(w, ret, n*sizeof(double));
  memcpy(lambda, ret+n, m*sizeof(double));
}

static void check_w0(double w[],
  const unsigned int n, const unsigned int meq,
  const double Aeq[], const double beq[])
{
  unsigned char should_reset_w = 0;
  double Aw_add_b[meq];
  memset(Aw_add_b, 0, sizeof(Aw_add_b));
  for (unsigned int i=0; i<meq; i++) {
    for (unsigned int j=0; j<n; j++)
      Aw_add_b[i] += Aeq[i*n+j]*w[j];
    Aw_add_b[i] += beq[i];
    if (fabs(Aw_add_b[i]) > FLT_EPSILON) {
      should_reset_w = 1;
      break;
    }
  }
  if (should_reset_w) {
    double Aeq_pinv[n*meq];
    dgepinv(Aeq_pinv, meq, n, Aeq);
    memset(w, 0, n*sizeof(double));
    for (unsigned int i=0; i<n; i++)
      for (unsigned int j=0; j<meq; j++)
        w[i] += Aeq_pinv[i*meq+j]*(-beq[j]);
  }
}

static void set_Ab(double A[], double b[],
  const unsigned int m, const unsigned int n, const unsigned int mineq,
  const double Ain[], const double bin[],
  const double lb[], const double ub[])
{
  memset(A, 0, m*n*sizeof(double));
  memset(b, 0, m*sizeof(double));
  memcpy(A, Ain, mineq*n*sizeof(double));
  memcpy(b, bin, mineq*sizeof(double));
  if (lb) {
    double nI[n*n];
    memset(nI, 0, sizeof(nI));
    for (unsigned int i=0; i<n; i++)
      nI[i*n+i]=-1;
    memcpy(&A[mineq*n], nI, sizeof(nI));
    memcpy(&b[mineq], lb, n*sizeof(double));
  }
  if (ub) {
    double pI[n*n];
    memset(pI, 0, sizeof(pI));
    for (unsigned int i=0; i<n; i++)
      pI[i*n+i]=1;
    if (lb) {
      memcpy(&A[(mineq+n)*n], pI, sizeof(pI));
      for (unsigned int i=0; i<n; i++)
        b[mineq+n+i] = -ub[i];
    }
    else {
      memcpy(&A[mineq*n], pI, sizeof(pI));
      for (unsigned int i=0; i<n; i++)
        b[mineq+i] = -ub[i];
    }
  }
}

static void init_W(unsigned int W[],
  const unsigned int n, const unsigned int m,
  const double w[], const double A[], const double b[])
{
  double Aw_add_b[m];
  memset(Aw_add_b, 0, sizeof(Aw_add_b));
  for (unsigned int i=0; i<m; i++) {
    for (unsigned int j=0; j<n; j++)
      Aw_add_b[i] += A[i*n+j]*w[j];
    Aw_add_b[i] += b[i];
  }
  for (unsigned int i=0; i<m; i++) {
    if (-FLT_EPSILON < Aw_add_b[i])
      W[i] = 1;
    else
      W[i] = 0;
  }
}

/**
 * @p[n]: the increment of the state.
 * @lambda[m]: the multiplier of effective equality constraints.
 * @n: dimension of the state / the degree of freedom.
 * @m: the number of effective equality constraints.
 * @w[n]: the current state.
 */
static void calc_p(double p[], double lambda[],
  const unsigned int n, const unsigned int m, const unsigned int meq,
  const unsigned int W[],
  const double w[],
  const double H[], const double c[],
  const double A[], const double b[],
  const double Aeq[], const double beq[])
{
  unsigned int m_act = 0;
  for (unsigned int i=0; i<m; i++)
    if (W[i])
      m_act++;
  unsigned int meq_eff = meq + m_act;

  /**
   * k - the index to track the union of eqlin set and ineqlin set
   *
   * Example: 2 eqlin constraints and 5 ineqlin constraints: meq=2, m=4.
   *
   *   eqlin: [0 1]
   * ineqlin:     [0 1 2 3 4]
   *               i=0
   *     all: [0 1 2 3 4 5 6]
   *               k=2
   *
   * "i,j" is usually for local addressing: i=0,...,m/meq-1, j=0,...,n-1.
   * "k" for "global" addressing: j=0,...,meq-1,meq,...,meq+m-1
   */
  unsigned int k;

  /* set effective eqlin constraints (eqlin + active ineqlin constraints) */
  double Aeq_eff[meq_eff*n];
  double beq_eff[meq_eff*1];
  memset(Aeq_eff, 0, sizeof(Aeq_eff));
  memset(beq_eff, 0, sizeof(beq_eff));
  if (meq) {
    memcpy(Aeq_eff, Aeq, meq*n*sizeof(double));
    memcpy(beq_eff, beq, meq*1*sizeof(double));
  }
  /**
   * ineqlin:     [0 1 2 3 4]
   *       W: [0 1 N N Y N Y]
   *  (initialize) k=2 |   |
   *                   k=2 |
   *                       k=3
   */
  k = meq;
  for (unsigned int i=0; (i<m); i++) {
    if (W[i]) {
      memcpy(&Aeq_eff[k*n], &A[i*n], n*sizeof(double));
      memcpy(&beq_eff[k*1], &b[i*1],   sizeof(double));
      k++;
    }
  }

  /* g = Hw+c */
  double g[n];
  memset(g, 0, sizeof(g));
  for (unsigned int i=0; i<n; i++) {
    for (unsigned int j=0; j<n; j++)
      g[i] += H[i*n+j]*w[j];
    g[i] += c[i];
  }

  /* h = Aeq_eff*w + beq_eff */
  double h[meq_eff];
  memset(h, 0, sizeof(h));
  for (unsigned int i=0; i<meq_eff; i++) {
    for (unsigned int j=0; j<n; j++)
      h[i] += Aeq_eff[i*n+j]*w[j];
    h[i] += beq_eff[i];
  }

  /* lambda_eff = lambda + mu */
  double lambda_eff[meq_eff];
  memset(lambda_eff, 0, sizeof(lambda_eff));
  qp_lagrange(p, lambda_eff, n, meq_eff, H, g, Aeq_eff, h);

  /* store lambda */
  if (meq)
    memcpy(lambda, lambda_eff, meq*sizeof(double));
  /**
   * lambda_eff: [5     -4                   -2            -1]
   *                                         k=0           k=1
   * lambda:     [5     -4      0      0     -2      0     -1]
   *                         i+meq=2                     i+meq=6
   */
  k = 0;
  for (unsigned int i=0; i<m; i++) {
    if (W[i]) {
      lambda[meq+i] = lambda_eff[meq+k];
      k++;
    }
    else
      lambda[i+meq] = 0;
  }
}

/**
 * calc_mu - calculate the multipliers of the ineqlin constraints
 */
static void calc_mu(unsigned int *mu_idx_, double *mu_,
  const unsigned int m, const unsigned int meq,
  const unsigned int W[], const double lambda[])
{
  unsigned int mu_idx = UINT32_MAX;
  double mu = DBL_MAX;
  for (unsigned int i=0; i<m; i++) {
    if ((W[i]!=0) && (lambda[meq+i] < mu)) {
      mu_idx = i;
      mu = lambda[meq+mu_idx];
    }
  }
  *mu_idx_ = mu_idx;
  *mu_ = mu;
}

/**
 * calc_alpha - calculate alpha
 * 
 * alpha = min{1, min{-(aw+b)/ap}}, where ap>0
 */
static void calc_alpha(unsigned int *alpha_idx_, double *alpha_,
  const unsigned int n, const unsigned int m,
  const unsigned int W[], const double w[],
  const double p[],
  const double A[], const double b[])
{
  unsigned int alpha_idx = UINT32_MAX;
  double alpha = DBL_MAX;
  for (unsigned int i=0; i<m; i++) {
    if (!W[i]) { /* i not in W, i.e. W[i]==0 */
      double a[n];
      memcpy(a, &A[i*n], n*sizeof(double));
      double ap=0;
      for (unsigned int j=0; j<n; j++)
        ap += a[j]*p[j];
      if (ap > 0) {
        double aw_add_b=0;
        for (unsigned int j=0; j<n; j++)
          aw_add_b += a[j]*w[j];
        aw_add_b += b[i];
        double temp = -aw_add_b/ap;
        if (temp < alpha) {
          alpha_idx = i;
          alpha = temp;
        }
      }
    }
  }
  if (alpha >= 1) {
    alpha = 1;
    alpha_idx = m;
    /* m is out of range of W[], just means that need not add constraints. */
  }
  *alpha_idx_ = alpha_idx;
  *alpha_ = alpha;
}

/**
 * quadprog - quadprog solver with active-set method.
 *
 * Solve the quadratic program 
 * min 0.5 x'*H*x + c'*x
 *  x
 * subject to
 * A*x + b <= 0,
 * Aeq*x + beq = 0,
 * lb <= x <= ub.
 * 
 * @w[n]: the state that minimize the objective function.
 * @lambda[meq+m]: multipliers of equality+inequality constraints.
 * @n: dimension of w
 * @m: number of inequality constraints
 * @meq: number of equality constraints
 * @H[n*n]: Hessian matrix
 * @c[n]:
 * @A[m*n]: coefficient matrix of inequality constraints
 * @b[m]: 
 * @Aeq[meq*n]: coefficient matrix of equality constraints
 * @beq[meq]:
 */
int quadprog(double w[],
  const unsigned int n, const unsigned int mineq, const unsigned int meq,
  const double H[], const double c[],
  const double Ain[], const double bin[],
  const double Aeq[], const double beq[],
  const double lb[], const double ub[])
{
  if (Aeq)
    check_w0(w, n, meq, Aeq, beq);
  unsigned int m = mineq;
  if (lb)
    m += n;
  if (ub)
    m += n;
  double A[m*n];
  double b[m];
  set_Ab(A, b, m, n, mineq, Ain, bin, lb, ub);
  double lambda[meq+m];
  memset(lambda, 0, sizeof(lambda));

  unsigned int W[m];
  init_W(W, n, m, w, A, b);
  double p[n];

  unsigned int iter = 0;
  for (; iter<HECTR_ITER_MAX; iter++) {
    calc_p(p, lambda, n, m, meq, W, w, H, c, A, b, Aeq, beq);
    double p_norm = 0;
    for (unsigned int i=0; i<n; i++)
      p_norm += p[i]*p[i];

    if (p_norm < HECTR_TOLERANCE) {
      unsigned int mu_idx = 0;
      double mu = 0;
      calc_mu(&mu_idx, &mu, m, meq, W, lambda);
      if (mu > 0)
        break;
      else
        W[mu_idx] = 0;
    }
    else {
      unsigned int alpha_idx = 0;
      double alpha = 0;
      calc_alpha(&alpha_idx, &alpha, n, m, W, w, p, A, b);
      for (unsigned int i=0; i<n; i++)
        w[i] += alpha*p[i];
      if (alpha_idx<m)
        W[alpha_idx] = 1;
    }
  }
  printf("iteration = %i.\n", iter);
#if 0
  printf("lambda = [");
  for (unsigned int i=0; i<m+meq; i++)
    printf("%8.4f ", lambda[i]);
  printf("]\n");
#endif
  return 0;
}

END_DECLS
