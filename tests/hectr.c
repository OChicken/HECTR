/*
 * Test file: homomorphic encryption, decryption, addition and multiplication.
 * Copyright (C) shouran.ma@rwth-aachen.de
 *
 * This file is part of GPQHE.
 *
 * GPQHE is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2.1 of
 * the License, or (at your option) any later version.
 *
 * GPQHE is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this program; if not, see <http://www.gnu.org/licenses/>.
 */

#define _GNU_SOURCE
#include "../hectr/hectr.h"
#include "../gpqhe/gpqhe.h"
#include "../pmu/pmu.h"

static void test_quadprog()
{
  TEST_BEGIN();
  int err = 0;

  /**
   * 二次规划（1）：Lagrange法_qcyfred的博客-CSDN博客_二次规划 拉格朗日
   * https://blog.csdn.net/qcyfred/article/details/71598807
   *
   *  min  0.5 [x1,x2]*[4 1]*[x1] + [-1,-1]*[x1] = 2x1^2 + x1*x2 + x2^2 - x1-x2
   * x1,x2             [1 2] [x2]           [x2]
   * s.t [1,1]*[x1] - 1 = 0 (x1+x2-1=0)
   *           [x2]
   * 显然当x1=0.25时有极小值, 对应x2=0.75.
   *
   * 此问题的参数如下:
   * H = [4, 1;       c = [-1;
   *      1, 2] (2x2)      -1] (2x1)
   *
   * G = [1, 1] (1x2) g = -1
   */
  TEST_DO("LangLagR");{
  const unsigned int n=2, m=0, meq=1;
  double w[n];
  memset(w, 0, sizeof(w));

  double H[n*n];
  memcpy(H, (double[]){
    4,1,
    1,2
  }, sizeof(H));

  double c[n];
  memcpy(c, (double[]){
    -1,
    -1
  }, sizeof(c));

  double Aeq[meq*n];
  memcpy(Aeq, (double[]){
    1, 1
  }, sizeof(Aeq));

  double beq[meq];
  memcpy(beq, (double[]){
    -1
  }, sizeof(beq));

  err |= quadprog(w, n, m, meq, H, c, NULL, NULL, Aeq, beq, NULL, NULL);

  printf("w = [");
  for (unsigned int i=0; i<n; i++)
    printf("%7.4f ", w[i]);
  printf("]\n");

  }TEST_DONE();

  /**
   * Quadratic programming - MATLAB quadprog - MathWorks Deutschland
   * https://de.mathworks.com/help/optim/ug/quadprog.html
   * Quadratic Program with Linear Constraints
   */
  TEST_DO("Matlab: Quadratic Program with Linear Constraints");{
  const unsigned int n=2, m=3, meq=0;
  double w[n];
  memset(w, 0, sizeof(w));

  double H[n*n];
  memcpy(H, (double[]){
     1, -1,
    -1,  2
  }, sizeof(H));

  double c[n];
  memcpy(c, (double[]){
    -2,
    -6
  }, sizeof(c));

  double A[m*n];
  memcpy(A, (double[]){
     1, 1,
    -1, 2,
     2, 1
  }, sizeof(A));

  double b[m*1];
  memcpy(b, (double[]){
    -2,
    -2,
    -3
  }, sizeof(b));

  for (unsigned int i=0; i<n; i++)
    w[i] = -3;
  err |= quadprog(w, n, m, meq, H, c, A, b, NULL, NULL, NULL, NULL);

  printf("w = [");
  for (unsigned int i=0; i<n; i++)
    printf("%7.4f ", w[i]);
  printf("]\n");

  }TEST_DONE();

  /**
   * Quadratic programming - MATLAB quadprog - MathWorks Deutschland
   * https://de.mathworks.com/help/optim/ug/quadprog.html
   * Quadratic Program with Linear Equality Constraint
   */
  TEST_DO("Matlab: Quadratic Program with Linear Equality Constraint");{
  const unsigned int n=2, m=0, meq=1;
  double w[n];
  memset(w, 0, sizeof(w));

  double H[n*n];
  memcpy(H, (double[]){
     1,-1,
    -1, 2
  }, sizeof(H));

  double c[n];
  memcpy(c, (double[]){
    -2,
    -6
  }, sizeof(c));

  double Aeq[meq*n];
  memcpy(Aeq, (double[]){
    1,1
  }, sizeof(Aeq));

  double beq[meq];
  memcpy(beq, (double[]){
    0
  }, sizeof(beq));

  err |= quadprog(w, n, m, meq, H, c, NULL, NULL, Aeq, beq, NULL, NULL);

  printf("w = [");
  for (unsigned int i=0; i<n; i++)
    printf("%6.3f ", w[i]);
  printf("]\n");

  }TEST_DONE();

  /**
   * Quadratic programming - MATLAB quadprog - MathWorks Deutschland
   * https://de.mathworks.com/help/optim/ug/quadprog.html
   * Quadratic Minimization with Linear Constraints and Bounds
   */
  TEST_DO("Matlab: Quadratic Minimization with Linear Constraints and Bounds");{
  const unsigned int n=3, m=0, meq=1;
  double w[n];
  memset(w, 0, sizeof(w));

  double H[n*n];
  memcpy(H, (double[]){
     1, -1,  1,
    -1,  2, -2,
     1, -2,  4
  }, sizeof(H));

  double c[n];
  memcpy(c, (double[]){
     2,
    -3,
     1
  }, sizeof(c));

  double lb[n], ub[n];
  for (unsigned int i=0; i<n; i++)
    lb[i]=0, ub[i]=1;

  double Aeq[meq*n];
  memcpy(Aeq, (double[]){
    1,1,1
  }, sizeof(Aeq));

  double beq[meq];
  memcpy(beq, (double[]){
    -0.5
  }, sizeof(beq));

  err |= quadprog(w, n, m, meq, H, c, NULL, NULL, Aeq, beq, lb, ub);

  printf("w = [");
  for (unsigned int i=0; i<n; i++)
    printf("%6.3f ", w[i]);
  printf("]\n");

  }TEST_DONE();

  /**
   * Quadratic programming - MATLAB quadprog - MathWorks Deutschland
   * https://de.mathworks.com/help/optim/ug/quadprog.html
   * Return quadprog Objective Function Value
   */
  TEST_DO("Matlab: Return quadprog Objective Function Value");{
  const unsigned int n=3, m=1, meq=0;
  double w[n];
  memset(w, 0, sizeof(w));

  double H[n*n];
  memcpy(H, (double[]){
     1, -1,  1,
    -1,  2, -2,
     1, -2,  4
  }, sizeof(H));

  double c[n];
  memcpy(c, (double[]){
    -7,
    -12,
    -15
  }, sizeof(c));

  double A[m*n];
  memcpy(A, (double[]){
    1,1,1
  }, sizeof(A));

  double b[m*1];
  memcpy(b, (double[]){
    -3
  }, sizeof(b));

  err |= quadprog(w, n, m, meq, H, c, A, b, NULL, NULL, NULL, NULL);

  printf("w = [");
  for (unsigned int i=0; i<n; i++)
    printf("%6.4f ", w[i]);
  printf("]\n");

  }TEST_DONE();

  /**
   * Quadratic programming - MATLAB quadprog - MathWorks Deutschland
   * https://de.mathworks.com/help/optim/ug/quadprog.html
   * Examine quadprog Optimization Process
   */
  TEST_DO("Matlab: Examine quadprog Optimization Process");{
  const unsigned int n=3, m=0, meq=0;
  double w[n];
  memset(w, 0, sizeof(w));

  double H[n*n];
  memcpy(H, (double[]){
     2,  1, -1,
     1,  3,0.5,
    -1,0.5,  5
  }, sizeof(H));

  double c[n];
  memcpy(c, (double[]){
     4,
    -7,
     12
  }, sizeof(c));

  double lb[n], ub[n];
  for (unsigned int i=0; i<n; i++)
    lb[i]=0, ub[i]=1;

  err |= quadprog(w, n, m, meq, H, c, NULL, NULL, NULL, NULL, lb, ub);

  printf("w = [");
  for (unsigned int i=0; i<n; i++)
    printf("%8.4f ", w[i]);
  printf("]\n");

  }TEST_DONE();

  /**
   * Quadratic programming - MATLAB quadprog - MathWorks Deutschland
   * https://de.mathworks.com/help/optim/ug/quadprog.html
   * Return quadprog Lagrange Multipliers
   */
  TEST_DO("Matlab: Return quadprog Lagrange Multipliers");{
  const unsigned int n=3, m=1, meq=0;
  double w[n];
  memset(w, 0, sizeof(w));

  double H[n*n];
  memcpy(H, (double[]){
     1, -1,  1,
    -1,  2, -2,
     1, -2,  4
  }, sizeof(H));

  double c[n];
  memcpy(c, (double[]){
    -7,
    -12,
    -15
  }, sizeof(c));

  double A[m*n];
  memcpy(A, (double[]){
     1, 1, 1
  }, sizeof(A));

  double b[m*1];
  memcpy(b, (double[]){
     -3
  }, sizeof(b));

  double lb[n];
  for (unsigned int i=0; i<n; i++)
    lb[i]=0;

  err |= quadprog(w, n, m, meq, H, c, A, b, NULL, NULL, lb, NULL);

  printf("w = [");
  for (unsigned int i=0; i<n; i++)
    printf("%8.4f ", w[i]);
  printf("]\n");

  }TEST_DONE();

  /**
   * Function Reference: quadprog
   * https://octave.sourceforge.io/optim/function/quadprog.html
   */
  TEST_DO("Octave quadprog");{
  const unsigned int n=4, meq=1, m=3, dimy=5;
  double w[n];
  memset(w, 0, sizeof(w));

  double C[dimy*n];
  memcpy(C, (double[]){
    0.9501, 0.7620, 0.6153, 0.4057,
    0.2311, 0.4564, 0.7919, 0.9354,
    0.6068, 0.0185, 0.9218, 0.9169,
    0.4859, 0.8214, 0.7382, 0.4102,
    0.8912, 0.4447, 0.1762, 0.8936
  }, sizeof((C)));

  double d[dimy];
  memcpy(d, (double[]){
    0.0578,
    0.3528,
    0.8131,
    0.0098,
    0.1388
  }, sizeof(d));

  double H[n*n];
  memset(H, 0, sizeof(H));
  for (unsigned int i=0; i<n*n; i++)
    for (unsigned int j=0; j<dimy; j++)
      H[i] += C[j*n+(i/n)] * C[j*n+(i%n)];

  double c[n];
  memset(c, 0, sizeof(c));
  for (unsigned int i=0; i<n; i++)
    for (unsigned int j=0; j<dimy; j++)
      c[i] += -C[j*n+i] * d[j];

  double A[m*n];
  memcpy(A, (double[]){
     0.2027, 0.2721, 0.7467, 0.4659,
     0.1987, 0.1988, 0.4450, 0.4186,
     0.6037, 0.0152, 0.9318, 0.8462
  }, sizeof(A));

  double b[m*1];
  memcpy(b, (double[]){
    -0.5251,
    -0.2026,
    -0.6721,
  }, sizeof(b));

  double lb[n], ub[n];
  for (unsigned int i=0; i<n; i++)
    lb[i]=-0.1, ub[i]=1;

  double Aeq[meq*n];
  memcpy(Aeq, (double[]){
    3,5,7,9
  }, sizeof(Aeq));

  double beq[meq*1];
  memcpy(beq, (double[]){
    -4
  }, sizeof(beq));

  err |= quadprog(w, n, m, meq, H, c, A, b, Aeq, beq, lb, ub);

  printf("w = [");
  for (unsigned int i=0; i<n; i++)
    printf("%8.5f ", w[i]);
  printf("]\n");

  }TEST_DONE();

  TEST_DO("PyAdvancedControl/mpc_tracking/pyecosqp.py/test1");{
  const unsigned int n=2, meq=0, m=5;
  double w[n];
  memset(w, 0, sizeof(w));
  double H[n*n];
  memset(H, 0, sizeof(H));
  H[0]=1;
  double c[n];
  c[0]=3, c[1]=4;
  double A[m*n];
  memcpy(A, (double[]){
    -1.0,  0.0,
     0. , -1.0,
    -1.0, -3.0,
     2.0,  5.0,
     3.0,  4.0
  }, sizeof(A));
  double b[m];
  memcpy(b, (double[]){
    0.0, 0.0, 15.0, -100.0, -80.0
  }, sizeof(b));
  err |= quadprog(w, n, m, meq, H, c, A, b, NULL, NULL, NULL, NULL);
  printf("w = [");
  for (unsigned int i=0; i<n; i++)
    printf("%8.5f ", w[i]);
  printf("]\n");
  }TEST_DONE();

  TEST_DO("PyAdvancedControl/mpc_tracking/pyecosqp.py/test2");{
  const unsigned int n=9, meq=6, m=0;
  double w[n];
  memset(w, 0, sizeof(w));
  double H[n*n];
  memset(H, 0, sizeof(H));
  for (unsigned int i=0; i<n; i++)
    H[i*n+i]=1;
  double c[n];
  memset(c, 0, sizeof(c));
  double Aeq[meq*n];
  memcpy(Aeq, (double[]){
     1.,  0.,  0.,  1. ,  0. ,  0. ,  0. , 0., 0.,
    -2., -0., -0.,  0. ,  1. ,  0. ,  0. , 0., 0.,
     0.,  1.,  0., -0.8, -1. ,  1. ,  0. , 0., 0.,
    -0., -2., -0.,  0. , -0.9,  0. ,  1. , 0., 0.,
     0.,  0.,  1.,  0. ,  0. , -0.8, -1. , 1., 0.,
    -0., -0., -2.,  0. ,  0. ,  0. , -0.9, 0., 1.
  }, sizeof(Aeq));
  double beq[meq];
  memset(beq, 0, sizeof(beq));
  beq[0]=-2.8, beq[1]=-1.8;
  err |= quadprog(w, n, m, meq, H, c, NULL, NULL, Aeq, beq, NULL, NULL);
  printf("w = [");
  for (unsigned int i=0; i<n; i++)
    printf("%8.5f ", w[i]);
  printf("]\n");
  }TEST_DONE();

  TEST_DO("PyAdvancedControl/mpc_tracking/pyecosqp.py/test3");{
  const unsigned int n=9, meq=6, m=6;
  double w[n];
  memset(w, 0, sizeof(w));
  double H[n*n];
  memset(H, 0, sizeof(H));
  for (unsigned int i=0; i<n; i++)
    H[i*n+i]=1;
  double c[n];
  memset(c, 0, sizeof(c));
  double Aeq[meq*n];
  memcpy(Aeq, (double[]){
     1.,  0.,  0.,  1. ,  0. ,  0. ,  0. , 0., 0.,
    -2., -0., -0.,  0. ,  1. ,  0. ,  0. , 0., 0.,
     0.,  1.,  0., -0.8, -1. ,  1. ,  0. , 0., 0.,
    -0., -2., -0.,  0. , -0.9,  0. ,  1. , 0., 0.,
     0.,  0.,  1.,  0. ,  0. , -0.8, -1. , 1., 0.,
    -0., -0., -2.,  0. ,  0. ,  0. , -0.9, 0., 1.
  }, sizeof(Aeq));
  double beq[meq];
  memset(beq, 0, sizeof(beq));
  beq[0]=-2.8, beq[1]=-1.8;
  double A[m*n];
  memcpy(A, (double[]){
     1.,  0.,  0., 0., 0., 0., 0., 0., 0.,
     0.,  1.,  0., 0., 0., 0., 0., 0., 0.,
     0.,  0.,  1., 0., 0., 0., 0., 0., 0.,
    -1., -0., -0., 0., 0., 0., 0., 0., 0.,
    -0., -1., -0., 0., 0., 0., 0., 0., 0.,
    -0., -0., -1., 0., 0., 0., 0., 0., 0.
  }, sizeof(A));
  double b[m];
  for (unsigned int i=0; i<m; i++)
    b[i] = -0.7;
  err |= quadprog(w, n, m, meq, H, c, A, b, Aeq, beq, NULL, NULL);
  printf("w = [");
  for (unsigned int i=0; i<n; i++)
    printf("%8.5f ", w[i]);
  printf("]\n");
  }TEST_DONE();

  TEST_END();
}

/* Steady-state values */
static const double cs = 0.878; /* kmol/m^3 */
static const double Ts = 324.5; /* K */
static const double hs = 0.659; /* m */
static const double Fs = 0.1;   /* m^3/min */
static const double Tcs = 300;  /* K */
static const double F0s = 0.1;  /* m^3/min */

static void test_cstr_ode()
{
  /* Steady State Initial Conditions for the States */
  double x_ode45[3]  = {cs,Ts,hs}; /* x_ss */
  double x_ode15s[3] = {cs,Ts,hs}; /* x_ss */
  /* Steady State Initial Condition for the Control */
  double u_ss = 290;
  /* Open Loop Step Change */
  double u[2] = {u_ss, Fs};
  double p[1] = {F0s};

  /* Final Time (sec) */
  double tf = 5;

  double dt = 1;
  unsigned int N = tf/dt;
  FILE *fd = fopen("results/cstr-ode.txt", "w");
  fprintf(fd, "%9.6f %9.6f %9.6f %9.6f %9.6f\n", 0.,
    x_ode45[0], x_ode45[1], x_ode15s[0], x_ode15s[1]);
  for (unsigned int i=1; i<=N; i++) {
    ode45(x_ode45, 3, x_ode45, u, p, cstr_ode, dt);
    ode15s(x_ode15s, 3, x_ode15s, u, p, cstr_ode, cstr_jacobian, dt);
    fprintf(fd, "%9.6f %9.6f %9.6f %9.6f %9.6f\n", i*dt,
      x_ode45[0], x_ode45[1], x_ode15s[0], x_ode15s[1]);
  }
  fclose(fd);
}

static void calc_xnew(double x[],
  const unsigned int n, const unsigned int m, const unsigned int N,
  const double A[], const double B[], const double x0[], const double u[])
{
  const unsigned int nNp1 = n*(N+1);
  memset(x, 0, nNp1*sizeof(double));
  for (unsigned int i=0; i<n; i++)
    x[i] = x0[i];
  double Ax[n],Bu[n];
  for (unsigned int k=0; k<N; k++) {
    memset(Ax, 0, sizeof(Ax));
    memset(Bu, 0, sizeof(Bu));
    for (unsigned int i=0; i<n; i++) {
      for (unsigned int j=0; j<n; j++)
        Ax[i] += A[i*n+j]*x[k*n+j];
      for (unsigned int j=0; j<m; j++)
        Bu[i] += B[i*m+j]*u[k*m+j];
      x[(k+1)*n+i] = Ax[i]+Bu[i];
    }
  }
}

static void test_mpc_tracking()
{
  TEST_BEGIN();
  const unsigned int n=2, m=1, l=n, N=30;
  double A[n*n];
  memcpy(A, (double[]){
    0.8, 1.0,
    0  , 0.9
  }, sizeof(A));
  double B[n*m];
  memcpy(B, (double[]){-1.0, 2.0}, sizeof(B));
  double C[l*n];
  memcpy(C, (double[]){
    1,0,
    0,1
  }, sizeof(C));
  double Q[n*n];
  memset(Q, 0, sizeof(Q));
  for (unsigned int i=0; i<n; i++)
    Q[i*n+i]=1;
  double R[m*m];
  memset(R, 0, sizeof(R));
  for (unsigned int i=0; i<m; i++)
    R[i*m+i]=1;
  double x0[n];
  memcpy(x0, (double[]){0.0, -1.0}, sizeof(x0));
  double u0[m];
  memcpy(u0, (double[]){-0.1}, sizeof(u0));
  double u[m*N];
  double y[l*(N+1)];
  double r[n*(N+1)];
  for (unsigned int i=0; i<N+1; i++)
    r[i*n+0]=1, r[i*n+1]=0.25;

  TEST_DO("PyAdvancedControl/mpc_tracking/mpc_tracking.py/test5");{
  ctr_mpc(u, l, n, m, N, A, B, C, Q, R, x0, u0, r, u0, NULL, NULL, NULL, NULL, NULL, NULL);
  calc_xnew(y, n, m, N, A, B, x0, u);
  FILE *fd = fopen("results/mpc-tracking-5.txt", "w");
  for (unsigned int k=0; k<N+1; k++)
    fprintf(fd, "%2u %12.8f %12.8f %12.8f\n", k, (k<N)? u[k *m+0] : u[(N-1)*m+0], y[k*l+0], y[k*l+1]);
  fclose(fd);
  }TEST_DONE();

  TEST_DO("PyAdvancedControl/mpc_tracking/mpc_tracking.py/test6");{
  double dumin[m], dumax[m];
  dumin[0]=-0.5, dumax[0]=+0.5;
  ctr_mpc(u, l, n, m, N, A, B, C, Q, R, x0, u0, r, u0, dumin, dumax, NULL, NULL, NULL, NULL);
  calc_xnew(y, n, m, N, A, B, x0, u);
  FILE *fd = fopen("results/mpc-tracking-6.txt", "w");
  for (unsigned int k=0; k<N+1; k++)
    fprintf(fd, "%2u %12.8f %12.8f %12.8f\n", k, (k<N)? u[k *m+0] : u[(N-1)*m+0], y[k*l+0], y[k*l+1]);
  fclose(fd);
  }TEST_DONE();

  TEST_DO("PyAdvancedControl/mpc_tracking/mpc_tracking.py/test7");{
  double dumin[m], dumax[m];
  dumin[0]=-0.3, dumax[0]=+0.2;
  ctr_mpc(u, l, n, m, N, A, B, C, Q, R, x0, u0, r, u0, dumin, dumax, NULL, NULL, NULL, NULL);
  calc_xnew(y, n, m, N, A, B, x0, u);
  FILE *fd = fopen("results/mpc-tracking-7.txt", "w");
  for (unsigned int k=0; k<N+1; k++)
    fprintf(fd, "%2u %12.8f %12.8f %12.8f\n", k, (k<N)? u[k *m+0] : u[(N-1)*m+0], y[k*l+0], y[k*l+1]);
  fclose(fd);
  }TEST_DONE();

  for (unsigned int i=0; i<N+1; i++)
    r[i*n+0]=0, r[i*n+1]=0;

  TEST_DO("PyAdvancedControl/mpc_tracking/mpc_tracking.py/test8");{
  double dumin[m], dumax[m];
  dumin[0]=-0.3, dumax[0]=+0.2;
  ctr_mpc(u, l, n, m, N, A, B, C, Q, R, x0, u0, r, u0, dumin, dumax, NULL, NULL, NULL, NULL);
  calc_xnew(y, n, m, N, A, B, x0, u);
  FILE *fd = fopen("results/mpc-tracking-8.txt", "w");
  for (unsigned int k=0; k<N+1; k++)
    fprintf(fd, "%2u %12.8f %12.8f %12.8f\n", k, (k<N)? u[k *m+0] : u[(N-1)*m+0], y[k*l+0], y[k*l+1]);
  fclose(fd);
  }TEST_DONE();

  TEST_DO("PyAdvancedControl/mpc_tracking/mpc_tracking.py/test9");{
  double umin[m], umax[m];
  umin[0]=-0.3, umax[0]=+0.1;
  ctr_mpc(u, l, n, m, N, A, B, C, Q, R, x0, u0, r, u0, NULL, NULL, umin, umax, NULL, NULL);
  calc_xnew(y, n, m, N, A, B, x0, u);
  FILE *fd = fopen("results/mpc-tracking-9.txt", "w");
  for (unsigned int k=0; k<N+1; k++)
    fprintf(fd, "%2u %12.8f %12.8f %12.8f\n", k, (k<N)? u[k *m+0] : u[(N-1)*m+0], y[k*l+0], y[k*l+1]);
  fclose(fd);
  }TEST_DONE();

  TEST_DO("PyAdvancedControl/mpc_tracking/mpc_tracking.py/test11");{
  double xmin[n], xmax[n];
  xmin[0]=-1.5, xmax[0]=+0.5;
  xmin[1]=-2.5, xmax[1]=+0.2;
  ctr_mpc(u, l, n, m, N, A, B, C, Q, R, x0, u0, r, u0, NULL, NULL, NULL, NULL, xmin, xmax);
  calc_xnew(y, n, m, N, A, B, x0, u);
  FILE *fd = fopen("results/mpc-tracking-11.txt", "w");
  for (unsigned int k=0; k<N+1; k++)
    fprintf(fd, "%2u %12.8f %12.8f %12.8f\n", k, (k<N)? u[k *m+0] : u[(N-1)*m+0], y[k*l+0], y[k*l+1]);
  fclose(fd);
  }TEST_DONE();

  TEST_DO("PyAdvancedControl/mpc_tracking/mpc_tracking.py/test12");{
  double dumin[m], dumax[m];
  dumin[0]=-0.5, dumax[0]=+0.5;
  double xmin[n], xmax[n];
  xmin[0]=-1.5, xmax[0]=+0.5;
  xmin[1]=-2.5, xmax[1]=+0.2;
  ctr_mpc(u, l, n, m, N, A, B, C, Q, R, x0, u0, r, u0, dumin, dumax, NULL, NULL, xmin, xmax);
  calc_xnew(y, n, m, N, A, B, x0, u);
  FILE *fd = fopen("results/mpc-tracking-12.txt", "w");
  for (unsigned int k=0; k<N+1; k++)
    fprintf(fd, "%2u %12.8f %12.8f %12.8f\n", k, (k<N)? u[k *m+0] : u[(N-1)*m+0], y[k*l+0], y[k*l+1]);
  fclose(fd);
  }TEST_DONE();

  TEST_END();
}

static void test_cstr_mpc()
{
  const unsigned int nx=3, nu=2, np=1, ny=nx;
  /* Steady state initial conditions for the states, control and system parameters */
  double xs[nx], us[nu], ps[np];
  memcpy(xs, (double[]){cs, Ts, hs}, sizeof(xs));
  memcpy(us, (double[]){Tcs, Fs}   , sizeof(us));
  memcpy(ps, (double[]){F0s}       , sizeof(ps));

  /* regulator */
  const double dt = 1;
  double A[nx*nx], B[nx*nu], Bp[nx*np], C[ny*nx];
  cstr_linearize(A, B, Bp, nx, nu, np, xs, us, ps, dt);
  memset(C, 0, sizeof(C));
  for (unsigned int i=0; i<nx; i++)
    C[i*nx+i]=1;

  /* disturbance model with offset */
  const unsigned int nd=2;
  double Bd[nx*nd], Cd[ny*nd];
  memset(Bd, 0, sizeof(Bd));
  memcpy(Cd, (double[]){
    1,0,
    0,0,
    0,1
  }, sizeof(Cd));

  /* selector */
  double Hr[nu*ny];
  memcpy(Hr, (double[]){
    1,0,0,
    0,0,1
  }, sizeof(Hr)); /* Control Concentration and height. */

  /* Closed-loop simulation */
  const unsigned int N = 40;
  double x[nx*(N+1)], u[nu*N];

  /* Disturbance in inlet flow. */
  double p[np*N];
  memset(p, 0, sizeof(p));
  for (unsigned int i=9; i<np*N; i++)
    p[i]=0.1*F0s;

  ctr_simulate(x, u, p, nx, nu, np, ny, nd, A, B, C, Bd, Cd, Hr,
    xs, us, ps, cstr_ode, cstr_jacobian, dt, N);
  /* save */
  FILE *fd = fopen("results/cstr-mpc.txt", "w");
  for (unsigned int k=0; k<N+1; k++)
    fprintf(fd, "%2u %13g %13g %13g %13g %13g\n",
      k, x[k*nx+0], x[k*nx+1], x[k*nx+2], (k<N)? u[k*nu+0]:u[(N-1)*nu+0], (k<N)? u[k*nu+1]:u[(N-1)*nu+1]);
  fclose(fd);
  fd = fopen("results/cstr-mpc.bin", "wb");
  for (unsigned int k=0; k<N+1; k++) {
    fwrite(&k, sizeof(unsigned int), 1, fd);
    fwrite(&x[k*nx], sizeof(double), nx, fd);
    fwrite((k<N)? &u[k*nu] : &u[(N-1)*nu], sizeof(double), nu, fd);
  }
  fclose(fd);
}

static void test_cstr_hempc()
{
  const unsigned int nx=3, nu=2, np=1, ny=nx;
  /* Steady state initial conditions for the states, control and system parameters */
  double xs[nx], us[nu], ps[np];
  memcpy(xs, (double[]){cs, Ts, hs}, sizeof(xs));
  memcpy(us, (double[]){Tcs, Fs}   , sizeof(us));
  memcpy(ps, (double[]){F0s}       , sizeof(ps));

  /* regulator */
  const double dt = 1;
  double A[nx*nx], B[nx*nu], Bp[nx*np], C[ny*nx];
  cstr_linearize(A, B, Bp, nx, nu, np, xs, us, ps, dt);
  memset(C, 0, sizeof(C));
  for (unsigned int i=0; i<nx; i++)
    C[i*nx+i]=1;

  /* disturbance model with offset */
  const unsigned int nd=2;
  double Bd[nx*nd], Cd[ny*nd];
  memset(Bd, 0, sizeof(Bd));
  memcpy(Cd, (double[]){
    1,0,
    0,0,
    0,1
  }, sizeof(Cd));

  /* selector */
  double Hr[nu*ny];
  memcpy(Hr, (double[]){
    1,0,0,
    0,0,1
  }, sizeof(Hr)); /* Control Concentration and height. */

  /* Closed-loop simulation */
  const unsigned int N = 40;
  double x[nx*(N+1)], u[nu*N];

  /* Disturbance in inlet flow. */
  double p[np*N];
  memset(p, 0, sizeof(p));
  for (unsigned int i=9; i<np*N; i++)
    p[i]=0.1*F0s;

  hectr_simulate(x, u, p, nx, nu, np, ny, nd, A, B, C, Bd, Cd, Hr,
    xs, us, ps, cstr_ode, cstr_jacobian, dt, N);
  /* save */
  FILE *fd = fopen("results/cstr-hempc.txt", "w");
  for (unsigned int k=0; k<N+1; k++)
    fprintf(fd, "%2u %13g %13g %13g %13g %13g\n",
      k, x[k*nx+0], x[k*nx+1], x[k*nx+2], (k<N)? u[k*nu+0]:u[(N-1)*nu+0], (k<N)? u[k*nu+1]:u[(N-1)*nu+1]);
  fclose(fd);
  fd = fopen("results/cstr-hempc.bin", "wb");
  for (unsigned int k=0; k<N+1; k++) {
    fwrite(&k, sizeof(unsigned int), 1, fd);
    fwrite(&x[k*nx], sizeof(double), nx, fd);
    fwrite((k<N)? &u[k*nu] : &u[(N-1)*nu], sizeof(double), nu, fd);
  }
  fclose(fd);
}

static void test_cstr_cmp()
{
  const unsigned int nx=3, nu=2, N=40;
  double xp[nx], xc[nx], xdiff[nx], up[nu], uc[nu], udiff[nu];
  FILE *fdp, *fdc, *fd;
  fdp = fopen("results/cstr-mpc.bin", "rb");
  fdc = fopen("results/cstr-hempc.bin", "rb");
  fd  = fopen("results/cstr-cmp.bin", "wb");
  for (unsigned int k=0; k<N+1; k++) {
    fread(&k, sizeof(unsigned int), 1, fdp);
    fread(&k, sizeof(unsigned int), 1, fdc);
    fread(xp, sizeof(double), nx, fdp);
    fread(xc, sizeof(double), nx, fdc);
    fread(up, sizeof(double), nu, fdp);
    fread(uc, sizeof(double), nu, fdc);
    for (unsigned int i=0; i<nx; i++)
      xdiff[i] = fabs(xp[i]-xc[i]);
    for (unsigned int i=0; i<nu; i++)
      udiff[i] = fabs(up[i]-uc[i]);
    fwrite(&k, sizeof(unsigned int), 1, fd);
    fwrite(xdiff, sizeof(double), nx, fd);
    fwrite(udiff, sizeof(double), nu, fd);
  }
  fclose(fdp);
  fclose(fdc);
  fclose(fd);
}

/**
 * A = [0, 1  , 0,0;
 *      0,-0.1, 3,0;
 *      0, 0  , 0,1;
 *      0,-0.5,30,0,];
 * B = [0;2;0;5];
 * C = [0,1,0,0;
 *      0,0,1,0];
 * Q = C'*C;
 * R = [1];
 * 
 */
static void test_inverted_pendulum_mpc_control()
{
  const double l_bar = 2.0; /* length of bar */
  const double mcar  = 1.0; /* mass of car   */
  const double mball = 0.3; /* mass of ball  */
  const double g     = 9.8; /* gravitation */
  const unsigned int n = 4;
  const unsigned int m = 1;
  const unsigned int l = 2;
  const unsigned int N = 30;
  const double dt = 0.1;
  const unsigned int ln = l*n;
  const unsigned int nn = n*n;
  const unsigned int nm = n*m;
  const unsigned int mN = m*N;
  const unsigned int nNp1 = n*(N+1);
  /* A */
  double A[nn], Ad[nn];
  memcpy(A, (double[]){
    0.0, 1.0,                         0.0, 0.0,
    0.0, 0.0,                mball*g/mcar, 0.0,
    0.0, 0.0,                         0.0, 1.0,
    0.0, 0.0, g*(mcar+mball)/(l_bar*mcar), 0.0,
  }, sizeof(A));
  /* B */
  double B[nm], Bd[nm];
  memcpy(B, (double[]){
    0,
    1/mcar,
    0,
    1/(l_bar*mcar)
  }, sizeof(B));
  double expmA[nn], expmB[nn];
  ctr_c2d(expmA, expmB, n, A, dt);
  memcpy(Ad, expmA, sizeof(expmA));
  memset(Bd, 0, n*m*sizeof(double));
  for (unsigned int i=0, all=n*m; i<all; i++)
    for (unsigned int j=0; j<n; j++)
      Bd[i] += expmB[(i/m)*n+j]*B[j*m+i%m];
  /* C */
  double C[ln];
  memcpy(C, (double[]){
    0,1,0,0,
    0,0,1,0
  }, sizeof(C));
  double Q[l*l];
  memset(Q, 0, sizeof(Q));
  for (unsigned int i=0; i<l; i++)
    Q[i*l+i]=1;
  double R[m*m];
  memcpy(R, (double[]){
    0.01
  }, sizeof(R));
  double x0[n];
  memcpy(x0, (double[]){
    0,0,0.3,0
  }, sizeof(x0));
  double x[nNp1],u[mN];
  double u0[m];
  u0[0]=0;
  double r[nNp1];
  for (unsigned int i=0; i<nNp1; i++)
    r[i]=0;
  ctr_mpc(u, l, n, m, N, Ad, Bd, C, Q, R, x0, u0, r, u0, NULL, NULL, NULL, NULL, NULL, NULL);
  calc_xnew(x, n, m, N, Ad, Bd, x0, u);
  FILE *fd = fopen("results/inverted-pendulum-mpc-control.txt", "w");
  for (unsigned int k=0; k<N+1; k++)
    fprintf(fd, "%2u %12.8f %12.8f %12.8f %12.8f %12.8f\n",
      k, (k<N)? u[k *m+0] : u[(N-1)*m+0], x[k*n+0], x[k*n+1], x[k*n+2], x[k*n+3]);
  fclose(fd);
}

int main(int argc, char *argv[])
{
  if (argc==1) {
    fprintf(stderr, "usage: %s [quadprog/] [sk/pk] "
      "--logn=num --logq=num --logp=num --iter=num\n", argv[0]);
    exit(1);
  }
  if (!strcmp(argv[1], "quadprog"))
    test_quadprog();
  if (!strcmp(argv[1], "cstr-ode")) {
    test_cstr_ode();
    FILE *fd = popen("cd results && gnuplot cstr-ode.gp", "w");
    pclose(fd);
  }
  if (!strcmp(argv[1], "mpc-tracking")) {
    test_mpc_tracking();
    FILE *fd = popen("cd results && gnuplot mpc-tracking.gp", "w");
    pclose(fd);
  }
  if (!strcmp(argv[1], "inverted-pendulum-mpc-control")) {
    test_inverted_pendulum_mpc_control();
    FILE *fd = popen("cd results && gnuplot inverted-pendulum-mpc-control.gp", "w");
    pclose(fd);
  }
  if (!strcmp(argv[1], "cstr-mpc")) {
    test_cstr_mpc();
    FILE *fd = popen("cd results && gnuplot cstr-mpc.gp", "w");
    pclose(fd);
  }
  if (!strcmp(argv[1], "cstr-hempc")) {
    test_cstr_hempc();
    FILE *fd = popen("cd results && gnuplot cstr-hempc.gp", "w");
    pclose(fd);
  }
  if (!strcmp(argv[1], "cstr-cmp")) {
    test_cstr_cmp();
    FILE *fd = popen("cd results && gnuplot cstr-cmp.gp", "w");
    pclose(fd);
  }
  return 0;
}
