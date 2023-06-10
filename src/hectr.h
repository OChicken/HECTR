/*
 * hectr.h - GNU post-quantum homomorphic-encryption library interface.
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

#ifndef HECTR_H
#define HECTR_H

#include "config.h"

#ifndef CONFIG_H
# error config.h must be included before types.h
#endif

#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <lapacke.h>
#include "../gpqhe/gpqhe.h"

BEGIN_DECLS

#define HECTR_TOLERANCE 1e-10
#define HECTR_SMALL 1e-5
#define HECTR_ITER_MAX 10000

void ode45(double xnext[],
  const unsigned int nx,
  const double x[], const double u[], const double p[],
  void (*const ode)(double[], const double[], const double[], const double[]),
  const double dt);
void ode15s(double xnext[],
  const unsigned int n,
  const double x[], const double u[], const double p[],
  void (*const ode     )(double[], const double[], const double[], const double[]),
  void (*const jacobian)(double[], const double[], const double[], const double[]),
  const double dt);

/* cstr.c */
void cstr_ode(double xdot[], const double x[], const double u[], const double p[]);
void cstr_jacobian(double J[], const double x[],
  __attribute__((unused))const double u[], const double p[]);
void cstr_linearize(double A[], double B[], double Bp[],
  const unsigned int nx, const unsigned int nu, const unsigned int np,
  const double xs[], const double us[], const double ps[], const double dt);

/* quadprog.c */
int quadprog(double w[],
  const unsigned int n, const unsigned int mineq, const unsigned int meq,
  const double H[], const double c[],
  const double Ain[], const double bin[],
  const double Aeq[], const double beq[],
  const double lb[], const double ub[]);

/* ctr.c */
void ctr_c2d(double expA[], double expB[],
  const unsigned int n, const double jacA[], const double dt);
void ctr_estimator(double Lx[], double Ld[],
  const unsigned int nx, const unsigned int nu, const unsigned int nd, const unsigned int ny,
  const double A[], const double B[], const double C[],
  const double Bd[], const double Cd[],
  const double xs[]);
void ctr_selector(double Ginv[],
  const unsigned int nx, const unsigned int nu, const unsigned int ny,
  const double A[], const double B[], const double C[], const double H[]);
void ctr_measure(double y[],
  const unsigned int ny, const unsigned int nx,
  const double C[], const double x[]);
void ctr_measure_forward(double xhat[], double dhat[],
  const unsigned int nx, const unsigned int nd, const unsigned int ny,
  const double C[],  const double Cd[],
  const double Lx[], const double Ld[],
  const double y[], const double xhatm[], const double dhatm[]);
void ctr_select(double xr[], double ur[],
  const unsigned int nx, const unsigned int nu, const unsigned int nd, const unsigned int ny,
  const double Bd[], const double Cd[],
  const double H[], const double Ginv[],
  const double dhat[], const double rsp[]);
void ctr_control(double u[],
  const unsigned int nx, const unsigned int nu,
  const double G[], const double xhat[], const double xr[], const double ur[]);
void ctr_estimate(double xhatm[], double dhatm[],
  const unsigned int nx, const unsigned int nu, const unsigned int nd,
  const double A[], const double B[], const double Bd[],
  const double xhat[], const double dhat[], const double u[]);
void ctr_actuate(double xnext[],
  const unsigned int nx, const unsigned int nu, const unsigned int np,
  const double x[], const double u[], const double p[],
  const double xs[], const double us[], const double ps[],
  void (*const ode     )(double[], const double[], const double[], const double[]),
  void (*const jacobian)(double[], const double[], const double[], const double[]),
  const double dt);
void ctr_simulate(double x[], double u[], double p[],
  const unsigned int nx, const unsigned int nu, const unsigned int np, /* state, control, parameters */
  const unsigned int ny, const unsigned int nd, /* observe, disturbance */
  const double A[], const double B[], const double C[],
  const double Bd[], const double Cd[], const double Hr[],
  const double xs[], const double us[], const double ps[],
  void (*const ode     )(double[], const double[], const double[], const double[]),
  void (*const jacobian)(double[], const double[], const double[], const double[]),
  const double dt, const unsigned int N);
void hectr_simulate(double x[], double u[], double p[],
  const unsigned int nx, const unsigned int nu, const unsigned int np, /* state, control, parameters */
  const unsigned int ny, const unsigned int nd, /* observe, disturbance */
  const double A[], const double B[], const double C[],
  const double Bd[], const double Cd[], const double Hr[],
  const double xs[], const double us[], const double ps[],
  void (*const ode     )(double[], const double[], const double[], const double[]),
  void (*const jacobian)(double[], const double[], const double[], const double[]),
  const double dt, const unsigned int N);

/* dare.c */
void dare(double X[],
  const unsigned int n, const unsigned int m,
  const double A[], const double B[], const double Q[], const double R[]);

/* dlqr.c */
void dlqr(double G[], double X[],
  const unsigned int n, const unsigned int m,
  const double A[], const double B[], const double Q[], const double R[]);

/* dlqe.c */
void dlqe(double G[],
  const unsigned int l, const unsigned int n,
  const double A[], const double C[], const double Q[], const double R[]);

/* matrices.c */
double eps(double a);
int dgeinv(double Ainv[], const unsigned int n, const double A[]);
int zgeinv(_Complex double Ainv[], const unsigned int n, const _Complex double A[]);
int dgepinv(double Ainv[], const unsigned int m, const unsigned int n, const double A[]);
int dexpm(double expA[], const unsigned int n, const double A[]);
void d2z_vector(_Complex double vz[],
  const unsigned int n, const unsigned int slots,
  const double vd[]);
void d2z_matrix(_Complex double Az[],
  const unsigned int row, const unsigned int col, const unsigned int slots,
  const double Ad[]);

/* mpc.c */
void ctr_mpc(double u[],
  const unsigned int l, const unsigned int n, const unsigned int m, const unsigned int N,
  const double A[], const double B[], const double C[], const double Q[], const double R[],
  const double xhat[], const double uhat[], const double xr[], const double ur[],
  const double dumin[], const double dumax[],
  const double  umin[], const double  umax[],
  const double  xmin[], const double  xmax[]);
void ctr_hempc(he_ct_t *up,
  const unsigned int l, const unsigned int n, const unsigned int m, const unsigned int N,
  const unsigned int slots,
  const double A[], const double B[], const double C[], const double Q[], const double R[],
  const he_ct_t *xhat, const he_ct_t *uhat, const he_ct_t *xr, const he_ct_t *ur,
  const he_evk_t *rk);

END_DECLS

#endif /* HECTR_H */
