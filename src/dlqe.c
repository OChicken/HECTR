/*
 * discrete-time linear-quadratic estimator.
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

BEGIN_DECLS

/**
 * dlqe - Linear-quadratic estimator for discrete-time systems.
 *
 * @G[m*n]: Kalman filter gain matrix
 * @l: dimension of y
 * @n: dimension of x
 * @m: dimension of u
 * @N: horizons
 * @A[n*n]: state matrix
 * @C[l*n]: output matrix
 * @Q[n*n]: state-cost weighted matrix
 * @R[l*l]: input-cost weighted matrix
 * @x0[n]: initial state
 */
void dlqe(double G[],
  const unsigned int l, const unsigned int n,
  const double A[], const double C[], const double Q[], const double R[])
{
  double AT[n*n], CT[n*l];
  for (unsigned int i=0, nn=n*n; i<nn; i++)
    AT[(i%n)*n+i/n] = A[i];
  for (unsigned int i=0, ln=l*n; i<ln; i++)
    CT[(i%n)*l+i/n] = C[i];
  /* Unique stabilizing solution of the discrete-time Riccati equation */
  double X[n*n];
  dare(X, n, l, AT, CT, Q, R);

  /* (X*C.T)[n*l] = X[n*n] * CT[n*l] */
  double XCT[n*l];
  memset(XCT, 0, sizeof(XCT));
  for (unsigned int i=0, all=n*l; i<all; i++)
    for (unsigned int j=0; j<n; j++)
      XCT[i] += X[(i/l)*n+j]*CT[j*l+i%l];

  /* (C*X*C.T+R)[l*l] = C[l*n] * (X*C.T)[n*l] + R[l*l] */
  double CXCTR[l*l];
  memset(CXCTR, 0, sizeof(CXCTR));
  for (unsigned int i=0, all=l*l; i<all; i++) {
    for (unsigned int j=0; j<n; j++)
      CXCTR[i] += C[(i/l)*n+j]*XCT[j*l+i%l];
    CXCTR[i] += R[i];
  }

  /* inv(C*X*C.T + R)[l*l] */
  double CXCTR_I[l*l];
  dgeinv(CXCTR_I, l, CXCTR);

  /* G[n*l] = (X*C.T)[n*l] * inv(C*X*C.T + R)[l*l] */
  memset(G, 0, n*l*sizeof(double));
  for (unsigned int i=0, all=n*l; i<all; i++)
    for (unsigned int j=0; j<l; j++)
      G[i] += XCT[(i/l)*l+j]*CXCTR_I[j*l+i%l];
}

END_DECLS
