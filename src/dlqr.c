/*
 * discrete-time linear-quadratic regulator.
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
 * dlqr - Linear-quadratic regulator for discrete-time systems.
 *
 * @G[m*n]: state feedback matrix
 * @n: dimension of x
 * @m: dimension of u
 * @A[n*n]: state matrix
 * @B[n*m]: input matrix
 * @Q[n*n]: state-cost weighted matrix
 * @R[m*m]: input-cost weighted matrix
 */
void dlqr(double G[], double X[],
  const unsigned int n, const unsigned m,
  const double A[], const double B[], const double Q[], const double R[])
{
  /* Unique stabilizing solution of the discrete-time Riccati equation */
  dare(X, n, m, A, B, Q, R);

  /* B.T*X[m*n(N+1)] = B.T[m*n(N+1)] * X[n(N+1)*n(N+1)] */
  double BTX[m*n];
  memset(BTX, 0, sizeof(BTX));
  for (unsigned int i=0, all=m*n; i<all; i++)
    for (unsigned int j=0; j<n; j++)
      BTX[i]+=B[j*m+i/n]*X[j*n+i%n];

  /* (B.T*X*B+R)[m*m] = B.T*X[m*n(N+1)] * B[n(N+1)*m] + R[m*m] */
  double BTXBR[m*m];
  memcpy(BTXBR, R, sizeof(BTXBR));
  for (unsigned int i=0, all=m*m; i<all; i++)
    for (unsigned int j=0; j<n; j++)
      BTXBR[i]+=BTX[(i/m)*n+j]*B[j*m+i%m];

  /* inv(B.T*X*B+R)[m * m] */
  double BTXBR_inv[m*m];
  dgeinv(BTXBR_inv, m, BTXBR);

  /* B.T*X*A[m*n] = B.T*X[m*n(N+1)] * A[n(N+1)*n] */
  double BTXA[m*n];
  memset(BTXA, 0, sizeof(BTXA));
  for (unsigned int i=0, all=m*n; i<all; i++)
    for (unsigned int j=0; j<n; j++)
      BTXA[i]+=BTX[(i/n)*n+j]*A[j*n+i%n];

  /* G[m*n] = inv(B.T*X*B+R)[m*m] * B.T*X*A[m*n] */
  memset(G, 0, m*n*sizeof(double));
  for (unsigned int i=0, all=m*n; i<all; i++)
    for (unsigned int j=0; j<m; j++)
      G[i]+=BTXBR_inv[(i/n)*m+j]*BTXA[j*n+i%n];
}

END_DECLS
