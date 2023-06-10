/*
 * discrete-time  Riccati equation solver.
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

BEGIN_DECLS

/**
 * dare - discrete-time Riccati equation solver:
 * X = A.T*X*A - A.T*X*B * inv(R+B.T*X*B) * B.T*X*A + Q
 *
 * @X[n*n]: unique stabilizing solution of the discrete-time Riccati equation
 * @n: dimension of x
 * @m: dimension of u
 * @A[n*n]: state matrix
 * @B[n*m]: input matrix
 * @Q[n*n]: state-cost weighted matrix
 * @R[m*m]: input-cost weighted matrix
 */
void dare(double X[],
  const unsigned int n, const unsigned int m,
  const double A[], const double B[], const double Q[], const double R[])
{
  const unsigned int nm = n*m;
  const unsigned int nn = n*n;
  const unsigned int mm = m*m;
  memcpy(X, Q, n*n*sizeof(double));

  double diff = 0;
  unsigned int iter;

  for (iter=0; iter<HECTR_ITER_MAX; iter++) {
    /* A.T*X[n*n] = A.T[n*n] * X[n*n] */
    double ATX[n*n];
    memset(ATX, 0, sizeof(ATX));
    for (unsigned int i=0; i<nn; i++)
      for (unsigned int j=0; j<n; j++)
        ATX[(i/n)*n+i%n]+=A[j*n+i/n]*X[j*n+i%n];

    /* A.T*X*A[n*n] = A.T*X[n*n] * A[n*n] */
    double ATXA[n*n];
    memset(ATXA, 0, sizeof(ATXA));
    for (unsigned int i=0; i<nn; i++)
      for (unsigned int j=0; j<n; j++)
        ATXA[(i/n)*n+i%n]+=ATX[(i/n)*n+j]*A[j*n+i%n];

    /* A.T*X*B[n*m] = A.T*X[n*n] * B[n*m] */
    double ATXB[n*n];
    memset(ATXB, 0, sizeof(ATXB));
    for (unsigned int i=0; i<nm; i++)
      for (unsigned int j=0; j<n; j++)
        ATXB[(i/m)*m+i%m]+=ATX[(i/m)*n+j]*B[j*m+i%m];

    /* B.T*X[m*n] = B.T[m*n] * X[n*n] */
    double BTX[m*n];
    memset(BTX, 0, sizeof(BTX));
    for (unsigned int i=0; i<nm; i++)
      for (unsigned int j=0; j<n; j++)
        BTX[(i/n)*n+i%n]+=B[j*m+i/n]*X[j*n+i%n];

    /* B.T*X*A[m*n] = B.T*X[m*n] * A[n*n] */
    double BTXA[m*n];
    memset(BTXA, 0, sizeof(BTXA));
    for (unsigned int i=0; i<nm; i++)
      for (unsigned int j=0; j<n; j++)
        BTXA[(i/n)*n+i%n]+=BTX[(i/n)*n+j]*A[j*n+i%n];

    /* B.T*X*B[m*m] = B.T*X[m*n] * B[n*m] */
    double BTXB[m*m];
    memset(BTXB, 0, sizeof(BTXB));
    for (unsigned int i=0; i<mm; i++)
      for (unsigned int j=0; j<n; j++)
        BTXB[(i/m)*m+i%m]+=BTX[(i/m)*n+j]*B[j*m+i%m];

    /* inv(B.T*X*B + R)[m*m] */
    double RpBTXBinv[m*m];
    memset(RpBTXBinv, 0, sizeof(RpBTXBinv));
    for (unsigned int i=0; i<m*m; i++)
      RpBTXBinv[i] = R[i]+BTXB[i];
    dgeinv(RpBTXBinv, m, RpBTXBinv);

    /* A.T*X*B[n*m] * inv(B.T*X*B + R)[m*m] */
    double ATXBmRpBTXBinv[n*m];
    memset(ATXBmRpBTXBinv, 0, sizeof(ATXBmRpBTXBinv));
    for (unsigned int i=0; i<nm; i++)
      for (unsigned int j=0; j<m; j++)
        ATXBmRpBTXBinv[(i/m)*m+i%m]+=ATXB[(i/m)*m+j]*RpBTXBinv[j*m+i%m];

    /* A.T*X*B*inv(B.T*X*B + R)[n*m] * B.T*X*A[m*n] */
    double ATXBmRpBTXBinvBTXA[n*n];
    memset(ATXBmRpBTXBinvBTXA, 0, sizeof(ATXBmRpBTXBinvBTXA));
    for (unsigned int i=0; i<nn; i++)
      for (unsigned int j=0; j<m; j++)
        ATXBmRpBTXBinvBTXA[(i/n)*n+i%n]+=ATXBmRpBTXBinv[(i/n)*m+j]*BTXA[j*n+i%n];

    /* Xnext = A.T*X*A - A.T*X*B*inv(R + BT*X * B)*BT*X*A + Q; */
    double Xnext[n*n];
    memset(Xnext, 0, sizeof(Xnext));
    for (unsigned int i=0; i<n*n; i++)
      Xnext[i] = ATXA[i] - ATXBmRpBTXBinvBTXA[i] + Q[i];

    diff = 0;
    for (unsigned int i=0; i<n*n; i++)
      if (diff < __builtin_fabs(Xnext[i]-X[i]))
        diff = __builtin_fabs(Xnext[i]-X[i]);

    memcpy(X, Xnext, sizeof(Xnext));

    if (diff < HECTR_TOLERANCE)
      break;
  }

  if (diff < HECTR_TOLERANCE)
    ;//printf("iteration = %i, diff = %e.\n", iter, diff);
  else
    printf("tolerance %e not reach, iterated %i, diff = %e.\n", HECTR_TOLERANCE, iter, diff);
}

END_DECLS
