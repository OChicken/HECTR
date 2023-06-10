/*
 * Continuous/Discrete Riccati equation solvers.
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
#include <float.h>

BEGIN_DECLS

double eps(double a)
{
  a = __builtin_fabs(a);
  double r;
  unsigned long tmp;
  //unsigned char buf[sizeof(double)];
  memcpy(&tmp, &a, sizeof(double));
  tmp += 1;
  memcpy(&r, &tmp, sizeof(double));
  return r-a;
}

int dgeinv(double Ainv[], const unsigned int n, const double A[])
{
  int err = 0;
  int ipiv[n];
  memcpy(Ainv, A, n*n*sizeof(double));
  err |= LAPACKE_dgetrf(LAPACK_ROW_MAJOR, n, n, Ainv, n, ipiv);
  err |= LAPACKE_dgetri(LAPACK_ROW_MAJOR, n, Ainv, n, ipiv);
  return err;
}

int zgeinv(_Complex double Ainv[], const unsigned int n, const _Complex double A[])
{
  int err = 0;
  int ipiv[n];
  memcpy(Ainv, A, n*n*sizeof(_Complex double));
  err |= LAPACKE_zgetrf(LAPACK_ROW_MAJOR, n, n, Ainv, n, ipiv);
  err |= LAPACKE_zgetri(LAPACK_ROW_MAJOR, n, Ainv, n, ipiv);
  return err;
}

int dgepinv(double Ainv[], const unsigned int m, const unsigned int n, const double A[])
{
  double a[m*n];
  memcpy(a, A, sizeof(a));

  unsigned int k = (m<n)? m:n;
  double s[k];
  memset(s, 0, sizeof(s));
  double superb[k];
  memset(superb, 0, sizeof(superb));

  double u[m*m];
  double vt[n*n];
  int err = LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', m, n, a,n, s, u,m, vt,n, superb);

  double eps_s0 = ((m>n)? m:n)* eps(s[0]);
  for (unsigned int i=0; i<k; i++)
    if (s[i]<eps_s0)
      {k=i; break;}

  if (m>n)
    for (unsigned int i=0; i<n*n; i++)
      vt[i] *= 1/s[i/n];
  else
    for (unsigned int i=0; i<m*m; i++)
      u[i] *= 1/s[i%m];

  memset(Ainv, 0, n*m*sizeof(double));
  for (unsigned int i=0, all=n*m; i<all; i++)
    for (unsigned int j=0; j<k; j++)
      Ainv[i] += vt[j*n+(i/m)] * u[(i%m)*m+j];

  return err;
}

int dexpm(double expA[], const unsigned int n, const double A[])
{
  unsigned int nn = n*n;
  _Complex double zA[nn], D[n], V[nn];
  for (unsigned int i=0; i<nn; i++)
    zA[i] = A[i];
  int err = LAPACKE_zgeev(LAPACK_ROW_MAJOR, 'N', 'V', n, zA, n, D, NULL, n, V, n);
  if (err)
    return err;

  _Complex double expD[n];
  for (unsigned int i=0; i<n; i++)
    expD[i] = cexp(D[i]);
  _Complex double VD[nn];
  for (unsigned int i=0; i<n; i++)
    for (unsigned int j=0; j<n; j++)
      VD[i*n+j] = V[i*n+j]*expD[j];
  _Complex double Vinv[nn];
  zgeinv(Vinv, n, V);

  _Complex double zexpA[nn];
  memset(zexpA, 0, sizeof(zexpA));
  for (unsigned int i=0; i<nn; i++) {
    for (unsigned int j=0; j<n; j++)
      zexpA[i] += VD[(i/n)*n+j]*Vinv[j*n+i%n];
    //assert(fabs(cimag(zexpA[i]))<FLT_EPSILON);
    expA[i] = creal(zexpA[i]);
  }
  return err;
}

void d2z_vector(_Complex double vz[],
  const unsigned int n, const unsigned int slots,
  const double vd[])
{
  memset(vz, 0, slots*sizeof(_Complex double));
  for (unsigned int i=0; i<n; i++)
    vz[i] = vd[i];
}

void d2z_matrix(_Complex double Az[],
  const unsigned int row, const unsigned int col, const unsigned int slots,
  const double Ad[])
{
  memset(Az, 0, slots*slots*sizeof(_Complex double));
  for (unsigned int i=0; i<row; i++)
    for (unsigned int j=0; j<col; j++)
      Az[i*slots+j]=Ad[i*col+j];
}

END_DECLS
