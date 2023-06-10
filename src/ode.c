/*
 * ode45 solver.
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

void ode45(double xnext[],
  const unsigned int nx,
  const double x[], const double u[], const double p[],
  void (*const ode)(double[], const double[], const double[], const double[]),
  const double dt)
{
  double k1[nx];
  ode(k1, x, u, p);
  for (unsigned int i=0; i<nx; i++)
    k1[i] *= dt;

  double k2[nx];
  double x_p_k1_d_2[nx]; /* x+k1/2 */
  for (unsigned int i=0; i<nx; i++)
    x_p_k1_d_2[i] = x[i] + k1[i]/2;
  ode(k2, x_p_k1_d_2, u, p);
  for (unsigned int i=0; i<nx; i++)
    k2[i] *= dt;

  double k3[nx];
  double x_p_k2_d_2[nx]; /* x+k2/2 */
  for (unsigned int i=0; i<nx; i++)
    x_p_k2_d_2[i] = x[i] + k2[i]/2;
  ode(k3, x_p_k2_d_2, u, p);
  for (unsigned int i=0; i<nx; i++)
    k3[i] *= dt;

  double k4[nx];
  double x_p_k3[nx]; /* x+k3 */
  for (unsigned int i=0; i<nx; i++)
    x_p_k3[i] = x[i] + k3[i];
  ode(k4, x_p_k3, u, p);
  for (unsigned int i=0; i<nx; i++)
    k4[i] *= dt;

  /* k1+2k2+2k3+k4 */
  for (unsigned int i=0; i<nx; i++)
    xnext[i] = x[i] + (k1[i]+2*k2[i]+2*k3[i]+k4[i])/6;
}

void ode15s(double xnext[],
  const unsigned int n,
  const double x[], const double u[], const double p[],
  void (*const ode     )(double[], const double[], const double[], const double[]),
  void (*const jacobian)(double[], const double[], const double[], const double[]),
  const double dt)
{
  const unsigned int nn = n*n;
  double J[nn];
  jacobian(J, x, u, p);
  double C[nn];
  for (unsigned int i=0; i<nn; i++) {
    C[i] = -J[i]*dt;
    if ((i/n)==(i%n))
      C[i] += 1;
  }
  double Cinv[nn];
  dgeinv(Cinv, n, C);
  double fx[n];
  ode(fx, x, u, p);

  memcpy(xnext, x, n*sizeof(double));
  double diff[n];
  memset(diff, 0, sizeof(diff));
  for (unsigned int i=0; i<n; i++) {
    for (unsigned int j=0; j<n; j++)
      diff[i] += Cinv[i*n+j]*fx[j];
    diff[i] *= dt;
    xnext[i] += diff[i];
  }
}

END_DECLS
