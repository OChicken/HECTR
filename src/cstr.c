/*
 * CSTR model.
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

/* Parameters: */
static double rho = 1000;  /* Density of A-B Mixture (kg/m^3) */
static double cp = 0.239;   /* Heat capacity of A-B Mixture (kJ/kg K) */
static double DeltaH = -5e4; /* Heat of reaction for A->B (kJ/mol) */
/* E - Activation energy in the Arrhenius Equation (J/mol)
 * R - Universal Gas Constant = 8.31451 J/mol-K */
static double EoverR = 8750; /* The activation energy E devided by the universal ideal gas constant R (K) */
static double k0 = 7.2e10; /* Pre-exponential factor (min^-1) */
/* U - Overall Heat Transfer Coefficient (W/m^2-K)
 * A - Area - this value is specific for the U calculation (m^2) */
static double U = 54.94;
static double c0 = 1;   /* Feed Concentration (kmol/m^3) */
static double T0 = 350; /* Feed Temperature (K) */
static double r  = 0.219; /* radius of the container (m) */

/**
 * CSTR model from
 * 
 * Michael A. Henson and Dale E. Seborg.  Nonlinear Process Control.
 *      Prentice Hall PTR, Upper Saddle River, New Jersey, 1997.
 * 
 * Description:
 * Continuously Stirred Tank Reactor with energy balance and reaction A->B.
 * The temperature of the cooling jacket is the control.
 */
void cstr_ode(double xdot[], const double x[], const double u[], const double p[])
{
  double c = x[0];  /* Concentration of A in CSTR (mol/m^3) */
  double T = x[1];  /* Temperature in CSTR (K) */
  double h = x[2];
  double Tc = u[0]; /* Temperature of cooling jacket (K) */
  double F  = u[1];
  double F0 = p[0];

  double kT = k0*exp(-EoverR/T);
  double S = M_PI*r*r;

  xdot[0] = F0*(c0 - c)/(S*h) - kT*c;
  xdot[1] = F0*(T0 - T)/(S*h) + -DeltaH/(rho*cp)*kT*c + 2*U/(r*rho*cp)*(Tc-T);
  xdot[2] = (F0 - F)/S;
}

void cstr_jacobian(double J[], const double x[],
  __attribute__((unused))const double u[], const double p[])
{
  double c = x[0];
  double T = x[1];
  double h = x[2];
  double F0 = p[0];

  double kT = k0*exp(-EoverR/T);
  double S = M_PI*r*r;

  memcpy(J, (double[]){
    -F0/(S*h) - kT,
    -kT*EoverR/(T*T)*c,
    -F0*(c0-c)/(S*h*h),
    (-DeltaH)/(rho*cp)*kT,
    -F0/(S*h) + (-DeltaH)/(rho*cp)*kT*EoverR/(T*T)*c + -2*U/(r*rho*cp),
    -F0*(T0-T)/(S*h*h),
    0,0,0
  }, 3*3*sizeof(double));
}

void cstr_linearize(double A[], double B[], double Bp[],
  const unsigned int nx, const unsigned int nu, const unsigned int np,
  const double xs[], const double us[], const double ps[], const double dt)
{
  /* init */
  double c = xs[0];
  double T = xs[1];
  double h = xs[2];

  double jacA [nx*nx];
  double jacB [nx*nu];
  double jacBp[nx*np];

  double S = M_PI*r*r;

  /* jacobian A, B, Bp */
  cstr_jacobian(jacA, xs, us, ps);
  memcpy(jacB, (double[]){
    0             ,0,
    2*U/(r*rho*cp),0,
    0             ,-1/S
  }, sizeof(jacB));
  memcpy(jacBp, (double[]){
    (c0-c)/(S*h),
    (T0-T)/(S*h),
    1/S
  }, sizeof(jacBp));

  double expmA[nx*nx], expmB[nx*nx];
  ctr_c2d(expmA, expmB, nx, jacA, dt);
  memcpy(A, expmA, sizeof(expmA));

  /* B[nx*nu] = expmB[nx*nx]*jacB[nx*nu] */
  memset(B, 0, nx*nu*sizeof(double));
  for (unsigned int i=0, all=nx*nu; i<all; i++)
    for (unsigned int j=0; j<nx; j++)
      B[i] += expmB[(i/nu)*nx+j]*jacB[j*nu+i%nu];

  /* Bp[nx*np] = expmB[nx*nx]*jacB[nx*np] */
  memset(Bp, 0, nx*np*sizeof(double));
  for (unsigned int i=0, all=nx*np; i<all; i++)
    for (unsigned int j=0; j<nx; j++)
      Bp[i] += expmB[(i/np)*nx+j]*jacBp[j*np+i%np];
}

END_DECLS
