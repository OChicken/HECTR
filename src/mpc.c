/*
 * Unconstrained linear model predictive control.
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
#include <errno.h>

BEGIN_DECLS

static void calc_horizon_matrices(
  double AA[], double BB[], double Theta[], double CC[], double QQ[], double RR[],
  const unsigned int l, const unsigned int n, const unsigned int m, const unsigned int N,
  const double A[], const double B[], const double C[], const double Q[], const double R[])
{
  const unsigned int ll = l*l;
  const unsigned int nn = n*n;
  const unsigned int nm = n*m;
  const unsigned int mm = m*m;
  const unsigned int ln = l*n;
  const unsigned int mN = m*N;
  const unsigned int nNp1 = n*(N+1);
  const unsigned int lNp1 = l*(N+1);
  double An[nn]; /* A^n, n=0,1,...,N */
  memset(An, 0, sizeof(An));
  double AnB[nNp1*m]; /* A^n*B, n=-1,0,1,...,N-1 */
  memset(AnB, 0, sizeof(AnB));
  memset(AA, 0, nNp1 * n   * sizeof(double));
  memset(BB, 0, nNp1 * m   * sizeof(double)); /* sum_{j=-1}^{i-1} A^j*B, i=0,...,N */
  memset(Theta, 0, nNp1 * mN * sizeof(double));
  memset(CC, 0, lNp1 * nNp1 * sizeof(double));
  memset(QQ, 0, lNp1 * lNp1 * sizeof(double));
  memset(RR, 0, mN * mN * sizeof(double));
  /* k=0 */
  for (unsigned int i=0; i<nn; i+=n)
    An[i++]=1;
  for (unsigned int j=0, all=nNp1*n; j<all; j++)
    AA[(j/n)*n+j%n]=An[j];
  for (unsigned int j=0; j<ll; j++)
    QQ[(j/l)*lNp1 + j%l] = Q[j];
  for (unsigned int j=0; j<mm; j++)
    RR[(j/m)*mN + j%m] = R[j];
  for (unsigned int j=0; j<ln; j++)
    CC[(j/n)*lNp1 + j%n] = C[j];
  /* k = 1 ~ N+1 */
  for (unsigned int k=1; k<N+1; k++) {
    /* set AnB */
    double tmpAnB[nm];
    memset(tmpAnB, 0, sizeof(tmpAnB));
    for (unsigned int i=0; i<nm; i++)
      for (unsigned int j=0; j<n; j++)
        tmpAnB[i]+=An[(i/m)*n+j]*B[j*m+i%m];
    memcpy(&AnB[k*nm], tmpAnB, sizeof(tmpAnB));
    /* set BB */
    for (unsigned int i=0; i<nm; i++)
      BB[k*nm+i] = BB[(k-1)*nm+i] + tmpAnB[i];
    /* set AA */
    double tmpAn[nn]; /* n*n */
    memset(tmpAn, 0, sizeof(tmpAn));
    for (unsigned int i=0; i<nn; i++)
      for (unsigned int j=0; j<n; j++)
        tmpAn[i]+=An[(i/n)*n+j]*A[j*n+i%n];
    memcpy(An, tmpAn, sizeof(tmpAn));
    memcpy(&AA[k*nn], An, sizeof(An));
    /* set Theta */
    for (unsigned int i=k; i<N+1; i++)
      for (unsigned int j=0; j<nm; j++)
        Theta[(i*n+j/m)*mN+(i-k)*m+j%m]=BB[k*nm+(j/m)*m+j%m];
    /* set QQ and RR */
    for (unsigned int j=0; j<ll; j++)
      QQ[(k*l+j/l)*lNp1 + k*l+j%l] = Q[j];
    if (k<N)
      for (unsigned int j=0; j<mm; j++)
        RR[(k*m+j/m)*mN + k*m+j%m]=R[j];
    /* set CC */
    for (unsigned int j=0; j<ln; j++)
      CC[(k*l+j/n)*nNp1 + k*n+j%n] = C[j];
  }
}

/**
 * calc_ef - calculate e(.|0) and f(.|0)
 * 
 * @e[l(N+1)]: e(.|0) = r(.|0) - CC*f(.|0)
 * @f[n(N+1)]: f(.|0) = AA*x0 + BB*u0
 * @l: dimension of r
 * @n: dimension of x
 * @m: dimension of u
 * @N: horizons
 * @AA[n(N+1)*n]: horizon extended A
 * @BB[n(N+1)*m]: sum_{j=-1}^{i-1} A^j*B, i=0,...,N
 * @CC[l(N+1)*n(N+1)]: diag(C,...,C)
 * @r[l(N+1)]: reference trajectory
 * @x0[n]: initial state
 * @ur[m]: initial control (more accurate, is u_{-1})
 */
static void calc_ef(double e[], double f[],
  const unsigned int l, const unsigned int n, const unsigned int m, const unsigned int N,
  const double AA[], const double BB[], const double CC[],
  const double xhat[], const double uhat[], const double xr[], const double ur[])
{
  const unsigned int nNp1 = n*(N+1);
  const unsigned int lNp1 = l*(N+1);
  memset(e, 0, lNp1*sizeof(double));
  memset(f, 0, nNp1*sizeof(double));
  double AAxhat[nNp1];
  double BBuhat[nNp1];
  double AAx[nNp1];
  double BBu[nNp1];
  memset(AAxhat, 0, sizeof(AAxhat));
  memset(BBuhat, 0, sizeof(BBuhat));
  memset(AAx,0,sizeof(AAx));
  memset(BBu,0,sizeof(BBu));
  for (unsigned int i=0; i<nNp1; i++){
    for (unsigned int j=0; j<n; j++)
      AAxhat[i] += AA[i*n+j]*xhat[j];
    for (unsigned int j=0; j<n; j++)
      AAx[i]    += AA[i*n+j]*(xhat[j]-xr[j]);
    for (unsigned int j=0; j<m; j++)
      BBuhat[i] += BB[i*m+j]*uhat[j];
    for (unsigned int j=0; j<m; j++)
      BBu[i]    += BB[i*m+j]*(uhat[j]-ur[j]);
    f[i] = AAxhat[i]+BBuhat[i];
  }
  for (unsigned int i=0; i<lNp1; i++)
    for (unsigned int j=0; j<nNp1; j++)
      e[i] += CC[i*nNp1+j]*(AAx[j]+BBu[j]);
}

/**
 * calc_Hc - calculate H, c, K
 * 
 * @H[mN*mN]: = Theta.T * CC.T * QQ * CC * Theta + RR
 * @c[mN]: = -Theta.T * CC.T * QQ * e
 * @K[mN*l(N+1)]: = inv(H) * Theta.T * CC.T * QQ
 * @l: dimension of r
 * @n: dimension of x
 * @m: dimension of u
 * @N: horizons
 * @Theta[n(N+1)*mN]: coefficient matrix for du[mN]
 * @CC[l(N+1)*n(N+1)]: diag(C,...,C)
 * @QQ[m(N+1)*m(N+1)]: diag(Q,...,Q)
 * @RR[mN*mN]: diag(R,...,R)
 */
static void calc_Hc(double H[], double c[],
  const unsigned int l, const unsigned int n, const unsigned int m, const unsigned int N,
  const double Theta[], const double CC[], const double QQ[], const double RR[], const double e[])
{
  const unsigned int mN = m*N;
  const unsigned int nNp1 = n*(N+1);
  const unsigned int lNp1 = l*(N+1);

  /* CC*Theta[l(N+1)*mN] = CC[l(N+1)*n(N+1)] * Theta[n(N+1)*mN] */
  double CCTheta[lNp1*mN];
  memset(CCTheta, 0, sizeof(CCTheta));
  for (unsigned int i=0, all=lNp1*mN; i<all; i++)
    for (unsigned int j=0; j<nNp1; j++)
      CCTheta[i]+=CC[(i/mN)*nNp1+j]*Theta[j*mN+i%mN];

  /**  Theta.T*CC.T*QQ[mN*l(N+1)]
   * = (CC*Theta).T[mN*l(N+1)] * QQ[l(N+1)*l(N+1)] */
  double ThetaTCCTQQ[mN*lNp1];
  memset(ThetaTCCTQQ, 0, sizeof(ThetaTCCTQQ));
  for (unsigned int i=0, all=mN*lNp1; i<all; i++)
    for (unsigned int j=0; j<lNp1; j++)
      ThetaTCCTQQ[i]+=CCTheta[j*mN+i/lNp1]*QQ[j*lNp1+i%lNp1];

  /**  H[mN*mN] = (Theta.T*CC.T*QQ*CC*Theta+RR)[mN*mN]
   * = (Theta.T*CC.T*QQ)[mN*l(N+1)] * (CC*Theta)[l(N+1)*mN] + RR[mN*mN] */
  memcpy(H, RR, mN * mN * sizeof(double));
  for (unsigned int i=0, all=mN*mN; i<all; i++)
    for (unsigned int j=0; j<lNp1; j++)
      H[i] += ThetaTCCTQQ[(i/mN)*lNp1+j]*CCTheta[j*mN+i%mN];

  /* c[mN] = (Theta.T*CC.T*QQ)[mN*l(N+1)] * e[l(N+1)] */
  memset(c, 0, mN*sizeof(double));
  for (unsigned int i=0; i<mN; i++)
    for (unsigned int j=0; j<lNp1; j++)
      c[i] += ThetaTCCTQQ[i*lNp1+j]*e[j];
}

static unsigned int calc_bnddim(
  const unsigned int n, const unsigned int m, const unsigned int N,
  const double dumin[], const double dumax[],
  const double  umin[], const double  umax[],
  const double  xmin[], const double  xmax[])
{
  const unsigned int nNp1 = n*(N+1);
  const unsigned int mN = m*N;
  unsigned int dim = 0;
  if (((dumin)&&(!dumax)) || ((!dumin)&&(dumax))) {
    errno = EINVAL;
    fprintf(stderr, "\033[1m\033[31merror:\033[0m \033[1m%s\033[0m. "
      "dumin and dumax should be set in pairs.\n", strerror(errno));
    abort();
  }
  else if (dumin&&dumax)
    dim+=2*mN;
  if (((umin)&&(!umax)) || ((!umin)&&(umax))) {
    errno = EINVAL;
    fprintf(stderr, "\033[1m\033[31merror:\033[0m \033[1m%s\033[0m. "
      "umin and umax should be set in pairs.\n", strerror(errno));
    abort();
  }
  else if (umin&&umax)
    dim+=2*mN;
  if (((xmin)&&(!xmax)) || ((!xmin)&&(xmax))) {
    errno = EINVAL;
    fprintf(stderr, "\033[1m\033[31merror:\033[0m \033[1m%s\033[0m. "
      "xmin and xmax should be set in pairs.\n", strerror(errno));
    abort();
  }
  else if (xmin&&xmax)
    dim+=2*nNp1;
  return dim;
}

static void calc_bnd_du(double A[], double b[],
  const unsigned int m, const unsigned int N,
  const double dumin[], const double dumax[])
{
  const unsigned int mN = m*N;
  double nIx[mN*mN], pIx[mN*mN];
  memset(nIx, 0, sizeof(nIx));
  memset(pIx, 0, sizeof(pIx));
  for (unsigned int k=0; k<N; k++) {
    for (unsigned int i=0; i<m; i++) {
      nIx[(k*m+i)*N + k*m+i]=-1;
      pIx[(k*m+i)*N + k*m+i]=+1;
      b[m* k   +i]=+dumin[i];
      b[m*(k+N)+i]=-dumax[i];
    }
  }
  memcpy(A      , nIx, sizeof(nIx));
  memcpy(A+mN*mN, pIx, sizeof(pIx));
}

static void calc_bnd_u(double A[], double b[],
  const unsigned int m, const unsigned int N,
  const double umin[], const double umax[],
  const double u0[])
{
  const unsigned int mN = m*N;
  double nIx[mN*mN], pIx[mN*mN];
  memset(nIx, 0, sizeof(nIx));
  memset(pIx, 0, sizeof(pIx));
  for (unsigned int k=0; k<N; k++) {
    for (unsigned int i=0; i<m; i++) {
      nIx[(k*m+i)*N + k*m+i]=-1;
      pIx[(k*m+i)*N + k*m+i]=+1;
      b[m* k   +i]=+umin[i]-u0[i];
      b[m*(k+N)+i]=-umax[i]+u0[i];
    }
  }
  memcpy(A      , nIx, sizeof(nIx));
  memcpy(A+mN*mN, pIx, sizeof(pIx));
}

/**
 * @A[n(N+1)*mN]:
 * @Theta[n(N+1)*mN]: coefficient matrix for du[mN]
 */
static void calc_bnd_x(double A[], double b[],
  const unsigned int n, const unsigned int m, const unsigned int N,
  const double xmin[], const double xmax[],
  const double Theta[], const double f[])
{
  const unsigned int nm = n*m;
  const unsigned int mN = m*N;
  const unsigned int nNp1 = n*(N+1);
  double nTheta[nNp1*mN];
  memset(nTheta, 0, sizeof(nTheta));
  double pTheta[nNp1*mN];
  memcpy(pTheta, Theta, sizeof(pTheta));
  for (unsigned int k=0; k<N+1; k++) {
    for (unsigned int i=k; i<N; i++)
      for (unsigned int j=0; j<nm; j++)
        nTheta[((i+1)*n+j/m)*mN+(i-k)*m+j%m] = -pTheta[((i+1)*n+j/m)*mN+(i-k)*m+j%m];
    for (unsigned int i=0; i<n; i++) {
      b[n* k     +i]=+xmin[i]-f[k*n+i];
      b[n*(k+N+1)+i]=-xmax[i]+f[k*n+i];
    }
  }
  memset(A, 0, 2*nNp1*mN*sizeof(double));
  memcpy(A        , nTheta, sizeof(nTheta));
  memcpy(A+nNp1*mN, pTheta, sizeof(pTheta));
}

static void calc_bnd(double A[], double b[],
  const unsigned int n, const unsigned int m, const unsigned int N,
  const unsigned int bnddim,
  const double u0[], const double Theta[], const double f[],
  const double dumin[], const double dumax[],
  const double  umin[], const double  umax[],
  const double  xmin[], const double  xmax[])
{
  const unsigned int mN = m*N;
  const unsigned int nNp1 = n*(N+1);
  if (bnddim==2*mN) {
    if ((dumin)&&(dumax))
      calc_bnd_du(A, b, m, N, dumin, dumax);
    if ((umin)&&(umax))
      calc_bnd_u(A, b, m, N, umin, umax, u0);
  }
  else if (bnddim==2*nNp1) {
    calc_bnd_x(A, b, n, m, N, xmin, xmax, Theta, f);
  }
  else if (bnddim==4*mN) {
    calc_bnd_du(A, b, m, N, dumin, dumax);
    calc_bnd_u(A+2*mN*mN, b+2*mN, m, N, umin, umax, u0);
  }
  else if (bnddim==2*mN+2*nNp1) {
    if ((dumin)&&(dumax))
      calc_bnd_du(A, b, m, N, dumin, dumax);
    if ((umin)&&(umax))
      calc_bnd_u(A, b, m, N, umin, umax, u0);
    calc_bnd_x(A+2*mN*mN, b+2*mN, n, m, N, xmin, xmax, Theta, f);
  }
  else if (bnddim==4*mN+2*nNp1) {
    calc_bnd_du(A, b, m, N, dumin, dumax);
    calc_bnd_u(A+2*mN*mN, b+2*mN, m, N, umin, umax, u0);
    calc_bnd_x(A+4*mN*mN, b+4*mN, n, m, N, xmin, xmax, Theta, f);
  }
  else {
    assert(bnddim==0);
    A=NULL, b=NULL;
  }
}

static void calc_u(double u[],
  const unsigned int m, const unsigned int N,
  const double du[], const double ur[])
{
  const unsigned int mN = m*N;
  double sum_du[mN];
  memcpy(sum_du, du, mN*sizeof(double));
  for (unsigned int k=0; k<N; k++) {
    for (unsigned int i=0; i<m; i++) {
      if (k>=1)
        sum_du[k*m+i] += sum_du[(k-1)*m+i];
      u[k*m+i]=ur[i]+sum_du[k*m+i];
    }
  }
}

/**
 * mpc - Model predictive controller for discrete-time systems.
 *
 * @x[n(N+1)]: predicted state
 * @u[mN]: optimal control
 * @l: dimension of r
 * @n: dimension of x
 * @m: dimension of u
 * @N: horizons
 * @A[n*n]: state matrix
 * @B[n*m]: input matrix
 * @C[l*n]: output matrix
 * @Q[n*n]: state-cost weighted matrix
 * @R[m*m]: input-cost weighted matrix
 * @r[l(N+1)]: reference trajectory
 * @xr[n]: initial state (more accurate, is x(k))
 * @ur[m]: initial control (more accurate, is u(k-1))
 */
void ctr_mpc(double u[],
  const unsigned int l, const unsigned int n, const unsigned int m, const unsigned int N,
  const double A[], const double B[], const double C[], const double Q[], const double R[],
  const double xhat[], const double uhat[], const double xr[], const double ur[],
  const double dumin[], const double dumax[],
  const double  umin[], const double  umax[],
  const double  xmin[], const double  xmax[])
{
  const unsigned int mN = m*N;
  const unsigned int nNp1 = n*(N+1);
  const unsigned int lNp1 = l*(N+1);
  double AA[nNp1 * n];
  double BB[nNp1 * m];
  double Theta[nNp1 * mN];
  double CC[lNp1 * nNp1];
  double QQ[lNp1 * lNp1];
  double RR[mN * mN];
  calc_horizon_matrices(AA, BB, Theta, CC, QQ, RR, l, n, m, N, A, B, C, Q, R);
  double e[lNp1*1];
  double f[nNp1*1];
  calc_ef(e, f, l, n, m, N, AA, BB, CC, xhat, uhat, xr, ur);
  double H[mN*mN];
  double c[mN*1];
  calc_Hc(H, c, l, n, m, N, Theta, CC, QQ, RR, e);
  /* constraints */
  unsigned int bnddim = calc_bnddim(n, m, N, dumin, dumax, umin, umax, xmin, xmax);
  double Ain[bnddim*mN], bin[bnddim];
  calc_bnd(Ain, bin, n, m, N, bnddim, uhat, Theta, f, dumin, dumax, umin, umax, xmin, xmax);
  double du[mN];
  memset(du, 0, sizeof(du));
  if (bnddim)
    quadprog(du, mN, bnddim, 0, H, c, Ain, bin, NULL, NULL, NULL, NULL);
  else {
    double Hinv[mN*mN];
    dgeinv(Hinv, mN, H);
    for (unsigned int i=0; i<mN; i++)
      for (unsigned int j=0; j<mN; j++)
        du[i] += Hinv[i*mN+j]*(-c[j]);
  }
  calc_u(u, m, N, du, uhat);
}

#if 0
  printf("f=\n");
  for (unsigned int i=0; i<nNp1; i++)
    printf("%3u %6.2f \n", i, f[i]);
  printf("e=\n");
  for (unsigned int i=0; i<lNp1; i++)
    printf("%3u %6.2f \n", i, e[i]);

  printf("BB=\n");
  for (unsigned int i=0, all=nNp1*m; i<all; i++) {
    printf("%9.6f ", BB[i]);
    if ((i+1)%m==0)
      printf("\n");
  }
  printf("bin=\n");
  for (unsigned int i=0; i<bnddim; i++)
    printf("%3u %6.2f \n", i, bin[i]);
  printf("Ain=\n");
  for (unsigned int i=0; i<bnddim; i++) {
    for (unsigned int j=0; j<mN; j++)
      printf("%8.4f ", Ain[i*mN+j]);
    printf("\n");
  }
  printf("Theta=\n");
  for (unsigned int i=0, all=nNp1*mN; i<all; i++) {
    printf("%9.6f ", Theta[i]);
    if ((i+1)%mN==0)
      printf("\n");
  }
  printf("QQ=\n");
  for (unsigned int i=0, all=lNp1*lNp1; i<all; i++) {
    printf("%.f ", QQ[i]);
    if ((i+1)%lNp1==0)
      printf("\n");
  }
  printf("RR=\n");
  for (unsigned int i=0, all=mN*mN; i<all; i++) {
    printf("%9.6f ", RR[i]);
    if ((i+1)%mN==0)
      printf("\n");
  }
  printf("Delta u=\n");
  for (unsigned int i=0; i<mN; i++)
    printf("%9.6f \n", du[i]);
  printf("f(.|0)=\n");
  for (unsigned int i=0; i<lNp1; i++)
    printf("%9.6f \n", f[i]);
  printf("e(.|0)=\n");
  for (unsigned int i=0; i<lNp1; i++)
    printf("%9.6f \n", e[i]);
  printf("H=\n");
  for (unsigned int i=0, all=mN*mN; i<all; i++) {
    printf("%9.4f ", H[i]);
    if ((i+1)%mN==0)
      printf("\n");
  }
  printf("c=\n");
  for (unsigned int i=0; i<mN; i++)
    printf("%9.6f \n", c[i]);
#endif

END_DECLS
