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
 * calc_Hc - calculate H, c, K
 * 
 * @H[mN*mN]: = Theta.T * CC.T * QQ * CC * Theta + RR
 * @c[mN]: = -Theta.T * CC.T * QQ * e
 * @K[mN*l(N+1)]: = inv(H) * Theta.T * CC.T * QQ
 * @l: dimension of r
 * @n: dimension of x
 * @m: dimension of u
 * @N: horizons
 * @AA[n(N+1)*n]: horizon extended A
 * @BB[n(N+1)*m]: sum_{j=-1}^{i-1} A^j*B, i=0,...,N
 * @Theta[n(N+1)*mN]: coefficient matrix for du[mN]
 * @CC[l(N+1)*n(N+1)]: diag(C,...,C)
 * @QQ[m(N+1)*m(N+1)]: diag(Q,...,Q)
 * @RR[mN*mN]: diag(R,...,R)
 * @xhat[l(N+1)]: reference trajectory
 * @x0[n]: initial state
 * @ur[m]: initial control (more accurate, is u_{-1})
 */
static void calc_coeff(_Complex double HinvThetaTCCTQQCCAAz[], _Complex double HinvThetaTCCTQQCCBBz[],
  const unsigned int l, const unsigned int n, const unsigned int m, const unsigned int N,
  const unsigned int slots,
  const double AA[], const double BB[],
  const double Theta[], const double CC[], const double QQ[], const double RR[])
{
  const unsigned int mN = m*N;
  const unsigned int nNp1 = n*(N+1);
  const unsigned int lNp1 = l*(N+1);

  /* (CC*AA)[l(N+1)*n] = CC[l(N+1)*n(N+1)]*AA[n(N+1)*n] */
  double CCAA[lNp1*n];
  memset(CCAA, 0, sizeof(CCAA));
  for (unsigned int i=0, all=lNp1*n; i<all; i++)
    for (unsigned int j=0; j<nNp1; j++)
      CCAA[i] += CC[(i/n)*nNp1+j]*AA[j*n+i%n];

  /* (CC*BB)[l(N+1)*m] = CC[l(N+1)*n(N+1)]*BB[n(N+1)*m] */
  double CCBB[lNp1*m];
  memset(CCBB, 0, sizeof(CCBB));
  for (unsigned int i=0, all=lNp1*m; i<all; i++)
    for (unsigned int j=0; j<nNp1; j++)
      CCBB[i] += CC[(i/m)*nNp1+j]*BB[j*m+i%m];

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

  /* (Theta.T*CC.T*QQ*CC*AA)[mN*n] = (Theta.T*CC.T*QQ)[mN*l(N+1)] * (CC*AA)[l(N+1)*n] */
  double ThetaTCCTQQCCAA[mN*n];
  memset(ThetaTCCTQQCCAA, 0, sizeof(ThetaTCCTQQCCAA));
  for (unsigned int i=0, all=mN*n; i<all; i++)
    for (unsigned int j=0; j<lNp1; j++)
      ThetaTCCTQQCCAA[i] += ThetaTCCTQQ[(i/n)*lNp1+j] * CCAA[j*n+i%n];

  /* (Theta.T*CC.T*QQ*CC*BB)[mN*m] = (Theta.T*CC.T*QQ)[mN*l(N+1)] * (CC*BB)[l(N+1)*m] */
  double ThetaTCCTQQCCBB[mN*m];
  memset(ThetaTCCTQQCCBB, 0, sizeof(ThetaTCCTQQCCBB));
  for (unsigned int i=0, all=mN*m; i<all; i++)
    for (unsigned int j=0; j<lNp1; j++)
      ThetaTCCTQQCCBB[i] += ThetaTCCTQQ[(i/m)*lNp1+j] * CCBB[j*m+i%m];

  /**  H[mN*mN] = (Theta.T*CC.T*QQ*CC*Theta+RR)[mN*mN]
   * = (Theta.T*CC.T*QQ)[mN*l(N+1)] * (CC*Theta)[l(N+1)*mN] + RR[mN*mN] */
  double H[mN*mN];
  memcpy(H, RR, mN * mN * sizeof(double));
  for (unsigned int i=0, all=mN*mN; i<all; i++)
    for (unsigned int j=0; j<lNp1; j++)
      H[i] += ThetaTCCTQQ[(i/mN)*lNp1+j]*CCTheta[j*mN+i%mN];

  double Hinv[mN*mN];
  dgeinv(Hinv, mN, H);

  /* (Hinv*Theta.T*CC.T*QQ*CC*AA)[mN*n] = Hinv[mN] * (Theta.T*CC.T*QQ*CC*AA)[mN*n] */
  double HinvThetaTCCTQQCCAA[mN*n];
  memset(HinvThetaTCCTQQCCAA, 0, mN*n*sizeof(double));
  for (unsigned int i=0, all=mN*n; i<all; i++)
    for (unsigned int j=0; j<mN; j++)
      HinvThetaTCCTQQCCAA[i] += Hinv[(i/n)*mN+j] * ThetaTCCTQQCCAA[j*n+i%n];
  d2z_matrix(HinvThetaTCCTQQCCAAz, mN, n, slots, HinvThetaTCCTQQCCAA);

  /* (Hinv*Theta.T*CC.T*QQ*CC*BB)[mN*m] = Hinv[mN] * (Theta.T*CC.T*QQ*CC*BB)[mN*m] */
  double HinvThetaTCCTQQCCBB[mN*m];
  memset(HinvThetaTCCTQQCCBB, 0, mN*m*sizeof(double));
  for (unsigned int i=0, all=mN*m; i<all; i++)
    for (unsigned int j=0; j<mN; j++)
      HinvThetaTCCTQQCCBB[i] += Hinv[(i/m)*mN+j] * ThetaTCCTQQCCBB[j*m+i%m];
  d2z_matrix(HinvThetaTCCTQQCCBBz, mN, m, slots, HinvThetaTCCTQQCCBB);
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
void ctr_hempc(he_ct_t *up,
  const unsigned int l, const unsigned int n, const unsigned int m, const unsigned int N,
  const unsigned int slots,
  const double A[], const double B[], const double C[], const double Q[], const double R[],
  const he_ct_t *xhat, const he_ct_t *uhat, const he_ct_t *xr, const he_ct_t *ur,
  const he_evk_t *rk)
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
  _Complex double BBz[slots*slots];
  d2z_matrix(BBz, nNp1, m, slots, BB);

  _Complex double HinvThetaTCCTQQCCAA[slots*slots];
  _Complex double HinvThetaTCCTQQCCBB[slots*slots];
  calc_coeff(HinvThetaTCCTQQCCAA, HinvThetaTCCTQQCCBB, l, n, m, N, slots, AA, BB, Theta, CC, QQ, RR);

  he_ct_t xdiff;
  he_ct_t udiff;
  he_ct_t du;
  he_ct_t uhat_copy;
  he_ct_t HinvThetaTCCTQQCCAAx;
  he_ct_t HinvThetaTCCTQQCCBBu;
  he_alloc_ct(&xdiff);
  he_alloc_ct(&udiff);
  he_alloc_ct(&du);
  he_alloc_ct(&uhat_copy);
  he_alloc_ct(&HinvThetaTCCTQQCCAAx);
  he_alloc_ct(&HinvThetaTCCTQQCCBBu);
  /* xhat-xr */
  he_sub(&xdiff, xhat, xr);
  /* uhat-ur */
  he_sub(&udiff, uhat, ur);
  /* Hinv*Theta.T*CC.T*QQ*CC*AA*(xhat-xr) */
  he_gemv(&HinvThetaTCCTQQCCAAx, HinvThetaTCCTQQCCAA, &xdiff, rk);
  /* Hinv*Theta.T*CC.T*QQ*CC*BB*(uhat-ur) */
  he_gemv(&HinvThetaTCCTQQCCBBu, HinvThetaTCCTQQCCBB, &udiff, rk);
  /* du = -Hinv*Theta.T*CC.T*QQ*CC*(AA*(xhat-xr)+BB*(uhat-ur)) */
  he_add(&du, &HinvThetaTCCTQQCCAAx, &HinvThetaTCCTQQCCBBu);
  he_neg(&du);
  /* u = uhat+du */
  he_copy_ct(&uhat_copy, uhat);
  he_moddown(&uhat_copy);
  he_add(up, &uhat_copy, &du);
  /* release */
  he_free_ct(&xdiff);
  he_free_ct(&udiff);
  he_free_ct(&du);
  he_free_ct(&uhat_copy);
  he_free_ct(&HinvThetaTCCTQQCCAAx);
  he_free_ct(&HinvThetaTCCTQQCCBBu);
}

END_DECLS
