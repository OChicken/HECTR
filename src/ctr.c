/*
 * Controller.
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
#include "../pmu/pmu.h"
#include "../gpqhe/gpqhe.h"
#include <stdio.h>

BEGIN_DECLS

void ctr_c2d(double expA[], double expB[],
  const unsigned int n, const double jacA[], const double dt)
{
  const unsigned int nn = n*n;
  const unsigned int n2 = n*2;
  double C[n2*n2], expC[n2*n2];
  memset(C, 0, sizeof(C));
  for (unsigned int i=0; i<nn; i++)
    C[(i/n)*n2+i%n] = jacA[i]*dt;
  for (unsigned int i=0; i<n; i++)
    C[i*n2 + n+i] = dt;
  unsigned int shift = n*n2;
  double eps1 = eps(1);
  for (unsigned int i=shift, all=n2*n2; i<all; i++)
    C[i] = eps1;
  dexpm(expC, n2, C);
  for (unsigned int i=0; i<nn; i++) {
    expA[i] = expC[(i/n)*n2 +   i%n];
    expB[i] = expC[(i/n)*n2 + n+i%n];
  }
}

void ctr_weighting_matrices(double Q[], double R[],
  const unsigned int nx, const unsigned int nu,
  const double xs[], const double us[])
{
  memset(Q, 0, nx*nx*sizeof(double));
  memset(R, 0, nu*nu*sizeof(double));
  for (unsigned int i=0; i<nx; i++)
    Q[i*nx+i] = 1/(xs[i]*xs[i]);
  for (unsigned int i=0; i<nu; i++)
    R[i*nu+i] = 1/(us[i]*us[i]);
}

void ctr_estimator(double Lx[], double Ld[],
  const unsigned int nx, const unsigned int nu, const unsigned int nd, const unsigned int ny,
  const double A[], const double B[], const double C[],
  const double Bd[], const double Cd[],
  const double xs[])
{
  assert(((nd!=0)&&Ld&&Bd&&Cd) || ((nd==0)&&(!Ld)&&(!Bd)&&(!Cd)));
  _Bool disturbance;
  if ((nd!=0)&&Ld&&Bd&&Cd)
    disturbance=1;
  else
    disturbance=0;
  const unsigned int na = nx+nd; /* 'a' stands for 'augmented' */

  /* Aaug[(nx+nd)*(nx+nd)] = [A[nx*nx], Bd[nx*nd]; zeros[nd*nx], eye[nd*nd]] */
  double Aaug[na*na];
  memset(Aaug, 0, sizeof(Aaug));
  for (unsigned int i=0, all=nx*nx; i<all; i++)
    Aaug[(i/nx)*na+i%nx] = A[i];
  if (disturbance)
    for (unsigned int i=0, all=nx*nd; i<all; i++)
      Aaug[(i/nd)*na+nx+i%nd] = Bd[i];
  for (unsigned int i=nx; i<na; i++)
    Aaug[i*na+i]=1;

  /* Baug[(nx+nd)*nu] = [B[nx*nu]; zeros[nd*nu]] */
  double Baug[na*nu];
  memset(Baug, 0, sizeof(Baug));
  memcpy(Baug, B, nx*nu*sizeof(double));

  /* Caug[ny*(ny+nd)] = [C[ny*ny], Cd[ny*nd]] */
  double Caug[ny*na];
  memset(Caug, 0, sizeof(Caug));
  for (unsigned int i=0, all=ny*nx; i<all; i++)
    Caug[(i/nx)*na+i%nx] = C[i];
  if (disturbance)
    for (unsigned int i=0, all=ny*nd; i<all; i++)
      Caug[(i/nd)*na+ny+i%nd] = Cd[i];

  /* Qw[(nx+nd)*(nx+nd)] */
  double Qw[na*na];
  memset(Qw, 0, sizeof(Qw));
  for (unsigned int i=0; i<na; i++)
    Qw[i*na+i] = HECTR_SMALL;
  Qw[na*na-1] = 1;

  /* Rv[ny*ny] */
  double Rv[ny*ny];
  memset(Rv, 0, sizeof(Rv));
  for (unsigned int i=0; i<ny; i++)
    Rv[i*ny+i] = HECTR_SMALL*xs[i]*xs[i];

  double L[na*ny];
  dlqe(L, ny, na, Aaug, Caug, Qw, Rv);
  memcpy(Lx, L      , nx*ny*sizeof(double));
  if (disturbance)
    memcpy(Ld, L+nx*ny, nd*ny*sizeof(double));
}

void ctr_selector(double Ginv[],
  const unsigned int nx, const unsigned int nu, const unsigned int ny,
  const double A[], const double B[], const double C[], const double H[])
{
  const unsigned int ldg = nx+nu;
  /* I-A */
  double IsA[nx*nx];
  for (unsigned int i=0, all=nx*nx; i<all; i++) {
    IsA[i] = -A[i];
    if ((i/nx)==(i%nx))
      IsA[i]+=1;
  }

  /* HC[nu*nx] = H[nu*ny] * C[ny*nx] */
  double HC[nu*nx];
  memset(HC, 0, sizeof(HC));
  for (unsigned int i=0, all=nu*ny; i<all; i++)
    for (unsigned int j=0; j<ny; j++)
      HC[i] += H[(i/nx)*ny+j]*C[j*nx+i%nx];

  /* G[(nx+nu)*(nx+nu)] = [I-A, -B; H*C, zeros(nu*nu)] */
  double G[(nx+nu)*(nx+nu)];
  memset(G, 0, sizeof(G));
  for (unsigned int i=0, all=nx*nx; i<all; i++)
    G[(i/nx)*ldg+i%nx] = IsA[i];
  for (unsigned int i=0, all=nx*nu; i<all; i++)
    G[(i/nu)*ldg+nx+i%nu] = -B[i];
  for (unsigned int i=0, all=nu*nx; i<all; i++)
    G[(nx+i/nx)*ldg+i%nx] = HC[i];

  /* Ginv = inv(G) */
  memset(Ginv, 0, (nx+nu)*(nx+nu)*sizeof(double));
  dgeinv(Ginv, nx+nu, G);
}

void ctr_measure(double y[],
  const unsigned int ny, const unsigned int nx,
  const double C[], const double x[])
{
  memset(y, 0, ny*sizeof(double));
  for (unsigned int i=0; i<ny; i++)
    for (unsigned int j=0; j<nx; j++)
      y[i] += C[i*nx+j]*x[i];
}

void ctr_measure_forward(double xhat[], double dhat[],
  const unsigned int nx, const unsigned int nd, const unsigned int ny,
  const double C[],  const double Cd[],
  const double Lx[], const double Ld[],
  const double y[], const double xhatm[], const double dhatm[])
{
  assert(((nd!=0)&&dhat&&Ld&&dhatm) || ((nd==0)&&(!dhat)&&(!Ld)&&(!dhatm)));

  double Cxhatm[ny];
  memset(Cxhatm, 0, sizeof(Cxhatm));
  double Cddhatm[ny];
  memset(Cddhatm, 0, sizeof(Cddhatm));
  double e[ny];
  memset(e, 0, sizeof(e));

  if ((nd!=0)&&dhat&&Ld&&dhatm)
    goto disturbance;
  else
    goto nodisturbance;

disturbance:
  /* (C*xhatm)[ny] = C[ny*nx] * xhatm[nx] */
  for (unsigned int i=0; i<ny; i++)
    for (unsigned int j=0; j<nx; j++)
      Cxhatm[i] += C[i*nx+j]*xhatm[j];
  /* (Cd*xhatm)[ny] = Cd[ny*nd] * dhatm[nd] */
  for (unsigned int i=0; i<ny; i++)
    for (unsigned int j=0; j<nd; j++)
      Cddhatm[i] += Cd[i*nd+j]*dhatm[j];
  /* e[ny] = y - C*xhatm - Cd*xhatm */
  for (unsigned int i=0; i<ny; i++)
    e[i] = y[i] - Cxhatm[i] - Cddhatm[i];
  /* xhat[nx] = xhatm[nx] + Lx[nx*ny]*e[ny] */
  memset(xhat, 0, nx*sizeof(double));
  for (unsigned int i=0; i<nx; i++) {
    for (unsigned int j=0; j<ny; j++)
      xhat[i] += Lx[i*ny+j]*e[j];
    xhat[i] += xhatm[i];
  }
  /* dhat[nd] = dhatm[nd] + Ld[nd*ny]*e[ny] */
  memset(dhat, 0, nd*sizeof(double));
  for (unsigned int i=0; i<nd; i++) {
    for (unsigned int j=0; j<ny; j++)
      dhat[i] += Ld[i*ny+j]*e[j];
    dhat[i] += dhatm[i];
  }
  return;

nodisturbance:
  /* (C*xhatm)[ny] = C[ny*nx] * xhatm[nx] */
  for (unsigned int i=0; i<ny; i++)
    for (unsigned int j=0; j<nx; j++)
      Cxhatm[i] += C[i*nx+j]*xhatm[j];
  /* e[ny] = y - C*xhatm - Cd*xhatm */
  for (unsigned int i=0; i<ny; i++)
    e[i] = y[i] - Cxhatm[i];
  /* xhat[nx] = xhatm[nx] + Lx[nx*ny]*e[ny] */
  memset(xhat, 0, nx*sizeof(double));
  for (unsigned int i=0; i<nx; i++) {
    for (unsigned int j=0; j<ny; j++)
      xhat[i] += Lx[i*ny+j]*e[j];
    xhat[i] += xhatm[i];
  }
}

void ctr_select(double xr[], double ur[],
  const unsigned int nx, const unsigned int nu, const unsigned int nd, const unsigned int ny,
  const double Bd[], const double Cd[],
  const double H[], const double Ginv[],
  const double dhat[], const double rsp[])
{
  assert(((nd!=0)&&Bd&&Cd&&dhat) || ((nd==0)&&(!Bd)&&(!Cd)&&(!dhat)));

  double Bddhat[nx];
  memset(Bddhat, 0, sizeof(Bddhat));
  double Cddhat[ny];
  memset(Cddhat, 0, sizeof(Cddhat));
  double rsp_sub_HCddhat[nu];
  memset(rsp_sub_HCddhat, 0, sizeof(rsp_sub_HCddhat));
  double pack[nx+nu];

  if ((nd==0)&&(!Bd)&&(!Cd)&&(!dhat))
    goto packing;

  /* (Bd*dhat)[nx] = Bd[nx*nd] * dhat[nd] */
  for (unsigned int i=0; i<nx; i++)
    for (unsigned int j=0; j<nd; j++)
      Bddhat[i] += Bd[i*nd+j]*dhat[j];

  /* Cddhat[ny] = Cd[ny*nd]*dhat[nd] */
  for (unsigned int i=0; i<ny; i++)
    for (unsigned int j=0; j<nd; j++)
      Cddhat[i] += Cd[i*nd+j]*dhat[j];

  /* rsp[nu] - H[nu*ny]*(Cd*dhat)[ny] */
  for (unsigned int i=0; i<nu; i++) {
    for (unsigned int j=0; j<ny; j++)
      rsp_sub_HCddhat[i] += -H[i*ny+j]*Cddhat[j];
    rsp_sub_HCddhat[i] += rsp[i];
  }

packing:
  /* pack[nx+nu] = [(Bd*dhat)[nx]; (rsp-H*Cd*dhat))[nu]] */
  memcpy(&pack[0] , Bddhat,          nx*sizeof(double));
  memcpy(&pack[nx], rsp_sub_HCddhat, nu*sizeof(double));

  /* r[nx+nu] = Ginv[nx+nu]*r[nx+nu] */
  double r[nx+nu];
  memset(r, 0, sizeof(r));
  for (unsigned int i=0, ldr=nx+nu; i<ldr; i++)
    for (unsigned int j=0; j<ldr; j++)
      r[i] += Ginv[i*ldr+j]*pack[j];
  memcpy(xr, &r[0],  nx*sizeof(double));
  memcpy(ur, &r[nx], nu*sizeof(double));
}

void ctr_control(double u[],
  const unsigned int nx, const unsigned int nu,
  const double G[], const double xhat[], const double xr[], const double ur[])
{
  memset(u, 0, nu*sizeof(double));
  for (unsigned int i=0; i<nu; i++) {
    for (unsigned int j=0; j<nx; j++)
      u[i] += -G[i*nx+j]*(xhat[j]-xr[j]);
    u[i] += ur[i];
  }
}

void ctr_estimate(double xhatm[], double dhatm[],
  const unsigned int nx, const unsigned int nu, const unsigned int nd,
  const double A[], const double B[], const double Bd[],
  const double xhat[], const double dhat[], const double u[])
{
  assert(((nd!=0)&&dhatm&&dhat) || ((nd==0)&&(!dhatm)&&(!dhat)));
  double Axhat[nx], Bu[nx], Bddhat[nx];
  memset(Axhat, 0, sizeof(Axhat));
  memset(Bu, 0, sizeof(Bu));
  memset(Bddhat, 0, sizeof(Bddhat));

  if ((nd!=0)&&dhatm&&dhat)
    goto disturbance;
  else
    goto nodisturbance;

disturbance:
  for (unsigned int i=0; i<nx; i++) {
    for (unsigned int j=0; j<nx; j++)
      Axhat[i] += A[i*nx+j]*xhat[j];
    for (unsigned int j=0; j<nu; j++)
      Bu[i] += B[i*nu+j]*u[j];
    for (unsigned int j=0; j<nd; j++)
      Bddhat[i] += Bd[i*nd+j]*dhat[j];
    xhatm[i] = Axhat[i]+Bu[i]+Bddhat[i];
  }
  for (unsigned int i=0; i<nd; i++)
    dhatm[i] = dhat[i];
  return;

nodisturbance:
  for (unsigned int i=0; i<nx; i++) {
    for (unsigned int j=0; j<nx; j++)
      Axhat[i] += A[i*nx+j]*xhat[j];
    for (unsigned int j=0; j<nu; j++)
      Bu[i] += B[i*nu+j]*u[j];
    xhatm[i] = Axhat[i]+Bu[i];
  }
}

void ctr_actuate(double xnext[],
  const unsigned int nx, const unsigned int nu, const unsigned int np,
  const double x[], const double u[], const double p[],
  const double xs[], const double us[], const double ps[],
  void (*const ode     )(double[], const double[], const double[], const double[]),
  void (*const jacobian)(double[], const double[], const double[], const double[]),
  const double dt)
{
  double xx[nx], uu[nu], pp[np];
  for (unsigned int i=0; i<nx; i++)
    xx[i] = x[i]+xs[i];
  for (unsigned int i=0; i<nu; i++)
    uu[i] = u[i]+us[i];
  for (unsigned int i=0; i<np; i++)
    pp[i] = p[i]+ps[i];
  double ddt = dt/2;
  for (unsigned int n=0, N=(int)dt/ddt; n<N; n++)
    ode15s(xx, nx, xx, uu, pp, ode, jacobian, ddt);
  for (unsigned int i=0; i<nx; i++)
    xnext[i] = xx[i]-xs[i];
}

/**
 * @nx: dimension of state
 * @nu: dimension of control
 * @np: dimension of controlled parameters
 * @ny: dimension of observables
 * @nd: dimension of disturbance variables
 */
void ctr_simulate(double x[], double u[], double p[],
  const unsigned int nx, const unsigned int nu, const unsigned int np, /* state, control, parameters */
  const unsigned int ny, const unsigned int nd, /* observe, disturbance */
  const double A[], const double B[], const double C[],
  const double Bd[], const double Cd[], const double Hr[],
  const double xs[], const double us[], const double ps[],
  void (*const ode     )(double[], const double[], const double[], const double[]),
  void (*const jacobian)(double[], const double[], const double[], const double[]),
  const double dt, const unsigned int N)
{
  const unsigned int horizon = N/10;

  /* weighting matrices */
  double Q[nx*nx], R[nu*nu];
  ctr_weighting_matrices(Q, R, nx, nu, xs, us);

  /* regulator */
  double G[nu*nx], X[nx*nx];
  dlqr(G, X, nx, nu, A, B, Q, R);

  /* estimator */
  double Lx[nx*ny], Ld[nd*ny];
  ctr_estimator(Lx, Ld, nx, nu, nd, ny, A, B, C, Bd, Cd, xs);

  /* selector */
  double Ginv[(nx+nu)*(nx+nu)];
  ctr_selector(Ginv, nx, nu, ny, A, B, C, Hr);

  /* measure and estimate */
  double y[ny*(N+1)];
  double xhatm[nx*(N+1)], dhatm[nd*(N+1)];
  double xhat[nx*N], dhat[nd*N];
  for (unsigned int i=0; i<nx; i++)
    x[i]=0, xhatm[i]=0;
  for (unsigned int i=0; i<nd; i++)
    dhatm[i]=0;

  /* set point */
  double xr[nx*N], ur[nu*N], rsp[nu*N];
  memset(xr,  0, sizeof(xr));
  memset(ur,  0, sizeof(ur));
  memset(rsp, 0, sizeof(rsp));

  TEST_DO("closed-loop simulate");
  for (unsigned int k=0; k<N+1; k++) {
    /* Take measurement */
    ctr_measure(&y[k*ny], ny, nx, C, &x[k*nx]);
    /* Advance state measurement */
    ctr_measure_forward(&xhat[k*nx], &dhat[k*nd],
      nx, nd, ny, C, Cd, Lx, Ld, &y[k*ny], &xhatm[k*nx], &dhatm[k*nd]);
    /* Stop if at last time */
    printf("\r%2u/%u completed.", k, N);
    fflush(stdout);
    if (k==N) break;
    /* Use steady-state target selector */
    ctr_select(&xr[k*nx], &ur[k*nu], nx, nu, nd, ny, Bd, Cd, Hr, Ginv, &dhat[k*nd], &rsp[k*nu]);
    if (k==0)
      for (unsigned int i=0; i<nu; i++)
        u[i]=ur[i];
    /* Apply control law */
    //ctr_control(&u[k*nu], nx, nu, G, &xhat[k*nx], &xr[k*nx], &ur[k*nu]);
    double up[nu*horizon];/* u predictive */
    ctr_mpc(up, ny, nx, nu, horizon, A, B, C, Q, R,
      &xhat[k*nx], (k==0)? &u[0] : &u[(k-1)*nu], &xr[k*nx], &ur[k*nu], 
      NULL, NULL, NULL, NULL, NULL, NULL);
    memcpy(&u[k*nu], up, nu*sizeof(double));
    /* Evolve plant. Our variables are deviation but cstrsim needs positional. */
    ctr_actuate(&x[(k+1)*nx], nx, nu, np, &x[k*nx], &u[k*nu], &p[k*np], xs, us, ps, ode, jacobian, dt);
    /* advance state estimates */
    ctr_estimate(&xhatm[(k+1)*nx], &dhatm[(k+1)*nd], nx, nu, nd, A, B, Bd, &xhat[k*nx], &dhat[k*nd], &u[k*nu]);
  }
  printf("\n");
  TEST_DONE();
  /* Convert to positional units */
  for (unsigned int k=0; k<N+1; k++)
    for (unsigned int i=0; i<nx; i++)
      x[k*nx+i]+=xs[i];
  for (unsigned int k=0; k<N; k++)
    for (unsigned int i=0; i<nu; i++)
      u[k*nu+i]+=us[i];
}

static void hectr_enc_states(he_ct_t *ct_up, he_ct_t *ct_xhat, he_ct_t *ct_uhat, he_ct_t *ct_xr, he_ct_t *ct_ur,
  const unsigned int n, const unsigned int m, const unsigned int slots,
  const double xhat[], const double uhat[], const double xr[], const double ur[],
  const he_pk_t *pk)
{
  _Complex double zeros[slots];
  _Complex double xhatz[slots];
  _Complex double uhatz[slots];
  _Complex double xrz[slots];
  _Complex double urz[slots];
  memset(zeros, 0, sizeof(zeros));
  d2z_vector(xhatz, n, slots, xhat);
  d2z_vector(uhatz, m, slots, uhat);
  d2z_vector(xrz,   n, slots, xr);
  d2z_vector(urz,   m, slots, ur);
  he_pt_t pt_up, pt_xhat, pt_uhat, pt_xr, pt_ur;
  he_alloc_pt(&pt_up);
  he_alloc_pt(&pt_xhat);
  he_alloc_pt(&pt_uhat);
  he_alloc_pt(&pt_xr);
  he_alloc_pt(&pt_ur);
  he_ecd(&pt_up,   zeros);
  he_ecd(&pt_xhat, xhatz);
  he_ecd(&pt_uhat, uhatz);
  he_ecd(&pt_xr,   xrz);
  he_ecd(&pt_ur,   urz);
  he_enc_pk(ct_up,   &pt_up, pk);
  he_enc_pk(ct_xhat, &pt_xhat, pk);
  he_enc_pk(ct_uhat, &pt_uhat, pk);
  he_enc_pk(ct_xr,   &pt_xr, pk);
  he_enc_pk(ct_ur,   &pt_ur, pk);
  he_free_pt(&pt_up);
  he_free_pt(&pt_xhat);
  he_free_pt(&pt_uhat);
  he_free_pt(&pt_xr);
  he_free_pt(&pt_ur);
}

static void hectr_dec_state(double u[],
  const unsigned int m, const unsigned int slots,
  const he_ct_t *ct_up, const poly_mpi_t *sk)
{
  he_pt_t pt_up;
  he_alloc_pt(&pt_up);
  he_dec(&pt_up, ct_up, sk);
  _Complex double uz[slots];
  memset(uz, 0, sizeof(uz));
  he_dcd(uz, &pt_up);
  for (unsigned int i=0; i<slots; i++)
    assert(cimag(uz[i])<HECTR_SMALL);
  for (unsigned int i=0; i<m; i++)
    u[i] = creal(uz[i]);
  he_free_pt(&pt_up);
}

void hectr_simulate(double x[], double u[], double p[],
  const unsigned int nx, const unsigned int nu, const unsigned int np, /* state, control, parameters */
  const unsigned int ny, const unsigned int nd, /* observe, disturbance */
  const double A[], const double B[], const double C[],
  const double Bd[], const double Cd[], const double Hr[],
  const double xs[], const double us[], const double ps[],
  void (*const ode     )(double[], const double[], const double[], const double[]),
  void (*const jacobian)(double[], const double[], const double[], const double[]),
  const double dt, const unsigned int N)
{
  const unsigned int horizon = N/10;
  const unsigned int slots = 1<<(32-__builtin_clz(nu*horizon));

  /* encryption settings */
  unsigned int logn = 12;
  MPI q = mpi_set_ui(NULL, 1);
  mpi_lshift(q, q, 109);
  uint64_t Delta = 1UL<<50;
  hectx_init(logn, q, slots, Delta);
  poly_mpi_t sk;
  he_pk_t pk;
  he_evk_t rk[slots];

  /* init keys */
  he_alloc_pk(&pk);
  he_alloc_sk(&sk);
  for (unsigned int rot=0; rot<slots; rot++)
    he_alloc_evk(&rk[rot]);
  TEST_DO("he_keypair");
  he_keypair(&pk, &sk);
  TEST_DONE();
  TEST_DO("he_genrk");
  he_genrk(rk, &sk);
  TEST_DONE();

  /* ct */
  he_ct_t ct_xhat, ct_uhat, ct_xr, ct_ur, ct_up;
  he_alloc_ct(&ct_xhat);
  he_alloc_ct(&ct_uhat);
  he_alloc_ct(&ct_xr);
  he_alloc_ct(&ct_ur);
  he_alloc_ct(&ct_up); /* u predictive */

  /* weighting matrices */
  double Q[nx*nx], R[nu*nu];
  ctr_weighting_matrices(Q, R, nx, nu, xs, us);

  /* estimator */
  double Lx[nx*ny], Ld[nd*ny];
  ctr_estimator(Lx, Ld, nx, nu, nd, ny, A, B, C, Bd, Cd, xs);

  /* selector */
  double Ginv[(nx+nu)*(nx+nu)];
  ctr_selector(Ginv, nx, nu, ny, A, B, C, Hr);

  /* measure and estimate */
  double y[ny*(N+1)];
  double xhatm[nx*(N+1)], dhatm[nd*(N+1)];
  double xhat[nx*N], dhat[nd*N];
  for (unsigned int i=0; i<nx; i++)
    x[i]=0, xhatm[i]=0;
  for (unsigned int i=0; i<nd; i++)
    dhatm[i]=0;

  /* set point */
  double xr[nx*N], ur[nu*N], rsp[nu*N];
  memset(xr,  0, sizeof(xr));
  memset(ur,  0, sizeof(ur));
  memset(rsp, 0, sizeof(rsp));

  TEST_DO("closed-loop simulate");
  for (unsigned int k=0; k<N+1; k++) {
    /* Take measurement */
    ctr_measure(&y[k*ny], ny, nx, C, &x[k*nx]);
    /* Advance state measurement */
    ctr_measure_forward(&xhat[k*nx], &dhat[k*nd],
      nx, nd, ny, C, Cd, Lx, Ld, &y[k*ny], &xhatm[k*nx], &dhatm[k*nd]);
    /* Stop if at last time */
    printf("\r%2u/%u completed.", k, N);
    fflush(stdout);
    if (k==N) break;
    /* Use steady-state target selector */
    ctr_select(&xr[k*nx], &ur[k*nu], nx, nu, nd, ny, Bd, Cd, Hr, Ginv, &dhat[k*nd], &rsp[k*nu]);
    if (k==0)
      for (unsigned int i=0; i<nu; i++)
        u[i]=ur[i];
    /* Apply control law */
    hectr_enc_states(&ct_up, &ct_xhat, &ct_uhat, &ct_xr, &ct_ur, nx, nu, slots,
      &xhat[k*nx], (k==0)? &u[0] : &u[(k-1)*nu], &xr[k*nx], &ur[k*nu], &pk);
    ctr_hempc(&ct_up, ny, nx, nu, horizon, slots, A, B, C, Q, R, &ct_xhat, &ct_uhat, &ct_xr, &ct_ur, rk);
    hectr_dec_state(&u[k*nu], nu, slots, &ct_up, &sk);
    /* Evolve plant. Our variables are deviation but cstrsim needs positional. */
    ctr_actuate(&x[(k+1)*nx], nx, nu, np, &x[k*nx], &u[k*nu], &p[k*np], xs, us, ps, ode, jacobian, dt);
    /* advance state estimates */
    ctr_estimate(&xhatm[(k+1)*nx], &dhatm[(k+1)*nd], nx, nu, nd, A, B, Bd, &xhat[k*nx], &dhat[k*nd], &u[k*nu]);
  }
  printf("\n");
  TEST_DONE();
  /* Convert to positional units */
  for (unsigned int k=0; k<N+1; k++)
    for (unsigned int i=0; i<nx; i++)
      x[k*nx+i]+=xs[i];
  for (unsigned int k=0; k<N; k++)
    for (unsigned int i=0; i<nu; i++)
      u[k*nu+i]+=us[i];

  /* release */
  mpi_release(q);
  he_free_sk(&sk);
  he_free_pk(&pk);
  for (unsigned int rot=0; rot<slots; rot++)
    he_free_evk(&rk[rot]);
  he_free_ct(&ct_xhat);
  he_free_ct(&ct_uhat);
  he_free_ct(&ct_xr);
  he_free_ct(&ct_ur);
  he_free_ct(&ct_up);
  hectx_exit();
}


#if 0
  printf("A=\n");
  for (unsigned int i=0, all=nx*nx; i<all; i++) {
    printf("%12.8f ", A[i]);
    if ((i+1)%nx==0) printf("\n");
  }
  printf("B=\n");
  for (unsigned int i=0, all=nx*nu; i<all; i++) {
    printf("%12.8f ", B[i]);
    if ((i+1)%nu==0) printf("\n");
  }
  printf("Bp=\n");
  for (unsigned int i=0, all=nx*np; i<all; i++) {
    printf("%12.8f ", Bp[i]);
    if ((i+1)%np==0) printf("\n");
  }
  printf("C=\n");
  for (unsigned int i=0; i<ny*ny; i++) {
    printf("%12.8f ", C[i]);
    if ((i+1)%ny==0)printf("\n");
  }
  printf("G=\n");
  for (unsigned int i=0; i<nu*nx; i++) {
    printf("%13e ", G[i]);
    if ((i+1)%nx==0)printf("\n");
  }
  printf("Lx=\n");
  for (unsigned int i=0, all=nx*ny; i<all; i++) {
    printf("%13e ", Lx[i]);
    if ((i+1)%ny==0) printf("\n");
  }
  printf("Ld=\n");
  for (unsigned int i=0, all=nd*ny; i<all; i++) {
    printf("%13e ", Ld[i]);
    if ((i+1)%ny==0) printf("\n");
  }
#endif

END_DECLS
