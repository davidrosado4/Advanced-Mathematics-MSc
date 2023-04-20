/*
this is a general purpose package to integrate ordinary differential
equations. the method used is based on a Runge-Kutta-Fehlberg algorithm
of orders 4 and 5, with automatic stepsize control.
Version 1.3.0, Angel Jorba <angel@maia.ub.es>, February 13 2019
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define E1 (1.e0/360.e0)
#define E2 (-128.e0/4275.e0)
#define E3 (-2197.e0/75240.e0)
#define E4 (1.e0/50.e0)
#define E5 (2.e0/55.e0)

#define RK4_1 (25.e0/216.e0)
#define RK4_2 (1408.e0/2565.e0)
#define RK4_3 (2197.e0/4104.e0)
#define RK4_4 (-1.e0/5.e0)

#define RK5_1 (16.e0/135.e0)
#define RK5_2 (6656.e0/12825.e0)
#define RK5_3 (28561.e0/56430.e0)
#define RK5_4 (-9.e0/50.e0)
#define RK5_5 (2.e0/55.e0)

#define C1 (12.e0/13.e0)
#define C2 (1932.e0/2197.e0)
#define C3 (-7200.e0/2197.e0)
#define C4 (7296.e0/2197.e0)
#define C5 (439.e0/216.e0)
#define C6 (3680.e0/513.e0)
#define C7 (-845.e0/4104.e0)
#define C8 (-8.e0/27.e0)
#define C9 (-3544.e0/2565.e0)
#define C10 (1859.e0/4104.e0)

static inline double error_rk(int n, double k[6][n]);
static inline double step_rk(double *x, int n, double k[6][n]);
static inline void compute_ks(double t, double *x, int n, double h,
       double k[6][n], void (*ode)(double,double*,int,double*));

#define EW 1    /* if 1, prints a message when the error is too large */
#define IS 3    /* number of iterates for the stepsize prediction */
#define FC 0.90 /* safety factor for the stepsize prediction */
#define RK 5    /* Runge-Kutta used to compute the new point */

int rkf45(double *at, double *x, int n, double *ah, int sc, double tol,
          double *atf, double *aer, void (*ode)(double,double*,int,double*))
/*
this routine performs one step of Runge-Kutta-Fehlberg 4(5).
the initial condition (at,x) is changed by a new one corresponding
to the same orbit. the error is controlled by the threshold tol.

parameters:
at:  time. input: time corresponding to the actual initial condition.
           output: new value corresponding to the new initial condition.
x:   position. input: actual initial condition.
               output: new position at time *at.
n:   dimension of the system of odes.
ah:  time step (it can be modified by the routine according to the
     given threshold). if negative, the integration goes backwards.
     on exit, it will contain the time step for the next rkf78 call
sc:  stepsize control.
     0: no stepsize control, the step *ah is used
     !=0: stepsize control according to the threshold tol
tol: threshold to control the integration error
atf: final integration time. if NULL, it is ignored. Otherwise, if the
     stepsize is too large (*at+*ah>*atf), it is reduced so that the new
     time at is exactly atf (in that case, the function returns 1)
aer: if NULL, the routine stops if the estimated error is larger than tol.
     if not NULL, the integration returns the estimated absolute error of
     the performed step. this allows, for instance, to integrate with
     a constant stepsize (sc == 0) and to know an estimate of the error.
ode: pointer to the the vectorfield. The parameters of the function are
     are: (t,x,n,f), where t is the time, x the position vector, n
     the dimension and f the value of the vectorfield at (t,x)

returned value:
     0: ok.
     1: ok, and at=tf.
*/
{
  double k[6][n],x0[n],ea,ee,t,h,hn,nr,s;
  int i,j,flag=0;

  t=*at;
  h=*ah;
  if (atf != NULL)
    {
       if (h > 0) {if (t+h > *atf) {h=*atf-t; flag=1;}}
       else {if (t+h < *atf) {h=*atf-t; flag=1;}}
    }
  compute_ks(t,x,n,h,k,ode);
  if (sc == 0)
    {
      if (aer != NULL) *aer=error_rk(n,k);
      step_rk(x,n,k);
      *at=t+h;
      return flag;
    }
  for (j=0; j<n; j++) x0[j]=x[j];
  nr=step_rk(x,n,k);
  ee=(1+nr)*tol/2;
  ea=error_rk(n,k)+1.e-16*nr;
  s=FC*pow(ee/ea,0.2);
  if (2.0 < s) s=2.0;
  hn=h*s;
  if (ea < ee)
    {
      if (aer != NULL) *aer=ea;
      *at=t+h;
      *ah=hn;
      return flag;
    }
  for (i=0; i<IS; i++)
  {
    for (j=0; j<n; j++) x[j]=x0[j];
    h=hn;
    compute_ks(t,x,n,h,k,ode);
    nr=step_rk(x,n,k);
    ee=(1+nr)*tol/2;
    ea=error_rk(n,k)+1.e-16*nr;
    s=FC*pow(ee/ea,0.2);
    if (1.0 < s) s=1.0;
    hn=h*s;
    if (ea < ee)
      {
        if (aer != NULL) *aer=ea;
        *at=t+h;
        *ah=hn;
        return 0;
      }
  }
  if (aer == NULL)
    {
      puts("rkf45 error 1.");
      puts("this message appears because aer is NULL.");
      puts("It means that the stepsize cannot be adjusted to");
      puts("match the required accuracy. To (try to) solve it:");
      puts("* check that the initial stepsize is reasonable (and != 0)");
      puts("* reduce your accuracy requirement");
      printf("time: %g\n",t);
      printf("rk estimated error: %e  actual threshold: %e\n",ea,ee);
      printf("last stepsize: %g\n",h);
      exit(1);
    }
#if EW == 1
  printf("rkf45: t=%e  est. abserr=%e  est. relerr=%e\n",t,ea,ea/nr);
#endif
  *at=t+h;
  *ah=hn;
  *aer=ea;
  return 0;
}
static inline double step_rk(double *x, int n, double k[6][n])
{
  double nr;
  int i;
  nr=0;
#if RK == 4
  for (i=0; i<n; i++)
  {
    x[i] += RK4_1*k[0][i]+RK4_2*k[2][i]+RK4_3*k[3][i]+RK4_4*k[4][i];
    nr += x[i]*x[i];
  }
#elif RK == 5
  for (i=0; i<n; i++)
  {
    x[i] += RK5_1*k[0][i]+RK5_2*k[2][i]+RK5_3*k[3][i]+RK5_4*k[4][i]+RK5_5*k[5][i];
    nr += x[i]*x[i];
  }
#else
  puts("rkf45 error. the value of RK '#defined' in the rkf45.c file");
  puts("has to be either 4 or 5.");
#endif
  nr /= n;
  return sqrt(nr);
}
static inline double error_rk(int n, double k[6][n])
/*
error of the Runge-Kutta 4
*/
{
  double e,s;
  int i;
  e=0;
  for (i=0; i<n; i++)
  {
    s=E1*k[0][i]+E2*k[2][i]+E3*k[3][i]+E4*k[4][i]+E5*k[5][i];
    e += s*s;
  }
  e /= n;
  return sqrt(e);
}
static inline void compute_ks(double t, double *x, int n, double h,
              double k[6][n], void (*ode)(double,double*,int,double*))
{
  double aux[n];
  int i;

  (*ode)(t,x,n,k[0]);
  for (i=0; i<n; i++)
  {
    k[0][i] *= h;
    aux[i]=x[i]+0.25e0*k[0][i];
  }
  (*ode)(t+0.25e0*h,aux,n,k[1]);
  for (i=0; i<n; i++)
  {
    k[1][i] *= h;
    aux[i]=x[i]+0.09375e0*k[0][i]+0.28125e0*k[1][i];
  }
  (*ode)(t+0.375e0*h,aux,n,k[2]);
  for (i=0; i<n; i++)
  {
    k[2][i] *= h;
    aux[i]=x[i]+C2*k[0][i]+C3*k[1][i]+C4*k[2][i];
  }
  (*ode)(t+C1*h,aux,n,k[3]);
  for (i=0; i<n; i++)
  {
    k[3][i] *= h;
    aux[i]=x[i]+C5*k[0][i]-8.e0*k[1][i]+C6*k[2][i]+C7*k[3][i];
  }
  (*ode)(t+h,aux,n,k[4]);
  for (i=0; i<n; i++)
  {
    k[4][i] *= h;
    aux[i]=x[i]+C8*k[0][i]+2*k[1][i]+C9*k[2][i]+C10*k[3][i]-0.275e0*k[4][i];
  }
  (*ode)(t+0.5e0*h,aux,n,k[5]);
  for (i=0; i<n; i++) k[5][i] *= h;

  return;
 }
