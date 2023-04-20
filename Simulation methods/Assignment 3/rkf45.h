int rkf45(double *at, double *x, int n, double *ah, int sc, double tol,
          double *atf, double *aer, void (*ode)(double,double*,int,double*));
