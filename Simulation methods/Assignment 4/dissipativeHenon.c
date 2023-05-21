/* David Rosado Rodríguez */
/* NIUB: 20194344 */

/* Solution of the fourth assignment of Simulation Methods */

/* Add needed incldudes */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <complex.h>

/* Definition of the parameters we are going to use */
#define a 1.4
#define b 0.3

/* Function definition */

void dynamics(double, double, int, const char*);
void H(double, double, double*, int);
double Euclidian_norm(double*, int);
void gaussEliminationLS(int, int, double**, double*, double*);
int Newton(double, double, double*, double);
double *solver_real(double, double);
double complex *solver_complex(double, double);
void find_eigenvalues(double, double);
void DH_stimes(double**, int, double*);
void matrix_mult(double**, double**);
double Lyapunov(double, double, int, double*, double);
void mat_times_vect(double**, double*, double*);



/* Main */
int main(void){
    int iter_max,i, iter_newton,s;
    const char* file_name[3] = {"dynamics_(0,0).txt", "dynamics_(1,1).txt", "dynamics_(0.5,0.5).txt"};
    double init_cond[6] = {0,0,1,1,0.5,0.5}, *res, tol, tr, det,*v, lyapunov_value;
    
    /* Allocate memory */
    res = (double *)calloc(2,sizeof(double));
    v = (double *)calloc(2,sizeof(double));
    if(res == NULL || v == NULL){
        printf("Problems with memory in the main\n");
        exit(1);
    }
    
    /*-------------------------------------------------------------------------------------------------------*/
    /* -------------------------------------------------EXERCISE 1-------------------------------------------*/
    /*-------------------------------------------------------------------------------------------------------*/
    
    /* Compute the dynamic of H*/
    
    /* Define parameters to use the dynamics function */
    iter_max = 20000;
    for(i=0; i<3;i++){
        dynamics(init_cond[2*i], init_cond[2*i + 1], iter_max, file_name[i]);
    }
    
    printf("\n\n\033[1;34m-------------------------------EXERCISE 1----------------------------------\033[0m\n\n");
    printf("Multiple files with name \033[1;1mdynamics_(x,y).txt\033[0m have been generated for study the dynamics of H.\nWe provide \033[1;1mthree different files\033[0m, corresponding to \033[1;1mthree different starting points.\033[0m\n");
    
    /*-------------------------------------------------------------------------------------------------------*/
    /* -------------------------------------------------EXERCISE 2-------------------------------------------*/
    /*-------------------------------------------------------------------------------------------------------*/
    
    /* Compute a fixed point near (0.63,0.19) and its stability */
    tol = 1e-15;
    iter_newton = Newton(0.63,0.19,res, tol);
    printf("\n\n\033[1;34m-------------------------------EXERCISE 2-----------------------------------\033[0m\n\n");
    printf("The fixed point found with a threshold of %.le is \033[1;1m(%f,%f)\033[0m.\nThe method has required  \033[1;1m%d iterations.\033[0m\n\n",tol,res[0], res[1], iter_newton);
    
    /* Compute its stability*/
    tr = -2*a*res[0]; det = -b;
    find_eigenvalues(tr,det);
    
    /*-------------------------------------------------------------------------------------------------------*/
    /* -------------------------------------------------EXERCISE 3-------------------------------------------*/
    /*-------------------------------------------------------------------------------------------------------*/
    printf("\n\n\033[1;34m-------------------------------EXERCISE 3-----------------------------------\033[0m\n\n");
    
    /* Define parameters to use Lyapunov function */
    v[0] = 1; v[1] = 1; v[0]/=Euclidian_norm(v,2); v[0]/=Euclidian_norm(v,2);
    s = 1;
    printf("Computation of Lyapunov exponents using \033[1;1mdifferent initial conditions\033[0m and \033[1;1mdifferent values of s\033[0m\n\n\n");
    for(i = 0; i<5 ; i++){
        lyapunov_value = Lyapunov(0,0,s,v,tol);
        printf("An approximaion of the maximal Lyapunov exponent for (x,y) = (0,0) and s = %d is: %f\n", s, lyapunov_value);
        printf("The second Lyapunov exponent is given by: %f.\n", log(b)-lyapunov_value);
        s++;
        
    }
    printf("\n\n");
    s = 1;
    for(i = 0; i<5 ; i++){
        lyapunov_value = Lyapunov(1,1,s,v,tol);
        printf("An approximaion of the maximal Lyapunov exponent for (x,y) = (1,1) and s = %d is: %f\n", s, lyapunov_value);
        printf("The second Lyapunov exponent is given by: %f.\n", log(b)-lyapunov_value);
        s++;
        
    }
    
    /* Free up memory */
    free(res); free(v);
    
    /* Style considerations and return */
    printf("\n\n");
    printf("\033[0;31mScript made by David Rosado Rodríguez.\033[0m\n\n");
    return 0;
}
/*--------------------------------------------------------------------------------------------------------*/
/* ----------------------------FUNCTIONS FOR SOLVING THE FIRST EXERCISE ----------------------------------*/
/*--------------------------------------------------------------------------------------------------------*/

void dynamics(double x, double y, int iter_max, const char* file_name){
    /*
     Function that creates a txt file to plot the dynamics of H given a initial condition (x,y).
     Parameters:
     x: First component of the starting point to compute the dynamics of H.
     y: Second component of the starting point to compute the dynamics of H.
     iter_max: Maximum iterates to compute for the dynamics of H.
     file_name: Name of the file we want to create.
     */
    
    FILE *fit;
    double *res;
    int i;
    
    /* Allocate memory and open the file */
    res = (double *)calloc(2,sizeof(double));
    fit = fopen(file_name, "w");
    if(fit == NULL || res == NULL){
        printf("Problems with memory in the dynamics function");
        exit(1);
    }
    
    /* Compute iterations of H and plot the results in the file */
    for(i = 0; i<iter_max;i++){
        /* Compute the image */
        H(x,y,res,1);
        
        /* Add points to the file */
        fprintf(fit, "%le %le\n", res[0],res[1]);
        
        /* Update variables */
        x = res[0]; y = res[1];
    }
    
    /* Free up memory */
    free(res); fclose(fit);
    
}

/*--------------------------------------------------------------------------------------------------------*/
/* --------------------------------FUNCTIONS FOR SOLVING THE SECOND EXERCISE----------------------------- */
/*--------------------------------------------------------------------------------------------------------*/

int Newton(double x0, double y0, double *res, double tol){
    /*
     Function that given an initial point (x0,y0), computes a Newton method of the function H(x,y)-(x,y).
     Parameters:
     x0: First coordinate of the point to apply the Newton's method.
     y0: Second coordinate of the point to apply the Newton's method.
     res: Pointer to which we assign the results.
     tol: Threshold to stop the Newton loop.
     Returns: The number of iterations needed to perform Newton's method.
     */
    int i, iter_newton = 0;
    double *vector_b, **DH, *dif_vect, *previous_x, difference;
    
    /* Allocate memory */
    vector_b = (double *)calloc(2,sizeof(double));
    dif_vect = (double *)calloc(2,sizeof(double));
    previous_x = (double *)calloc(2,sizeof(double));
    DH = (double **)calloc(2,sizeof(double *));
    for(i=0 ;i<2;i++){
        DH[i] = (double *)calloc(3,sizeof(double));
    }
    if(vector_b == NULL || dif_vect == NULL || previous_x == NULL || DH == NULL){
        printf("Problems with memory in the Newton function\n");
        exit(1);
    }
    
    /* Start Newton loop */
    /* Prepare some variables*/
    difference = 1;
    previous_x[0] = x0; previous_x[1] = y0;
    while(difference > tol && iter_newton < 50){
        
        /* Prepare the Jacobian matrix*/
        DH[0][0] = -2*a*previous_x[0] - 1; DH[0][1] = 1;
        DH[1][0] = b; DH[1][1] = - 1;
        
        /* Prepare the vector b to solve DHx=b */
        H(previous_x[0],previous_x[1],res,1);
        vector_b[0] = -res[0] + previous_x[0]; vector_b[1] = -res[1] + previous_x[1];
        
        /* Gaussian elimination and store the result in res */
        gaussEliminationLS(2,3,DH,vector_b,res);
        
        /* Final solution */
        res[0] += previous_x[0]; res[1] += previous_x[1];
        
        /* Compute the difference between iterates */
        dif_vect[0] = res[0] - previous_x[0]; dif_vect[1] = res[1] - previous_x[1];
        difference = Euclidian_norm(dif_vect,2);
        
        /* Update variables */
        previous_x[0] = res[0]; previous_x[1] = res[1];
        
        iter_newton++;
    }
    
    /* Free up memory */
    free(previous_x); free(vector_b); free(dif_vect);
    for(i = 0; i<2; i++){
        free(DH[i]);
    }
    free(DH);
    
    return iter_newton;
}

/*--------------------------------------------------------------------------------------------------------*/
/* ---------------------------------FUNCTIONS FOR SOLVING THE LAST EXERCISE ------------------------------*/
/*--------------------------------------------------------------------------------------------------------*/
double Lyapunov(double x0, double y0, int s, double *v, double tol){
    /*
     Function that given an initial condition (x0,y0), computes an approximation of the Lyapunov exponent
     of the attractor.
     Parameters:
     x0: First coordinate of the initial point to compute the approximation.
     y0: Second coordinate of the initial point to compute the approximation.
     s: Positive integer.
     v: Vector of norm 1.
     tol: Threshold for the stop criteira.
     Returns: An approximation of the Lyapunov exponent.
     */
    int i,k;
    double *point,*w,**DH, previous_x, previous_y, *sol, alpha, lyapunov = 0, lambda_k = 0, lambda_previous;
    
    /* Allocate memory */
    point = (double*)calloc(2,sizeof(double));
    sol = (double*)calloc(2,sizeof(double));
    w = (double*)calloc(2,sizeof(double));
    DH = (double **)calloc(2,sizeof(double *));
    for(i=0 ;i<2;i++){
        DH[i] = (double *)calloc(3,sizeof(double));
    }
    if(point == NULL || w == NULL || DH == NULL || sol == NULL){
        printf("Problems with memory in the Lyapunov function\n");
        exit(1);
    }
    
    /* Step 1. Take an initial point in the basin of attraction --> (x0,y0)*/
    
    /* We are going to make Step 2 and Step 3 together */
    
    /* Step 2. 10.000 iterations of H given the initial point (x0,y0). Assume Assume (xn,yn) is in the basin of attraction */
    
    /* Step 3. Compute 10.000 iterates more, to see if they fill the attractor */
    H(x0,y0,point,20000);
    
    /* Step 4. Power method +  computation of the Lyapunov exponent */
    
    /* Prepare variables */
    w[0] = v[0]; w[1] = v[1];
    for(k=1;k<50000;k++){
        
        /* Update stop criteria */
        lambda_previous = lambda_k;
        
        /* Compute alpha_{k} */
        
        /* Point to evualte the matrix */
        if(k!=1){
            previous_x = point[0]; previous_y = point[1];
            H(previous_x,previous_y,point,s);
        }
        
        /* Compute the matrix D_{H^s} evaluated in the previous point */
        DH_stimes(DH,s,point);
        mat_times_vect(DH,w,sol);
        
        /* Value of alpha_{k}*/
        alpha = Euclidian_norm(sol,2);
        
        /* Compute w_{k}*/
        w[0] = sol[0]/alpha; w[1] = sol[1]/alpha;
        
        /* Compute Lyapunov values*/
        lyapunov +=log(alpha);
        
        /* Stop criteria*/
        lambda_k = (1./(k*s))*lyapunov;
        if(k!=1 && fabs(lambda_k-lambda_previous)<tol){
            break;
        }
    }
    
    /* Free up memory*/
    free(sol); free(w); free(point);
    for (i = 0; i < 2; i++) {
        free(DH[i]);
    }
    free(DH);
    return lambda_k;
}



/*--------------------------------------------------------------------------------------------------------*/
/*---------------------------------------- GENERAL FUNCTIONS ---------------------------------------------*/
/*--------------------------------------------------------------------------------------------------------*/
void H(double x, double y, double* image, int m){
    /*
     Function that computes m iterates of the Henon map, H, at (x,y).
     When m = 1, the function returns the image of the H at (x,y).
     Parameters:
     x: First coordinate of the point to apply the iterates of H.
     y: Second coordinate of the point to apply the iterates of H.
     image: Pointer to which we will assign the result of the iterates of H.
     m: Iterates of the function that we desire to perform, m > 0.
     */
    
    /* If m < 1, breaks. Otherwise, perform the iterations */
    if( m < 1){
        printf("Error in the number of iterations requested for H\n");
        exit(1);
    }else{
        int i;
        for(i = 0; i<m ; i++){
            image[0] = 1 + y - a*x*x; image[1] = b*x;
            x = image[0]; y = image[1];
        }
    }
}

double Euclidian_norm(double* vector, int n){
    /*
     Function that computes the Euclidan norm of a vector.
     Parameters:
     vector: Pointer to the vector for which we want to calculate the norm.
     n: length of the vector.
     Returns: The Euclidian norm of the vector.
     */
    int i;
    double square_sum = 0;
    
    for(i = 0; i < n; i++){
        square_sum += vector[i] * vector[i];
        }
    return sqrt(square_sum);
}

void gaussEliminationLS(int m, int n, double **mat, double *bv, double *x){
    /*
     Function that solves a linear system Ax=b using Gauss Elimination with partial pivoting.
     The result is stored in x.
     The function recieves a matrix of dimensions m x (n-1) and add the vector bv to the last
     column of the matrix, to start the Guassian Elimination process.
     Paramters:
     m: Rows dimensions of the matrix
     n: Columns + 1 dimensions of the matrix.
     mat: The matrix of the linear system.
     bv: The vector b to solve Ax=b.
     x: A pointer to which we assign the result.
     */
    int i,j,k;
    /* Append the vector to the matrix */
    for(i=0; i<m;i++){
        mat[i][n-1] = bv[i];
    }
    /* Start Gauss */
    for(i=0; i<m-1; i++) {
        /* Partial Pivoting */
        int max_row = i;
        for(k=i+1; k<m; k++) {
            /*If diagonal element(absolute value) is smaller than any of the terms below it */
            if(fabs(mat[i][i]) < fabs(mat[k][i])) {
                /*Track the row interchanges */
                max_row = k;
            }
        }
        if(max_row != i) {
            /* Swap the rows */
            for(j=0; j<n; j++) {
                double temp = mat[i][j];
                mat[i][j] = mat[max_row][j];
                mat[max_row][j] = temp;
            }
        }
        /* Begin Gauss Elimination */
        for(k=i+1; k<m; k++) {
            double term = mat[k][i] / mat[i][i];
            for(j=0; j<n; j++) {
                mat[k][j] = mat[k][j] - term*mat[i][j];
            }
        }
    }
    /* Begin Back-substitution */
    for(i=m-1; i>=0; i--) {
        x[i] = mat[i][n-1];
        for(j=i+1; j<n-1; j++) {
            x[i] = x[i] - mat[i][j]*x[j];
        }
        x[i] = x[i] / mat[i][i];
    }
}

double *solver_real(double tr, double det){
    /*
     Function that computes the eigenvalues of a matrix given its trace and determinant if the eigenvalues are real.
     Parameters:
     tr: Trace of the matrix.
     det: Determinant of the matrx.
     Returns: A pointer to both eigenvalues
     */
    
    double *vaps;
    
    vaps = (double *)calloc(2,sizeof(double));
    if(vaps == NULL){
        printf("No memory\n");
        exit(1);
    }
    vaps[0] = (tr + sqrt(tr*tr - 4*det))/2;
    vaps[1] = (tr - sqrt(tr*tr - 4*det))/2;
    return vaps;
}

double complex *solver_complex(double tr, double det){
    /*
     Function that computes the eigenvalues of a matrix given its trace and determinant if the eigenvalues are complex.
     Parameters:
     tr: Trace of the matrix.
     det: Determinant of the matrx.
     Returns: A pointer to both eigenvalues
     */
    
    double complex *vaps;
    
    vaps = (double complex *)calloc(2,sizeof(double complex));
    if(vaps == NULL){
        printf("No memory\n");
        exit(1);
    }
    vaps[0] = (tr + csqrt(tr*tr - 4*det))/2;
    vaps[1] = (tr - csqrt(tr*tr - 4*det))/2;
    return vaps;
}

void find_eigenvalues(double tr, double det){
    /*
     Function that, given the trace and determinant of a matrix, compute if the eigenvalues are real or complex and apply previous functions to determine its eigenvalues.
     Parameters:
     tr: Trace of the matrix.
     det: Determinant of the matrx.
     */
    if( (tr*tr-4*det) >= 0 ){
        double *vaps;
        vaps = solver_real(tr,det);
        printf("The eigenvalues are \033[1;1mreal\033[0m. They are the following: (%f,%f).\n",vaps[0],vaps[1]);
        printf("The modulus of the eigenvalues is: (%f, %f).\n",fabs(vaps[0]),fabs(vaps[1]));
        free(vaps);
    }else{
        double complex *vaps;
        vaps = solver_complex(tr,det);
        printf("The eigenvalues are \033[1;1mcomplex\033[0m. They are the following\n");
        printf("%.2f %+.2fi\n", creal(vaps[0]), cimag(vaps[0]));
        printf("%.2f %+.2fi\n\n", creal(vaps[1]), cimag(vaps[1]));
        printf("The modulus of the eigenvalues is (%f, %f)\n", cabs(vaps[0]),cabs(vaps[1]));
        free(vaps);
    }
}

void DH_stimes(double **A, int s, double *eval_point){
    /*
     Function that computes the matrix D_{H^s} evualted at eval_point. The result is stored in the pointer A.
     Parameters:
     A: The matrix to which we want to compute its derivative.
     s: A positive integer.
     eval_point: The evualtion point.
     */
    if( s < 1 ){
        printf("Error in the positive integer. Check DH_stimes function\n");
        exit(1);
    }
    if(s == 1){
        A[0][0] = -2*a*eval_point[0]; A[0][1] = 1; A[1][0] = b; A[1][1] = 0;
    }else{
        int i,j,n;
        double **aux, previous_x = eval_point[0], previous_y = eval_point[1], *image;
        
        /* Allocate memory */
        image = (double *)calloc(2,sizeof(double));
        aux = (double **)calloc(2,sizeof(double *));
        for(i=0 ;i<2;i++){
            aux[i] = (double *)calloc(3,sizeof(double));
        }
        if(aux == NULL || image == NULL){
            printf("Problems with memory in the DH_stimes function\n");
            exit(1);
        }
        /* Perform a matrix multiplication */
        for(n=0; n<s ; n++){
            /* If it is the first one, we just store the result */
            if(n == 0){
                A[0][0] = -2*a*eval_point[0]; A[0][1] = 1; A[1][0] = b; A[1][1] = 0;
                for(i=0;i<2;i++){
                    for(j=0;j<2;j++){
                        aux[i][j] = A[i][j];
                    }
                }
                /* Otherwise, we perform the matrix multiplication */
            }else{
                H(previous_x,previous_y,image,1);
                previous_x = image[0]; previous_y = image[1];
                A[0][0] = -2*a*image[0]; A[0][1] = 1; A[1][0] = b; A[1][1] = 0;
                matrix_mult(A,aux);
            }
        }
        
        /* Return the matrix in the desire pointer*/
        for(i=0;i<2;i++){
            for(j=0;j<2;j++){
                A[i][j] = aux[i][j];
            }
        }
        
        /* Free up memory */
        free(image);
        for(i = 0; i<2; i++){
            free(aux[i]);
        }
        free(aux);
    }
}

void matrix_mult(double **A, double **B){
    /*
     Function that performs the matrix multiplication AxB and store the result in B.
     Parameters:
     A: The first matrix.
     B: The second matrix.
     */
    
    int i, j, k;
    double C[2][2];
    
    /* Perform matrix multiplication */
    for (i = 0; i < 2; i++) {
        for (j = 0; j < 2; j++) {
            C[i][j] = 0;
            for (k = 0; k < 2; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
    
    /* Copy the result back to B*/
    for (i = 0; i < 2; i++) {
        for (j = 0; j < 2; j++) {
            B[i][j] = C[i][j];
        }
    }
}

void mat_times_vect(double **A, double *vect, double *sol){
    /*
     Function that perform the multiplication of a 2x2 matrix times a vector and store the results in sol.
     Parameters:
     A: The matrix to perform operations.
     vect: The vector to perform operations.
     sol: A pointer to which we assign the solution of the operation.
     */
    sol[0] = A[0][0]*vect[0] + A[0][1]*vect[1];
    sol[1] = A[1][0]*vect[0] + A[1][1]*vect[1];
}
