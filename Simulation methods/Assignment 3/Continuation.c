/* David Rosado Rodríguez */
/* NIUB: 20194344 */


#include "rkf45.h" /* Includes the given Runge-Kutta */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>

#define w_square1  2
#define w_square2  0.5
#define w_forwards 0.7081
#define w_backwards 1.40
#define eps 1e-2
#define pi acos(-1)

/* Function definition */
void gaussEliminationLS(int m, int n, double **a, double *b, double *x);
void ode_eq1(double t, double *sol, int n, double *vf);
void ode_eq2(double t, double *sol, int n, double *vf);
int Newton_ex1_ex2(double x0, double y0, double *res, int choice);
double norm_2(double *vect);
double norm_3(double *vect);
void plot_orbit(double *sol,double at, int n, double ah, double tol, double atf, void *ode, char *file_name);
void ode_eq3(double t, double *sol, int n, double *vf);
int Newton_ex3_continuation(double x0, double y0, double w0, double delta, double x_new, double y_new, double w_new, double *res);
int Newton_ex3_firstStep(double x0, double y0, double *res, int choice);
void continuation(double x0, double y0, int choice, double total_file_size);
void ode_eq_forwards(double t, double *sol, int n, double *vf);
void ode_eq_backrwards(double t, double *sol, int n, double *vf);
double complex *solver_complex(double tr, double det);
double *solver_real(double tr, double det);

int main(void){
	double *res, previous_1,previous_2;
	int iterations;
	
	/* Allocate memory to store the result */
	res = (double *)calloc(2,sizeof(double));
    if(res == NULL){
        printf("No memory\n");
        exit(1);
    }
	/*-----------------------EXERCISE 1--------------------------------*/
    
    /* Compute periodic orbit of period 2π near the origin */
    printf("\n\n\033[1;34m---------------------EXERCISE 1--------------------------\033[0m\n\n");
    
    iterations = Newton_ex1_ex2(0,0,res,1);
    printf("\033[1;1mNewton's method has done %d iterations\033[0m\n", iterations);
    printf("\033[1;31mThe fixed point by the Poincaré map is\033[0m: (%le ,%le)\n",res[0],res[1]);
    previous_1 = res[0]; previous_2 = res[1]; /* Store the points */
    
    /*-----------------------EXERCISE 2--------------------------------*/
    res[0] = 0; res[1] = 0;
    
    /* Compute periodic orbit of period 2π near the origin*/
    printf("\n\n\033[1;34m---------------------EXERCISE 2--------------------------\033[0m\n\n");
    
    iterations = Newton_ex1_ex2(0,0,res,2);
    printf("\033[1;1mNewton's method has done %d iterations\033[0m\n", iterations);
    printf("\033[1;31mThe fixed point by the Poincaré map is\033[0m: (%le ,%le)\n",res[0],res[1]);
    
    
    /*-----------------------EXERCISE 2--------------------------------*/
    printf("\n\n\033[1;34m---------------------EXERCISE 3--------------------------\033[0m\n\n");
    /* Printing 200 KB in each file. Change this variable for another option */
    double KB_MAX = 200;
    
    /* Forwards continuation */
    continuation(res[0],res[1],1,KB_MAX);
    
    /* Backwards continuation */
    continuation(previous_1,previous_2,2,KB_MAX);
    
	printf("\n\n"); free(res);
	return 0;
}

/* FUNCTIONS FOR EXERCISE 1 AND 2 */

/* Function that computes the vector field of the ODE we want to solve for w_square1 */
void ode_eq1(double t, double *sol, int n, double *vf){
	
	/* Assign values */
	vf[0] = sol[1];
	vf[1] = -w_square1*sin(sol[0]) + eps*sin(t);
	vf[2] = sol[4];
	vf[3] = sol[5];
	vf[4] = -w_square1*cos(sol[0])*sol[2];
	vf[5] = -w_square1*cos(sol[0])*sol[3];
}

/* Function that computes the vector field of the ODE we want to solve for w_sqaure2 */
void ode_eq2(double t, double *sol, int n, double *vf){
	
	/* Assign values */
	vf[0] = sol[1];
	vf[1] = -w_square2*sin(sol[0]) + eps*sin(t);
	vf[2] = sol[4];
	vf[3] = sol[5];
	vf[4] = -w_square2*cos(sol[0])*sol[2];
	vf[5] = -w_square2*cos(sol[0])*sol[3];
}

/* Function that performs a Newton method of the Poincaré map */
int Newton_ex1_ex2(double x0, double y0, double *res, int choice){
	
	/* If choice == 1, we compute exercise 1 (w = sqrt(2)), otherwise, we compute exercise 2 */
	
	
	double *sol,*vect,*aux,*b,**Jacob;
	int iter,rungekutta, n = 6,i;
	double diff;
	/* Variables for the integration method */
	double at = 0,ah = 1e-5,tol = 1e-15,atf = 2 * pi, *vf;
	void *ode;
	
	/* Allocate memory */
		
	/* Store the solution of the ode */
	sol = (double *)calloc(n,sizeof(double));
	/* Store the vectorfield */
	vf = (double *)calloc(n,sizeof(double));
	/* Auxilar variables */
	vect = (double *)calloc(2,sizeof(double));
	aux = (double *)calloc(2,sizeof(double));
	/* Variables to solve the linear system */
	b = (double *)calloc(2,sizeof(double));
	Jacob = (double **)calloc(2, sizeof(double *));
	
	for (i = 0; i < 2; i++) {
		Jacob[i] = (double *)calloc(3, sizeof(double));
	}
	
	if(sol == NULL || vf == NULL || vect == NULL || aux == NULL || b == NULL  || Jacob == NULL){
		printf("No memory\n");
		exit(1);
	}
	
	/* Initialize variables */
	sol[0]= x0; sol[1] = y0; sol[2] = 1; sol[3] = 0; sol[4]= 0; sol[5] = 1;
	diff = 1; iter = 0; /* For enter in the loop and counts iterations */
	
	
	/* Newton loop */
	while(diff>1e-15 && iter<=50){
		
		/* Initiliaze variables for the Runge-Kutta */
		sol[2] = 1; sol[3] = 0; sol[4]= 0; sol[5] = 1; at = 0;
		
		/* Choose the vectorfield*/
		if(choice == 1){
			ode = &ode_eq1;
		}else{
			ode = &ode_eq2;
		}
		/* Numerical integration using Runge-Kutta */
		do{
			rungekutta = rkf45(&at, sol, n, &ah, 1, tol, &atf, NULL, ode);
		}while(rungekutta != 1);
		
		/* Perform the Newton method using a linear system */
		
		/* Prepare the Jacobian matrix */
		Jacob[0][0] = sol[2]-1; Jacob[0][1] = sol[3]; Jacob[1][0] = sol[4]; Jacob[1][1] = sol[5]-1;
		
		/* Prepare the vector b to solve Jx=b */
		b[0] = -sol[0] + aux[0]; b[1] = -sol[1] + aux[1];
		
		/* Gaussian elimination and store the solution in res*/
		gaussEliminationLS(2,3,Jacob,b,res);
		
		/* Final solution */
		res[0] +=aux[0]; res[1]+=aux[1];
		
		/*Compute the difference between iterations*/
		vect[0] = res[0]-aux[0];
		vect[1] = res[1]-aux[1];
		diff = norm_2(vect);
		
		/*Update variables*/
		sol[0] = res[0]; sol[1] = res[1];
		aux[0] = sol[0]; aux[1] = sol[1];
		iter++;	
	}
	
	/* Plot the final orbit and find the eigenvalues of the Jacobian of Pat that point */
	if(choice == 1){
		double det,tr;
		plot_orbit(sol, at, n, ah, tol, atf, ode, "periodic_ex1.txt");
		
		/* Find the determinant of the matrix */
		det = sol[2]*sol[5] - sol[3]*sol[4];
		tr = sol[2] + sol[5];
		
		/* Find the eigenvalues */
		if( (tr*tr-4*det) >= 0 ){
			double *vaps;
			vaps = solver_real(tr,det);
			printf("The eigenvalues are real. They are:\n");
			printf("%le %le\n\n",vaps[0],vaps[1]);
			printf("The modulus of the eigenvalues are (%le, %le)\n",fabs(vaps[0]),fabs(vaps[1]));
			free(vaps);
		}else{
			double complex *vaps;
			vaps = solver_complex(tr,det);
			printf("The eigenvalues are compelx. They are:\n");
			printf("%.2f %+.2fi\n", creal(vaps[0]), cimag(vaps[0]));
			printf("%.2f %+.2fi\n\n", creal(vaps[1]), cimag(vaps[1]));
			
			printf("The modulus of the eigenvalues are (%le, %le)\n", cabs(vaps[0]),cabs(vaps[1]));
			free(vaps);
		}
	}else{
		double det,tr;
		plot_orbit(sol, at, n, ah, tol, atf, ode, "periodic_ex2.txt");
		/* Find the determinant of the matrix */
		det = sol[2]*sol[5] - sol[3]*sol[4];
		tr = sol[2] + sol[5];
		
		/* Find the eigenvalues */
		if( (tr*tr-4*det) >= 0 ){
			double *vaps;
			vaps = solver_real(tr,det);
			printf("The eigenvalues are real. They are:\n");
			printf("%le %le\n\n",vaps[0],vaps[1]);
			printf("The modulus of the eigenvalues are (%le, %le)\n",fabs(vaps[0]),fabs(vaps[1]));
			free(vaps);
		}else{
			double complex *vaps;
			vaps = solver_complex(tr,det);
			printf("The eigenvalues are compelx. They are:\n");
			printf("%.2f %+.2fi\n", creal(vaps[0]), cimag(vaps[0]));
			printf("%.2f %+.2fi\n\n", creal(vaps[1]), cimag(vaps[1]));
			
			printf("The modulus of the eigenvalues are (%le, %le)\n", cabs(vaps[0]),cabs(vaps[1]));
			free(vaps);
		}
	}
	
	/* Free up memory */
	free(vect); free(sol); free(vf);  free(aux); free(b);
	for (int i = 0; i < 2; i++) {
		free(Jacob[i]);
	}
	free(Jacob);
	
	return iter;
}

/* FUNCTIONS FOR EXERCISE 3*/

/* Function that computes the vector field of the ODE we want to solve (for forward continuation) the fist step */
void ode_eq_forwards(double t, double *sol, int n, double *vf){
	
	/* Assign values */
	vf[0] = sol[1];
	vf[1] = -w_forwards*w_forwards*sin(sol[0]) + eps*sin(t);
	vf[2] = sol[4];
	vf[3] = sol[5];
	vf[4] = -w_forwards*w_forwards*cos(sol[0])*sol[2];
	vf[5] = -w_forwards*w_forwards*cos(sol[0])*sol[3];
}

/* Function that computes the vector field of the ODE we want to solve (for backwards continuation) the first step */
void ode_eq_backrwards(double t, double *sol, int n, double *vf){
	
	/* Assign values */
	vf[0] = sol[1];
	vf[1] = -w_backwards*w_backwards*sin(sol[0]) + eps*sin(t);
	vf[2] = sol[4];
	vf[3] = sol[5];
	vf[4] = -w_backwards*w_backwards*cos(sol[0])*sol[2];
	vf[5] = -w_backwards*w_backwards*cos(sol[0])*sol[3];
}

/* Function that computes the vector field of the ODE we want to solve for the continuation method */
void ode_eq3(double t, double *sol, int n, double *vf){
	/* Assign values */
	vf[0] = sol[1];
	vf[1] = -sol[2]*sol[2]*sin(sol[0]) + eps*sin(t);
	vf[2] = 0;
	vf[3] = sol[6];
	vf[4] = sol[7];
	vf[5] = sol[8];
	vf[6] = -sol[2]*sol[2]*cos(sol[0])*sol[3] - 2*sol[2]*sin(sol[0])*sol[9];
	vf[7] = -sol[2]*sol[2]*cos(sol[0])*sol[4] - 2*sol[2]*sin(sol[0])*sol[10];
	vf[8] = -sol[2]*sol[2]*cos(sol[0])*sol[5] - 2*sol[2]*sin(sol[0])*sol[11];
	vf[9] = 0;
	vf[10] = 0;
	vf[11] = 0;
}

/* Function that performs a Newton method of the Poincaré map for the continuation method*/
int Newton_ex3_continuation(double x0, double y0, double w0, double delta, double x_new, double y_new, double w_new, double *res){
	
	double *sol,*vect,*aux,*b,**Jacob;
	int iter,rungekutta, n = 12,i;
	double diff;
	/* Variables for the integration method */
	double at = 0,ah = 1e-5,tol = 1e-15,atf = 2 * pi, *vf;
	void *ode;
	
	/* Allocate memory */
		
	/* Store the solution of the ode */
	sol = (double *)calloc(n,sizeof(double));
	/* Store the vectorfield */
	vf = (double *)calloc(n,sizeof(double));
	/* Auxilar variables */
	vect = (double *)calloc(3,sizeof(double));
	aux = (double *)calloc(3,sizeof(double));
	/* Variables to solve the linear system */
	b = (double *)calloc(3,sizeof(double));
	Jacob = (double **)calloc(3, sizeof(double *));
	
	for (i = 0; i < 3; i++) {
		Jacob[i] = (double *)calloc(4, sizeof(double));
	}
	
	if(sol == NULL || vf == NULL || vect == NULL || aux == NULL || b == NULL || Jacob == NULL){
		printf("No memory\n");
		exit(1);
	}
	
	/* Initialize variables */
	sol[0]= x_new; sol[1] = y_new; sol[2] = w_new; 
	aux[0] = x_new; aux[1] = y_new; aux[2] = w_new; 
	diff = 1; iter = 0; /* For enter in the loop and counts iterations */
	
	/* Newton loop */
	while(diff>1e-15 && iter<=50){
		
		/* Initiliaze variables for the Runge-Kutta */
		for(i=3;i<n;i++){
			sol[i] = 0;
		}
		sol[3] = 1; sol[7] = 1;
		ode = &ode_eq3;
		at = 0;
		/* Numerical integration using Runge-Kutta */
		do{
			rungekutta = rkf45(&at, sol, n, &ah, 1, tol, &atf, NULL, ode);
		}while(rungekutta != 1);
		/* Perform the Newton method using a linear system */
		/* Prepare the Jacobian matrix */
		Jacob[0][0] = sol[3]-1; Jacob[0][1] = sol[4]; Jacob[0][2] = sol[5];
		Jacob[1][0] = sol[6]; Jacob[1][1] = sol[7]-1; Jacob[1][2] = sol[8];
		Jacob[2][0] = 2*(aux[0]-x0); Jacob[2][1] = 2*(aux[1]-y0); Jacob[2][2] = 2*(aux[2]-w0); 
		
		/* Prepare the vector b to solve Jx=b */
		b[0] = -sol[0] + aux[0]; b[1] = -sol[1]+aux[1]; b[2] = -(aux[0]-x0)*(aux[0]-x0) - (aux[1]-y0)*(aux[1]-y0) - (aux[2]-w0)*(aux[2]-w0)+ delta*delta;
		
		/* Gaussian elimination and store the solution in res*/
		gaussEliminationLS(3,4,Jacob,b,res);
		
		/* Final solution */
		res[0] +=aux[0]; res[1]+=aux[1]; res[2]+=aux[2];
		
		/*Compute the difference between iterations*/
		vect[0] = res[0]-aux[0];
		vect[1] = res[1]-aux[1];
		vect[2] = res[2]-aux[2];
		diff = norm_3(vect);
		
		/*Update variables*/
		sol[0] = res[0]; sol[1] = res[1]; sol[2] = res[2];
		aux[0] = sol[0]; aux[1] = sol[1]; aux[2] = sol[2];
		iter++;	
	}
	/* Free up memory */
	free(vect); free(sol); free(vf);  free(aux); free(b);
	for (int i = 0; i < 3; i++) {
		free(Jacob[i]);
	}
	free(Jacob);
	
	return iter;
}

/* Function for computing the first step in the continuation method */
int Newton_ex3_firstStep(double x0, double y0, double *res, int choice){
	
	double *sol,*vect,*aux,*b,**Jacob;
	int iter,rungekutta, n = 6,i;
	double diff;
	/* Variables for the integration method */
	double at = 0,ah = 1e-5,tol = 1e-15,atf = 2 * pi, *vf;
	void *ode;
	
	/* Allocate memory */
		
	/* Store the solution of the ode */
	sol = (double *)calloc(n,sizeof(double));
	/* Store the vectorfield */
	vf = (double *)calloc(n,sizeof(double));
	/* Auxilar variables */
	vect = (double *)calloc(2,sizeof(double));
	aux = (double *)calloc(2,sizeof(double));
	/* Variables to solve the linear system */
	b = (double *)calloc(2,sizeof(double));
	Jacob = (double **)calloc(2, sizeof(double *));
	
	for (i = 0; i < 2; i++) {
		Jacob[i] = (double *)calloc(3, sizeof(double));
	}
	
	if(sol == NULL || vf == NULL || vect == NULL || aux == NULL || b == NULL || Jacob == NULL){
		printf("No memory\n");
		exit(1);
	}
	
	/* Initialize variables */
	sol[0]= x0; sol[1] = y0; sol[2] = 1; sol[3] = 0; sol[4]= 0; sol[5] = 1;
	diff = 1; iter = 0; /* For enter in the loop and counts iterations */
	
	
	/* Newton loop */
	while(diff>1e-15 && iter<=50){
		
		/* Initiliaze variables for the Runge-Kutta */
		sol[2] = 1; sol[3] = 0; sol[4]= 0; sol[5] = 1; at = 0;
		if(choice == 1){
			ode = &ode_eq_forwards;
		}else{
			ode = &ode_eq_backrwards;
		}
		
		/* Numerical integration using Runge-Kutta */
		do{
			rungekutta = rkf45(&at, sol, n, &ah, 1, tol, &atf, NULL, ode);
		}while(rungekutta != 1);
		
		/* Perform the Newton method using a linear system */
		
		/* Prepare the Jacobian matrix */
		Jacob[0][0] = sol[2]-1; Jacob[0][1] = sol[3]; Jacob[1][0] = sol[4]; Jacob[1][1] = sol[5]-1;
		
		/* Prepare the vector b to solve Jx=b */
		b[0] = -sol[0] + aux[0]; b[1] = -sol[1] + aux[1];
		
		/* Gaussian elimination and store the solution in res*/
		gaussEliminationLS(2,3,Jacob,b,res);
		
		/* Final solution */
		res[0] +=aux[0]; res[1]+=aux[1];
		
		/*Compute the difference between iterations*/
		vect[0] = res[0]-aux[0];
		vect[1] = res[1]-aux[1];
		diff = norm_2(vect);
		
		/*Update variables*/
		sol[0] = res[0]; sol[1] = res[1];
		aux[0] = sol[0]; aux[1] = sol[1];
		iter++;	
	}
	/* Free up memory */
	free(vect); free(sol); free(vf);  free(aux); free(b); 
	for (int i = 0; i < 2; i++) {
		free(Jacob[i]);
	}
	free(Jacob);
	
	return iter;
}

/* Continuation method */
void continuation(double x0, double y0, int choice,double total_file_size){
	int iterations;
	double *res,*res_first, delta = 0.001, *v, *v_ant, prod_esc,*aux,*aux_new, norm_vect;
	FILE *fit;
	long size;
	long long Size;
	double percentage;
	
	/* Create a file to plot the continuation method*/
	if(choice == 1){
		fit=fopen("continuation_forwards.txt","w");
	}else{
		fit=fopen("continuation_backwards.txt","w");
	}
	if(fit == NULL){
		printf("Problems with the file\n");
		exit(1);
	}
	/* Allocate memory */
	res = (double *)calloc(3,sizeof(double));
	res_first = (double *)calloc(2,sizeof(double));
	v = (double *)calloc(3,sizeof(double));
	v_ant = (double *)calloc(3,sizeof(double));
	aux = (double *)calloc(3,sizeof(double));
	aux_new = (double *)calloc(3,sizeof(double));
    if(res == NULL || v == NULL || v_ant == NULL || res_first== NULL || aux == NULL || aux_new == NULL){
        printf("No memory\n");
        exit(1);
    }
	/* Perform a first Newton step to start the continuation method */
	if(choice == 1){
		iterations = Newton_ex3_firstStep(x0,y0,res_first,1);
	}else{
		iterations = Newton_ex3_firstStep(x0,y0,res_first,2);
	}
	
	/* Compute v and normalize for computing the next point*/
	v[0] = res_first[0]-x0; v[1] = res_first[1]-y0; v[2] = delta;
	norm_vect = norm_3(v);
	v[0] = v[0]/norm_vect; v[1] = v[1]/norm_vect; v[2] = v[2]/ norm_vect;
	
	/* Prepare the previous and the next point*/
	if(choice == 1){
		aux[2] = 1/sqrt(2);
		aux_new[2] = aux[2] + delta*v[2];
	}else{
		v[2]= -v[2];
		aux[2] = sqrt(2);
		aux_new[2] = aux[2] + delta*v[2];
	}
	aux[0] = x0; aux[1] = y0; 
	aux_new[0] = aux[0] + delta*v[0]; aux_new[1] = aux[1] + delta*v[1]; 
	
	/* Start the continuation*/
	if(choice == 1){
		printf("Computing the forwards continuation...\n");
		printf("This part may take some time (about 2 minutes in a MacOS M1)\n");
		while(aux_new[2]<=sqrt(2)){
			v_ant[0] = v[0]; v_ant[1] = v[1]; v_ant[2] = v[2];
			
			/* Newton iteration */
			iterations = Newton_ex3_continuation(aux[0],aux[1],aux[2],delta,aux_new[0],aux_new[1],aux_new[2],res);
			fprintf(fit,"%le %le\n", res[2],res[1]); /* axis x-->w, axis y-->y */
			
			/* Prepare next vector v */
			v[0] = res[0]-aux[0]; v[1]= res[1]- aux[1]; v[2] = res[2]- aux[2];
			norm_vect = norm_3(v);
			v[0] = v[0]/norm_vect; v[1] = v[1]/norm_vect; v[2] = v[2]/ norm_vect;
			
			/* Check if we have to change the direction */
			prod_esc = v[0]*v_ant[0] + v[1]*v_ant[1] + v[2]*v_ant[2];
			if(prod_esc < 0){
				v[0] = -v[0]; v[1] = -v[1]; v[2] = -v[2];
			}
			
			/* Compute the new delta */
			if(iterations < 10){
				if(delta<=0.001){
					delta = delta * 1.5;
				}
			}else{
				delta = delta/2;
			}
			/* Prepare the new point */
			aux[0] = res[0]; aux[1] = res[1]; aux[2]=res[2];
			aux_new[0] = aux[0] + delta*v[0]; aux_new[1] = aux[1] + delta*v[1]; aux_new[2] = aux[2] + delta*v[2];
			res[0] = 0; res[1] = 0; res[2] = 0;
			
			/* Control de size of the file to finish the program */
			
			fseek(fit, 0, SEEK_END); /* move file pointer to end of file */
			size = ftell(fit); /* get current position of file pointer (size of file in bytes) */
			size /= 1024; /* convert to kilobytes */
			if(size >= total_file_size){
				printf("\rPercentage done: %.2f %%",100.);
				printf("\n%.1f KB memory full for the continuation forwards. Last omega: %le\n",total_file_size, aux_new[2]);
				break;
			}
			
			/* Print controls... */
			Size = ftell(fit);
			percentage = (double)Size/1024 *100/total_file_size;
			printf("\rPercentage done: %.2f %%",percentage);
			fflush(stdout);
		}
	}else{
		printf("\n\nComputing the backwards continuation...\n");
		printf("This part may take some time (about 2 minutes in a MacOS M1)\n");
		while(aux_new[2]>=1/sqrt(2)){
			v_ant[0] = v[0]; v_ant[1] = v[1]; v_ant[2] = v[2];
			
			/* Newton iteration */
			iterations = Newton_ex3_continuation(aux[0],aux[1],aux[2],delta,aux_new[0],aux_new[1],aux_new[2],res);
			fprintf(fit,"%le %le\n", res[2],res[1]); /* axis x-->w, axis y-->y */
			
			/* Prepare next vector v */
			v[0] = res[0]-aux[0]; v[1]= res[1]- aux[1]; v[2] = res[2]- aux[2];
			norm_vect = norm_3(v);
			v[0] = v[0]/norm_vect; v[1] = v[1]/norm_vect; v[2] = v[2]/ norm_vect;
			
			/* Check if we have to change the direction */
			prod_esc = v[0]*v_ant[0] + v[1]*v_ant[1] + v[2]*v_ant[2];
			if(prod_esc < 0){
				v[0] = -v[0]; v[1] = -v[1]; v[2] = -v[2];
			}
			
			/* Compute the new delta */
			if(iterations < 10){
				if(delta<=0.001){
					delta = delta * 1.5;
				}
			}else{
				delta = delta/2;
			}
			/* Prepare the new point */
			aux[0] = res[0]; aux[1] = res[1]; aux[2]=res[2];
			aux_new[0] = aux[0] + delta*v[0]; aux_new[1] = aux[1] + delta*v[1]; aux_new[2] = aux[2] + delta*v[2];
			res[0] = 0; res[1] = 0; res[2] = 0;
			
			/* Control de size of the file to finish the program */
			
			fseek(fit, 0, SEEK_END); /* move file pointer to end of file */
			size = ftell(fit); /* get current position of file pointer (size of file in bytes) */
			size /= 1024; /* convert to kilobytes */
			if(size >= total_file_size){
				printf("\rPercentage done: %.2f %%",100.);
				printf("\n%.1f KB memory full for the continuation backwards. Last omega: %le\n",total_file_size, aux_new[2]);
				break;
			}
			/* Print controls... */
			Size = ftell(fit);
			percentage = (double)Size/1024 *100/total_file_size;
			printf("\rPercentage done: %.2f %%",percentage);
			fflush(stdout);
		}
		
	}
	/* Free up memory */
	free(res); free(res_first); free(v); free(aux); free(v_ant); free(aux_new); fclose(fit);
}

/* GENERAL FUNCTIONS */

/*Function that computes the euclidian norm of a vector in R^2*/
double norm_2(double *vect){
	double res;
	
	res = sqrt(vect[0]*vect[0] + vect[1]*vect[1]);
	return res;
}

/*Function that computes the euclidian norm of a vector in R^3*/
double norm_3(double *vect){
	double res;
	
	res = sqrt(vect[0]*vect[0] + vect[1]*vect[1] + vect[2]*vect[2]);
	return res;
}

/* Function that performs Gaussian elimination to solve a linear system */
void gaussEliminationLS(int m, int n, double **a, double *b, double *x) {
    int i,j,k;
    /* Append the vector to the matrix */
    for(i=0; i<m;i++){
		a[i][n-1] = b[i];
	}
	/* Start Gauss */
    for(i=0; i<m-1; i++) {
        /* Partial Pivoting */
        int max_row = i;
        for(k=i+1; k<m; k++) {
            /*If diagonal element(absolute value) is smaller than any of the terms below it */
            if(fabs(a[i][i]) < fabs(a[k][i])) {
                /*Track the row interchanges */
                max_row = k;
            }
        }
        if(max_row != i) {
            /* Swap the rows */
            for(j=0; j<n; j++) {                
                double temp = a[i][j];
                a[i][j] = a[max_row][j];
                a[max_row][j] = temp;
            }
        }
        /* Begin Gauss Elimination */
        for(k=i+1; k<m; k++) {
            double term = a[k][i] / a[i][i];
            for(j=0; j<n; j++) {
                a[k][j] = a[k][j] - term*a[i][j];
            }
        }  
    }
    /* Begin Back-substitution */
    for(i=m-1; i>=0; i--) {
        x[i] = a[i][n-1];
        for(j=i+1; j<n-1; j++) {
            x[i] = x[i] - a[i][j]*x[j];
        }
        x[i] = x[i] / a[i][i];
    }
}

/* Function that plots an orbit using Runge-Kutta integration method*/
void plot_orbit(double *sol,double at, int n, double ah, double tol, double atf, void *ode, char *file_name){
	int rungekutta;
	/* Create a file */	
	FILE *fit;
	fit=fopen(file_name, "w");
	if(fit == NULL){
		printf("Failed to open file: %s\n", file_name);
		printf("Problems with the file\n");
		exit(1);
	}
	/* Perform an integration method to compute the periodic orbit and print points */
	sol[2] = 1; sol[3] = 0; sol[4]= 0; sol[5] = 1; at = 0;
	do{
		rungekutta = rkf45(&at, sol, n, &ah, 1, tol, &atf, NULL, ode);
		fprintf(fit,"%le %le\n",sol[0],sol[1]);
	}while(rungekutta != 1);
	/* Close the file*/
	fclose(fit); 
}

/* Function that solves the equation x^2 - tr*x + det=0 to find the eigenvalues of the Jacobian of P */
/* Real case */
double *solver_real(double tr, double det){
	double *vap;
	
	vap = (double *)calloc(2,sizeof(double));
    if(vap == NULL){
        printf("No memory\n");
        exit(1);
    }
	vap[0] = (tr + sqrt(tr*tr - 4*det))/2;
	vap[1] = (tr - sqrt(tr*tr - 4*det))/2;
	return vap;
}

/* Function that solves the equation x^2 - tr*x + det=0 to find the eigenvalues of the Jacobian of P */
/* Complex case */
double complex *solver_complex(double tr, double det){
	double complex *vap;
	
	vap = (double complex *)calloc(2,sizeof(double complex));
    if(vap == NULL){
        printf("No memory\n");
        exit(1);
    }
	vap[0] = (tr + csqrt(tr*tr - 4*det))/2;
	vap[1] = (tr - csqrt(tr*tr - 4*det))/2;
	return vap;
}







