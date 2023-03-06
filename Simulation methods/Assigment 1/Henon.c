/* David Rosado Rodríguez*/
/* NIUB: 20194344*/

/* In this script we study the dynamics of the conservative Henon map */
/* We are using the library pgplot to make it interactive */


#include <stdio.h>
#include <cpgplot.h>
#include <math.h>
#include <stdlib.h>
#include <complex.h>

void dynamics_plot(void);
float *f_f(double x, double y);
double *f(double x, double y);
double *F(double x, double y);
double **Jacobian(double x, double y);
double **Jacobian_5(double x0, double y0);
double norm(double *vect);
double **inv(double **mat);
double *Newton(double x0, double y0);
double *solver_real(double tr, double det);
double complex *solver_complex(double tr, double det);


int main(void) {
	double *res,**jacob,tr,det;
	int i,j,k;
	/* Let us start plotting the dynamics of the Henon map */
	dynamics_plot();
	printf("You have clicked 20 times. Let us continue with the script\n\n");
	
	/*-----------------------------------------------------------------*/
	/*-------FIND THE HYPERBOLIC PERIODIC ORBIT OF PERIOD FIVE---------*/
	/*-----------------------------------------------------------------*/
	
	 
	/* Let us do Newton method using a list of initial condition */
	
	double initial_cond_hyper[10] = {0.35,-0.5,0.55,0.2,0.35,0.48,-0.38,0.35,-0.35,-0.3};
	printf("HYPERBOLIC PERIODIC ORBIT OF PERIOD FIVE\n\n\n");
	for(i=0;i<5;i++){
		res = Newton(initial_cond_hyper[2*i],initial_cond_hyper[2*i+1]);
		
		printf(" Fixed point by f^5:\n");
		printf("%le %le\n",res[0],res[1]);
		
		/* Finally, let us compute the eigenvalues of the Jacobian of f^5 in this points */
		/* This is, solving x^2 - tr*x + det=0, where tr is the trace and det the determinant */
		jacob = Jacobian_5(res[0],res[1]);
		/* We have to sum the identity in order to compute the Jacobian of f^5 instead of the Jacobian of f^5(x,y)-(x,y) */
		jacob[0][0] = jacob[0][0] + 1;
		jacob[1][1] = jacob[1][1] + 1;
		tr = jacob[0][0] + jacob[1][1];
		det = jacob[0][0]*jacob[1][1] - jacob[0][1]*jacob[1][0];
		if( (tr*tr-4*det) >= 0 ){
			double *vaps;
			vaps = solver_real(tr,det);
			printf("The eigenvalues are real. There are:\n");
			printf("%le %le\n\n",vaps[0],vaps[1]);
			free(vaps);
		}else{
			double complex *vaps;
			vaps = solver_complex(tr,det);
			printf("The eigenvalues are compelx. There are:\n");
			printf("%.2f %+.2fi\n", creal(vaps[0]), cimag(vaps[0]));
			printf("%.2f %+.2fi\n\n", creal(vaps[1]), cimag(vaps[1]));
			free(vaps);
		}
	}	
	
	/*-----------------------------------------------------------------*/
	/*-------FIND THE ELLIPTIC PERIODIC ORBIT OF PERIOD FIVE-----------*/
	/*-----------------------------------------------------------------*/
	
	
	printf("\n\n\n");
	double initial_cond_elipt[10] = {0,-0.6,0.6,-0.2,0.6,0.48,0,0.6,-0.6,0.2};
	printf("ELLIPTIC PERIODIC ORBIT OF PERIOD FIVE\n\n\n");
	for(i=0;i<5;i++){
		res = Newton(initial_cond_elipt[2*i],initial_cond_elipt[2*i+1]);
		
		printf(" Fixed point by f^5:\n");
		printf("%le %le\n",res[0],res[1]);
		
		/* Finally, let us compute the eigenvalues of the Jacobian of f^5 in this points */
		/* This is, solving x^2 - tr*x + det=0, where tr is the trace and det the determinant */
		jacob = Jacobian_5(res[0],res[1]);
		/* We have to sum the identity in order to compute the Jacobian of f^5 instead of the Jacobian of f^5(x,y)-(x,y) */
		jacob[0][0] = jacob[0][0] + 1;
		jacob[1][1] = jacob[1][1] + 1;
		tr = jacob[0][0] + jacob[1][1];
		det = jacob[0][0]*jacob[1][1] - jacob[0][1]*jacob[1][0];
		
		if( (tr*tr-4*det) >= 0 ){
			double *vaps;
			vaps = solver_real(tr,det);
			printf("The eigenvalues are real. There are:\n");
			printf("%le %le\n\n",vaps[0],vaps[1]);
			free(vaps);
		}else{
			double complex *vaps;
			vaps = solver_complex(tr,det);
			printf("The eigenvalues are compelx. There are:\n");
			printf("%.2f %+.2fi\n", creal(vaps[0]), cimag(vaps[0]));
			printf("%.2f %+.2fi\n\n", creal(vaps[1]), cimag(vaps[1]));
			free(vaps);
		}
	}	
	
	for(i=0;i<2;i++){
        free(jacob[i]);
	}
	free(jacob);
	free(res);
    return 0;
}


/* Function that plots in an interactive way the  dynamics of the Henon map */
void dynamics_plot(void){
	float xmin = -1 , xmax = 1, ymin = -1, ymax = 1; /* Define the x and y ranges for the plot */
	float x, y; /* Mouse clicks positions */
	char click;
	int iter=0,i, boolean,device;
	float *image;
	
	/* Initialize the plot */
	device = cpgopen("/XWIN");
	
	/* Set up the plot */
	cpgenv(xmin, xmax, ymin, ymax, 1, 1); 
	
	/* Labels for x-axis, y-axis, and top of plot */
	cpglab("X", "Y", "Dynamics of the Henon map, alpha = 1.33");
	
	boolean = cpgcurs(&x, &y, &click); /* Get mouse click*/
	
	while(boolean == 1 && iter <= 20){ /* Maximum 20 clicks, then we continue with the script*/
		for(i=0;i<1000;i++){ /* Iteration of the Henon map */
			image = f_f(x,y);
			cpgpt(1, &image[0], &image[1], -1);
			x = image[0]; y = image[1];
		}
		boolean = cpgcurs(&x, &y, &click); 
		iter++;
	}
	/* Close pgplot */
	cpgend();
	

	free(image);
}


/* Function that computes the image of Henon map (precision float) */
/*This is for use it in the dynamics plot, we do not need maximum precision*/
float *f_f(double x, double y){
	float *res;
	
	res = (float *)calloc(2,sizeof(float));
    if(res==NULL){
        printf("No memory\n");
        exit(1);
    }
    res[0] = x*cos(1.33) - y*sin(1.33) + x*x*sin(1.33);
    res[1] = x*sin(1.33) + y*cos(1.33) - x*x*cos(1.33);
    
    return res;
}


/* Function that computes the image of the Henon map (precision double) */
/*This is for use it in the Newton's method part, we need maximum precision*/
double *f(double x, double y){
	double *res;
	
	res = (double *)calloc(2,sizeof(double));
    if(res == NULL){
        printf("No memory\n");
        exit(1);
    }
    res[0] = x*cos(1.33) - y*sin(1.33) + x*x*sin(1.33);
    res[1] = x*sin(1.33) + y*cos(1.33) - x*x*cos(1.33);
    
    return res;
}


/* Function that for a given (x,y), computes the image of f^5(x,y)-(x,y) */
/* This is the function that we will apply Newton's method */
double *F(double x, double y){
	double *res, aux1, aux2;
	int i;
	
	res = (double *)calloc(2,sizeof(double));
    if(res == NULL){
        printf("No memory\n");
        exit(1);
    }
    /* Firstly, we compute f^5(x,y) */
    aux1 = x; aux2 = y;
    for(i=0;i<5;i++){
		res = f(aux1,aux2);
		aux1 = res[0];
		aux2 = res[1];
	}
	/* f^5(x,y) - (x,y) */
	res[0] = res[0] - x; res[1] = res[1] - y;
	return res;
}


/* Function that computes the Jacobian matrix of f */
double **Jacobian(double x, double y){
	double **mat;
	int i;
	
	mat = (double **)calloc(2,sizeof(double *));
    if(mat == NULL){
        printf("No memory\n");
        exit (1);
    }
    for(i=0;i<2;i++){
        mat[i] = (double *)calloc(2,sizeof(double));
		if(mat[i]==NULL){
			printf("No memory\n");
			exit (1);
		}
    }
	mat[0][0] = cos(1.33) + 2*x*sin(1.33);
	mat[1][0] = sin(1.33) - 2*x*cos(1.33);
	mat[0][1] = -sin(1.33);
	mat[1][1] = cos(1.33);
	return mat;
}


/*Function that computes the Jacobian matrix of f^5(x,y)-(x,y)*/
double **Jacobian_5(double x0, double y0){
	double **mat1,**mat2,**mat3,**mat4,**mat5, x, y, *aux;
	int i,j,k;
	double **mult_aux;
	
	mult_aux = (double **)calloc(2,sizeof(double *));
    if(mult_aux == NULL){
        printf("No memory\n");
        exit (1);
    }
    for(i=0;i<2;i++){
        mult_aux[i] = (double *)calloc(2,sizeof(double));
		if(mult_aux[i]==NULL){
			printf("No memory\n");
			exit (1);
		}
    }
	/* Firstly, we compute the Jacobian of f^5(x,y) using the chain rule */
	/* Notice that we only have to multiply 5 matrices */
	mat5 = Jacobian(x0,y0);
	aux = f(x0,y0);
	mat4 = Jacobian(aux[0],aux[1]);
	aux = f(aux[0],aux[1]);
	mat3 = Jacobian(aux[0],aux[1]);
	aux = f(aux[0],aux[1]);
	mat2 = Jacobian(aux[0],aux[1]);
	aux = f(aux[0],aux[1]);
	mat1 = Jacobian(aux[0],aux[1]);
	
	for(i=0;i<2;i++){
		for(j=0;j<2;j++){
			for(k=0;k<2;k++){
				mult_aux[i][j]+=mat1[i][k] * mat2[k][j];
			}
		}
	}
	for(i=0;i<2;i++){
		for(j=0;j<2;j++){
			mat1[i][j]=0;
			for(k=0;k<2;k++){
				mat1[i][j]+=mult_aux[i][k] * mat3[k][j];
			}
		}
	}
	for(i=0;i<2;i++){
		for(j=0;j<2;j++){
			mult_aux[i][j]=0;
			for(k=0;k<2;k++){
				mult_aux[i][j]+=mat1[i][k] * mat4[k][j];
			}
		}
	}
	for(i=0;i<2;i++){
		for(j=0;j<2;j++){
			mat1[i][j]=0;
			for(k=0;k<2;k++){
				mat1[i][j]+=mult_aux[i][k] * mat5[k][j];
			}
		}
	}
	/* Finally, we substract the identity matrix since the Jacobian of (x,y) is I */
	mat1[0][0] = mat1[0][0] -1;
	mat1[1][1] = mat1[1][1] -1;
	for(i=0;i<2;i++){
        free(mat2[i]);
        free(mat3[i]);
        free(mat4[i]);
        free(mat5[i]);
        free(mult_aux[i]);
    }
    free(mat2);free(mat3);free(mat4);free(mat5);free(mult_aux);
    free(aux);
    /* Return the Jacobian of f^5(x,y) - (x,y) */
	return mat1;
}


/* Newton's method */
double *Newton(double x0, double y0){
	double diference, *res, **jacob,*vect,det,*image;
	int i,iter;
	iter=0;
	
	res = (double *)calloc(2,sizeof(double));
    if(res == NULL){
        printf("No memory\n");
        exit(1);
    }
    vect = (double *)calloc(2,sizeof(double));
    if(vect == NULL){
        printf("No memory\n");
        exit(1);
    }
	diference = 1; /* For entering the loop*/
	/*Newton's method */
	while(diference>1e-15 && norm(F(x0,y0))>1e-15 && iter<=50){
		/*Computes the Jacobian and its inverse*/
		jacob = Jacobian_5(x0,y0);
		jacob = inv(jacob);
		image = F(x0,y0);
		/* Newton iteration */
		res[0] = x0 - (jacob[0][0]*image[0] + jacob[0][1]*image[1]);
		res[1] = y0 - (jacob[1][0]*image[0] + jacob[1][1]*image[1]);
		/*Compute the difference between iterations*/
		vect[0] = res[0]-x0;
		vect[1] = res[1]-y0;
		diference = norm(vect);
		/*Update variables*/
		x0 = res[0];
		y0 = res[1];
		iter++;
	}
	printf("Newton's method has done %d iterations.", iter);
	free(vect); free(image);
	for(i=0;i<2;i++){
        free(jacob[i]);
	}
	free(jacob);
	return res;
}


/*Function that computes the euclidian norm of a vector*/
double norm(double *vect){
	double res;
	
	res = sqrt(vect[0]*vect[0] + vect[1]*vect[1]);
	return res;
}


/*Function that computes the inverse of a 2x2 matrix */
double **inv(double **mat){
	double aux,det;
	det = mat[0][0] * mat[1][1] - mat[0][1]*mat[1][0];
	if(fabs(det)<1e-12){
		printf("The determinant is 0. We have to stop the script\n");
		exit(1);
	}
	aux = mat[0][0];
	mat[0][0] = mat[1][1]/det;
	mat[1][1] = aux/det;
	mat[1][0] = -mat[1][0]/det;
	mat[0][1] = -mat[0][1]/det;
	return mat;
}


/* Function that solves the equation x^2 - tr*x + det=0 to find the eigenvalues of the Jacobian of F */
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


/* Function that solves the equation x^2 - tr*x + det=0 to find the eigenvalues of the Jacobian of F */
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
		
	
	








