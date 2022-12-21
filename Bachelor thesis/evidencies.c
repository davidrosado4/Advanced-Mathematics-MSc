#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
double complex p(double complex z,double complex a);
double complex derp(double complex z,double complex a);
double modul(double complex z);
int newton2(double complex z0,double complex a,double tol);

int main(void){
     int i,j,result,contador,*vector,true,imprimir;
     FILE *comprovacio;
     double tol=1e-5,r,w,pi,theta,h=0;
     double complex a,ci;
     pi=acos(-1);
    vector=(int *)calloc(200,sizeof(int));
    if(vector==NULL){
        printf("Error al guardar memoria pel vector\n");
    }
    
    comprovacio=fopen("comprovacio.txt","w");
    if(comprovacio==NULL){

        printf("problemes amb el fitxer de sortida\n");
        exit(1);
    }
    w=0.05;
    while(w<=1.01){
        for(theta=0;theta<1;theta=theta+0.005){
            a=w*cexp(2*pi*I*theta);
            r=2.18149;
            contador=0;
            imprimir=0;
            h=0;
            for(i=0;i<=56;i++){
                ci=r+(-7+h)*I;
                h=h+0.25;
                printf("%f+%f*I\n",creal(ci),cimag(ci));
                result=newton2(ci,a,tol);
                vector[i]=result;
                true=0;
                /*printf("%d\n",result);*/
                for(j=0;j<i;j++){
                    if(vector[j]==result && result<=8 && true==0){
                        true=1;
                    }
                }
                if(true==0 && result!=10 && result!=20){
                    contador++;
                }
                if(contador==3 && imprimir==0){
                    /*fprintf(comprovacio,"Per a= %f + %f*I hem trobat totes les arrels amb %d cond in \n", creal(a),cimag(a),i+1);*/
                    imprimir=1;
                }
            }
            if(contador<3){
                fprintf(comprovacio,"Per a= %f + %f*I  NO hem trobat totes les arrels amb %d cond in \n", creal(a),cimag(a),i+1);
            }
        }
        w=w+0.05;
    }
   
    
    
    fclose(comprovacio);
    free(vector);
    return 0;
}

int newton2(double complex z0,double complex a,double tol){
     int maxiter;
    double diferencia;
    double complex z1;
    z1=z0;
    maxiter=100;
    diferencia=1;
    
    
    while(diferencia >tol && maxiter!=0 && modul(derp(z1,a)+p(z1,a))>tol && modul(p(z1,a))>tol){
        z0=z1;
        z1=z0-p(z0,a)/(derp(z0,a)+p(z0,a));
        maxiter--;
        diferencia=modul(z1-z0);
    }
    if(maxiter==0){
        return 20;
    }
     if(maxiter!=0  ){
        if(modul(z1)<tol){
            return 1;
        }
        if(modul(z1-1)<tol){
            return 2;
        }
        if(modul(z1-a)<tol){
            return 3;
        }
     }
    return 10;
}
double complex p(double complex z,double complex a){
    return z*(z-1)*(z-a);
}
double complex derp(double complex z,double complex a){
    return (z-1)*(z-a) + z*(z-a) + z*(z-1);
}
double modul(double complex z){
    double valor;
    valor= sqrt(creal(z)*creal(z) + cimag(z)*cimag(z));
    return valor;
}
