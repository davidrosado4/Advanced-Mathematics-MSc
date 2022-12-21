#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
void newton(double complex z0,double tol, FILE *fitxer1,FILE *fitxer2,FILE *fitxer3,FILE *fitxer4);
double complex p(double complex z);
double complex derp(double complex z);
double modul(double complex z);

int main(void){
     FILE *fitxer1,*fitxer2,*fitxer3,*fitxer4;
     double tol=1e-5,h,m;
     double complex z0;
    fitxer1=fopen("fitxer1.txt","w");
    if(fitxer1==NULL){

        printf("problemes amb el fitxer de sortida\n");
        exit(1);
    }
    fitxer2=fopen("fitxer2.txt","w");
    if(fitxer2==NULL){

        printf("problemes amb el fitxer de sortida\n");
        exit(1);
    }
    fitxer3=fopen("fitxer3.txt","w");
    if(fitxer3==NULL){

        printf("problemes amb el fitxer de sortida\n");
        exit(1);
    }
    fitxer4=fopen("fitxer4.txt","w");
    if(fitxer4==NULL){

        printf("problemes amb el fitxer de sortida\n");
        exit(1);
    }
    m=2;
    h=-2;
    while(m>-2){
        h=-2;
        m=m-0.005;
        while(h<2){
            h=h+0.005;
            z0=h+m*I;
            newton(z0,tol,fitxer1,fitxer2,fitxer3,fitxer4);
        }
    }
    
    fclose(fitxer1);
    fclose(fitxer2);
    fclose(fitxer4);
    fclose(fitxer3);
    return 0;
}
void newton(double complex z0,double tol, FILE *fitxer1,FILE *fitxer2,FILE *fitxer3,FILE *fitxer4){
    int maxiter;
    double diferencia;
    double complex condinicial,z1;
    z1=z0;
    condinicial=z0;
    maxiter=50;
    diferencia=1;
    
    
    while(diferencia >tol && maxiter!=0 && modul(derp(z1))>tol && modul(p(z1))>tol){
        z0=z1;
        z1=z0-p(z0)/derp(z0);
        maxiter--;
        diferencia=modul(z1);
    }
    if(maxiter==0){
        fprintf(fitxer1,"(%f, %f) (%f,%f)\n",creal(condinicial),cimag(condinicial),0.0,4.0);
    }
     if(maxiter!=0  ){
        if(modul(z1+1.7693)<tol){
            fprintf(fitxer2,"(%f, %f) (%f,%f)\n",creal(condinicial),cimag(condinicial),-1.7693,0.0);
        }
        if(modul(z1-(0.88465-0.58974*I))<tol){
            fprintf(fitxer3,"(%f, %f) (%f,%f)\n",creal(condinicial),cimag(condinicial),0.88465,-0.58974);
        }
        if(modul(z1-(0.88465+0.58974*I))<tol){
            fprintf(fitxer4,"(%f, %f) (%f,%f)\n",creal(condinicial),cimag(condinicial),0.88465,0.58974);
        }
        
     }
    return;
}
double complex p(double complex z){
    return (z+1.7693)*(z-(0.88465-0.58974*I))*(z-(0.88465+0.58974*I));
}
double complex derp(double complex z){
    return (z-(0.88465-0.58974*I))*(z-(0.88465+0.58974*I)) + (z+1.7693)*(z-(0.88465+0.58974*I)) + (z+1.7693)*(z-(0.88465-0.58974*I));
}
double modul(double complex z){
    double valor;
    valor= sqrt(creal(z)*creal(z) + cimag(z)*cimag(z));
    return valor;
}
