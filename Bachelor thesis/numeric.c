#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
void newton(double complex z0,double tol, FILE *fitxer);
double complex p(double complex z);
double complex derp(double complex z);
double modul(double complex z);

int main(void){
     int i;
     FILE *fitxer,*condin,*comprovacio;
     double tol=1e-5,h,m,r,v,pi;
     double complex z0,ci;
     pi=acos(-1);
    fitxer=fopen("fitxer.txt","w");
    if(fitxer==NULL){

        printf("problemes amb el fitxer de sortida\n");
        exit(1);
    }
    condin=fopen("condin.txt","w");
    if(condin==NULL){

        printf("problemes amb el fitxer de sortida\n");
        exit(1);
    }
    comprovacio=fopen("comprovacio.txt","w");
    if(comprovacio==NULL){

        printf("problemes amb el fitxer de sortida\n");
        exit(1);
    }
    r=2.28322;
    for(i=0;i<67;i++){
        v=(2*pi*i)/67;
        ci=cexp(-pi*I/2)*r*cexp(v*I);
        newton(ci,tol,comprovacio);
        fprintf(condin,"%f %f \n",creal(ci),cimag(ci));
    }
 
    m=4;
    h=-4;
    while(m>-4){
        h=-4;
        m=m-0.005;
        while(h<4){
            h=h+0.005;
            z0=h+m*I;
            newton(z0,tol,fitxer);
        }
    }
    
    fclose(fitxer);
    fclose(condin);
    fclose(comprovacio);
    return 0;
}
void newton(double complex z0,double tol, FILE *fitxer){
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
        fprintf(fitxer,"(%f, %f) (%f,%f) [0]\n",creal(condinicial),cimag(condinicial),0.0,4.0);
    }
     if(maxiter!=0  ){
        if(modul(z1+0.585063)<tol){
            fprintf(fitxer,"(%f, %f) (%f,%f) [1]\n",creal(condinicial),cimag(condinicial),-0.585063,0.0);
        }
        if(modul(z1-(-0.12853-0.93325*I))<tol){
            fprintf(fitxer,"(%f, %f) (%f,%f) [2]\n",creal(condinicial),cimag(condinicial),-0.12853,-0.93325);
        }
        if(modul(z1-(-0.12853+0.93325*I))<tol){
            fprintf(fitxer,"(%f, %f) (%f,%f) [3]\n",creal(condinicial),cimag(condinicial),-0.12853,0.93325);
        }
        if(modul(z1-(0.42106-0.49397*I))<tol){
            fprintf(fitxer,"(%f, %f) (%f,%f) [4]\n",creal(condinicial),cimag(condinicial),0.42106,-0.49397);
        }
        if(modul(z1-(0.42106+0.49397*I))<tol){
            fprintf(fitxer,"(%f, %f) (%f,%f) [5]\n",creal(condinicial),cimag(condinicial),0.42106,0.49397);
        }
        
     }
    return;
}
double complex p(double complex z){
    return 32*z*z*z*z*z +24*z*z*z+7;
}
double complex derp(double complex z){
    return 160*z*z*z*z +72*z*z;
}
double modul(double complex z){
    double valor;
    valor= sqrt(creal(z)*creal(z) + cimag(z)*cimag(z));
    return valor;
}
