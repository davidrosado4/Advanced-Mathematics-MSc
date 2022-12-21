#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
void newton1(double complex z0,double tol, FILE *fitxer0,FILE *fitxer1,FILE *fitxer2,FILE *fitxer3,FILE *fitxer4,FILE *fitxer5);
double complex p(double complex z);
double complex derp(double complex z);
double modul(double complex z);
void newton2(double complex z0,double tol, FILE *comprovacio);

int main(void){
     double i;
     FILE *fitxer0,*fitxer1,*fitxer2,*fitxer3,*fitxer4,*fitxer5,*condin,*comprovacio;
     double tol=1e-5,h,m,r;
     double complex z0,ci;
    fitxer0=fopen("fitxer0.txt","w");
    if(fitxer0==NULL){

        printf("problemes amb el fitxer de sortida\n");
        exit(1);
    }
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
    fitxer5=fopen("fitxer5.txt","w");
    if(fitxer5==NULL){

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
    for(i=0;i<=18;i=i+0.75){
        ci=r+(-9+i)*I;
        newton2(ci,tol,comprovacio);
        fprintf(condin,"%f %f \n",creal(ci),cimag(ci));
    }
    m=10;
    h=-10;
    while(m>-10){
        h=-10;
        m=m-0.01;
        while(h<10){
            h=h+0.01;
            z0=h+m*I;
            newton1(z0,tol,fitxer0,fitxer1,fitxer2,fitxer3,fitxer4,fitxer5);
        }
    }
    
    fclose(fitxer0);
    fclose(fitxer1);
    fclose(fitxer2);
    fclose(fitxer3);
    fclose(fitxer4);
    fclose(fitxer5);
    fclose(condin);
    fclose(comprovacio);
    return 0;
}
void newton1(double complex z0,double tol, FILE *fitxer0,FILE *fitxer1,FILE *fitxer2,FILE *fitxer3,FILE *fitxer4,FILE *fitxer5){
    int maxiter;
    double diferencia;
    double complex condinicial,z1;
    z1=z0;
    condinicial=z0;
    maxiter=50;
    diferencia=1;
    
    
    while(diferencia >tol && maxiter!=0 && modul(derp(z1)+p(z1))>tol && modul(p(z1))>tol){
        z0=z1;
        z1=z0-p(z0)/(derp(z0)+p(z0));
        maxiter--;
        diferencia=modul(z1);
    }
    if(maxiter==0){
        fprintf(fitxer0,"(%f, %f) (%f,%f) [0]\n",creal(condinicial),cimag(condinicial),0.0,4.0);
    }
    if(maxiter!=0  ){
        if(modul(z1+0.585063)<tol){
            fprintf(fitxer1,"(%f, %f) (%f,%f) [1]\n",creal(condinicial),cimag(condinicial),-0.585063,0.0);
        }
        if(modul(z1-(-0.12853-0.93325*I))<tol){
            fprintf(fitxer2,"(%f, %f) (%f,%f) [2]\n",creal(condinicial),cimag(condinicial),-0.12853,-0.93325);
        }
        if(modul(z1-(-0.12853+0.93325*I))<tol){
            fprintf(fitxer3,"(%f, %f) (%f,%f) [3]\n",creal(condinicial),cimag(condinicial),-0.12853,0.93325);
        }
        if(modul(z1-(0.42106-0.49397*I))<tol){
            fprintf(fitxer4,"(%f, %f) (%f,%f) [4]\n",creal(condinicial),cimag(condinicial),0.42106,-0.49397);
        }
        if(modul(z1-(0.42106+0.49397*I))<tol){
            fprintf(fitxer5,"(%f, %f) (%f,%f) [5]\n",creal(condinicial),cimag(condinicial),0.42106,0.49397);
        }
        
     }
    return;
}
void newton2(double complex z0,double tol, FILE *comprovacio){
     int maxiter;
    double diferencia;
    double complex condinicial,z1;
    z1=z0;
    condinicial=z0;
    maxiter=50;
    diferencia=1;
    
    
    while(diferencia >tol && maxiter!=0 && modul(derp(z1)+p(z1))>tol && modul(p(z1))>tol){
        z0=z1;
        z1=z0-p(z0)/(derp(z0)+p(z0));
        maxiter--;
        diferencia=modul(z1);
    }
    if(maxiter==0){
        fprintf(comprovacio,"(%f, %f) (%f,%f) [0]\n",creal(condinicial),cimag(condinicial),0.0,4.0);
    }
     if(maxiter!=0  ){
        if(modul(z1+0.585063)<tol){
            fprintf(comprovacio,"(%f, %f) (%f,%f) [1]\n",creal(condinicial),cimag(condinicial),-0.585063,0.0);
        }
        if(modul(z1-(-0.12853-0.93325*I))<tol){
            fprintf(comprovacio,"(%f, %f) (%f,%f) [2]\n",creal(condinicial),cimag(condinicial),-0.12853,-0.93325);
        }
        if(modul(z1-(-0.12853+0.93325*I))<tol){
            fprintf(comprovacio,"(%f, %f) (%f,%f) [3]\n",creal(condinicial),cimag(condinicial),-0.12853,0.93325);
        }
        if(modul(z1-(0.42106-0.49397*I))<tol){
            fprintf(comprovacio,"(%f, %f) (%f,%f) [4]\n",creal(condinicial),cimag(condinicial),0.42106,-0.49397);
        }
        if(modul(z1-(0.42106+0.49397*I))<tol){
            fprintf(comprovacio,"(%f, %f) (%f,%f) [5]\n",creal(condinicial),cimag(condinicial),0.42106,0.49397);
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
