#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
void newton1(double complex z0,double tol, FILE *fitxer0,FILE *fitxer1,FILE *fitxer2,FILE *fitxer3,FILE *fitxer4,FILE *fitxer5,FILE *fitxer6,FILE *fitxer7,FILE *fitxer8);
double complex p(double complex z);
double complex derp(double complex z);
double modul(double complex z);
void newton2(double complex z0,double tol, FILE *comprovacio);

int main(void){
     int i;
     FILE *fitxer0,*fitxer1,*fitxer2,*fitxer3,*fitxer4,*fitxer5,*fitxer6,*fitxer7,*fitxer8,*condin,*comprovacio;
     double tol=1e-5,h,m,r,v,pi;
     double complex z0,ci;
     pi=acos(-1);
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
    fitxer6=fopen("fitxer6.txt","w");
    if(fitxer6==NULL){

        printf("problemes amb el fitxer de sortida\n");
        exit(1);
    }
    fitxer7=fopen("fitxer7.txt","w");
    if(fitxer7==NULL){

        printf("problemes amb el fitxer de sortida\n");
        exit(1);
    }
    fitxer8=fopen("fitxer8.txt","w");
    if(fitxer8==NULL){

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
    r=2.33495;
    for(i=0;i<139;i++){
        v=(2*pi*i)/139;
        ci=cexp(-pi*I/2)*r*cexp(v*I);
        newton2(ci,tol,comprovacio);
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
            newton1(z0,tol,fitxer0,fitxer1,fitxer2,fitxer3,fitxer4,fitxer5,fitxer6,fitxer7,fitxer8);
        }
    }
    
    fclose(fitxer0);
    fclose(fitxer1);
    fclose(fitxer2);
    fclose(fitxer3);
    fclose(fitxer4);
    fclose(fitxer5);
    fclose(fitxer6);
    fclose(fitxer7);
    fclose(fitxer8);
    fclose(condin);
    fclose(comprovacio);
    return 0;
}
void newton1(double complex z0,double tol, FILE *fitxer0,FILE *fitxer1,FILE *fitxer2,FILE *fitxer3,FILE *fitxer4,FILE *fitxer5,FILE *fitxer6,FILE *fitxer7,FILE *fitxer8){
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
        fprintf(fitxer0,"(%f, %f) (%f,%f) [0]\n",creal(condinicial),cimag(condinicial),0.0,4.0);
    }
     if(maxiter!=0  ){
        if(modul(z1-(-0.44981-0.19030*I))<tol){
            fprintf(fitxer1,"(%f, %f) (%f,%f) [1]\n",creal(condinicial),cimag(condinicial),-0.44981,-0.19030);
        }
        if(modul(z1-(-0.39964+0.23107*I))<tol){
            fprintf(fitxer2,"(%f, %f) (%f,%f) [2]\n",creal(condinicial),cimag(condinicial),-0.39964,0.23107);
        }
        if(modul(z1-(-0.05373-0.52765*I))<tol){
            fprintf(fitxer3,"(%f, %f) (%f,%f) [3]\n",creal(condinicial),cimag(condinicial),-0.05373,-0.52765);
        }
        if(modul(z1-(-0.03414+0.43680*I))<tol){
            fprintf(fitxer4,"(%f, %f) (%f,%f) [4]\n",creal(condinicial),cimag(condinicial),-0.03414,0.43680);
        }
        if(modul(z1-(0.08187+0.87686*I))<tol){
            fprintf(fitxer5,"(%f, %f) (%f,%f) [5]\n",creal(condinicial),cimag(condinicial),0.08187,0.87686);
        }
        if(modul(z1-(0.09381-0.86392*I))<tol){
            fprintf(fitxer6,"(%f, %f) (%f,%f) [6]\n",creal(condinicial),cimag(condinicial),0.09381,-0.86392);
        }
        if(modul(z1-(0.35470+0.24024*I))<tol){
            fprintf(fitxer7,"(%f, %f) (%f,%f) [7]\n",creal(condinicial),cimag(condinicial),0.35470,0.24024);
        }
         if(modul(z1-(0.40694-0.20310*I))<tol){
            fprintf(fitxer8,"(%f, %f) (%f,%f) [8]\n",creal(condinicial),cimag(condinicial),0.40694,-0.20310);
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
    
    
    while(diferencia >tol && maxiter!=0 && modul(derp(z1))>tol && modul(p(z1))>tol){
        z0=z1;
        z1=z0-p(z0)/derp(z0);
        maxiter--;
        diferencia=modul(z1);
    }
    if(maxiter==0){
        fprintf(comprovacio,"(%f, %f) (%f,%f) [0]\n",creal(condinicial),cimag(condinicial),0.0,4.0);
    }
     if(maxiter!=0  ){
        if(modul(z1-(-0.44981-0.19030*I))<tol){
            fprintf(comprovacio,"(%f, %f) (%f,%f) [1]\n",creal(condinicial),cimag(condinicial),-0.44981,-0.19030);
        }
        if(modul(z1-(-0.39964+0.23107*I))<tol){
            fprintf(comprovacio,"(%f, %f) (%f,%f) [2]\n",creal(condinicial),cimag(condinicial),-0.39964,0.23107);
        }
        if(modul(z1-(-0.05373-0.52765*I))<tol){
            fprintf(comprovacio,"(%f, %f) (%f,%f) [3]\n",creal(condinicial),cimag(condinicial),-0.05373,-0.52765);
        }
        if(modul(z1-(-0.03414+0.43680*I))<tol){
            fprintf(comprovacio,"(%f, %f) (%f,%f) [4]\n",creal(condinicial),cimag(condinicial),-0.03414,0.43680);
        }
        if(modul(z1-(0.08187+0.87686*I))<tol){
            fprintf(comprovacio,"(%f, %f) (%f,%f) [5]\n",creal(condinicial),cimag(condinicial),0.08187,0.87686);
        }
        if(modul(z1-(0.09381-0.86392*I))<tol){
            fprintf(comprovacio,"(%f, %f) (%f,%f) [6]\n",creal(condinicial),cimag(condinicial),0.09381,-0.86392);
        }
        if(modul(z1-(0.35470+0.24024*I))<tol){
            fprintf(comprovacio,"(%f, %f) (%f,%f) [7]\n",creal(condinicial),cimag(condinicial),0.35470,0.24024);
        }
         if(modul(z1-(0.40694-0.20310*I))<tol){
            fprintf(comprovacio,"(%f, %f) (%f,%f) [8]\n",creal(condinicial),cimag(condinicial),0.40694,-0.20310);
        }
        
     }
    return;
}
double complex p(double complex z){
    return 256*z*z*z*z*z*z*z*z + 192*z*z*z*z*z*z + 32*z*z*z*z*z + I*2*z + 2;
}
double complex derp(double complex z){
    return 2048*z*z*z*z*z*z*z + 1152*z*z*z*z*z + 160*z*z*z*z + I*2;
}
double modul(double complex z){
    double valor;
    valor= sqrt(creal(z)*creal(z) + cimag(z)*cimag(z));
    return valor;
}
