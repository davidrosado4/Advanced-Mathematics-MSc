#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<complex.h>
/*Metode de Newton per p(z)=z^3-1*/
void newton(double complex z0,double tol, FILE *arrel1, FILE *arrel2,FILE *arrel3);
double complex p(double complex z);
double complex derp(double complex z);
double modul(double complex z);

int main(void){
    double tol=1e-5,h,m;
    FILE *arrel1,*arrel2,*arrel3;
    double complex z0;
    
    /*obro tres fitxers de sortida, un per a cada arrel*/
    arrel1=fopen("archiu1.txt","w");
    if(arrel1==NULL){

        printf("problemes amb el fitxer de sortida\n");
        exit(1);
    }
    arrel2=fopen("archiu2.txt","w");
    if(arrel2==NULL){

        printf("problemes amb el fitxer de sortida\n");
        exit(1);
    }
    arrel3=fopen("archiu3.txt","w");
    if(arrel3==NULL){

        printf("problemes amb el fitxer de sortida\n");
        exit(1);
    }
   
    /*xarxa equiespaiada de punts per fer el newton en el quadrat [-5,5]^2*/
    m=2;
    h=-2;
    while(m>-2){
        h=-2;
        m=m-0.003;
        while(h<2){
            h=h+0.003;
            z0=h+m*I;
            newton(z0,tol,arrel1,arrel2,arrel3);
        }
    }
    fclose(arrel1);
    fclose(arrel2);
    fclose(arrel3);
    return 0;
}

void newton(double complex z0,double tol, FILE *arrel1, FILE *arrel2,FILE *arrel3){
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
    /*si el metode ha convergit a alguna de les tres arrels, passo la condinicial inicial al fitxer*/
    /*Despres dibuixarem el tres fitxers*/
    if(maxiter!=0 && modul(derp(z1))>tol){
        if(modul(z1-1)<tol){
            fprintf(arrel1,"%f %f\n",creal(condinicial),cimag(condinicial));
        }
        if(modul(z1-(-0.5+(sqrt(3)/2)*I))<tol){
            fprintf(arrel2,"%f %f\n",creal(condinicial),cimag(condinicial));
        }
        if(modul(z1-(-0.5-(sqrt(3)/2)*I))<tol){
            fprintf(arrel3,"%f %f\n",creal(condinicial),cimag(condinicial));
        }
    
    }
    return;
}

double complex p(double complex z){
    return (z-1)*(z-(-0.5+(sqrt(3)/2)*I))*(z-(-0.5-(sqrt(3)/2)*I));
}
double complex derp(double complex z){
    return (z-(-0.5+(sqrt(3)/2)*I))*(z-(-0.5-(sqrt(3)/2)*I)) + (z-1)*(z-(-0.5-(sqrt(3)/2)*I)) + (z-1)*(z-(-0.5+(sqrt(3)/2)*I));
}
double modul(double complex z){
    double valor;
    valor= sqrt(creal(z)*creal(z) + cimag(z)*cimag(z));
    return valor;
}
