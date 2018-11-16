//  main_implicit.c
//  
//
//  Resolucion de la la ecuacion de calor de una barra 1-d dependiente del tiempo
//  mediante el metodo de crank-michelson
//
//  Created by Edgardo Dörner on 9/3/18.
//  Copyright © 2018 Edgardo Dörner. All rights reserved.
//  Modified by Michael Haacke and Vicente Salinas



#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Parametros del problema 
const double L = 10.0;  // largo de la barra
const double T0 = 0.0;  // temperatura extremo izquierdo
const double TN = 2.0;  // temperatura extremo derecho
double omega = -0.4;    // altura/profundidad de la fuente/sumidero
double l = 1.0;         // ancho de la fuente/sumidero
double kappa = 1.0;      // difusividad termica
const double h = 0.5;       // ancho espacial
const int nstep = 300;
const double xi = 0.0;      // frontera inicial
const double xf  = L;       // frontera final
const int n = 21;   //1 +  (int)(((xf - xi)/h));
double dt;

// Declaracion de las funciones a utilizar.
void F(double y[n], double dydt[n]);
void tridag(double r[], double u[], int ne);

// Funcion principal
int main(int argc, const char * argv[]) {
   // Parametros del sistema.
    
    dt = 0.1;

    printf("*************************************************\n");
    printf(" Resolucion de la ecuacion de calor  \n");
    printf(" mediante diferencias finitas                    \n");
    printf("*************************************************\n");
    
    // Discretizacion del eje espacial y temporal.
    double time_space[nstep][n];
    for (int i = 1; i < n-1; ++i)
    {
        time_space[0][i] = 0;
    }
    time_space[0][0] = T0;
    time_space[0][n-1] = TN;

    
    // Resolucion de la ecuacion diferencial mediante el metodo de crank-michelson, con este metodo implicito
    // se logra crear un sistema de ecuaciones tridiagonales, las cuales son resuletas con un algoritmo obtenido
    // del libro "Numerical recipes in C.” Cambridge University Press 1 (1988), con esto se obtiene los resultados
    // deseados forzando siempre las condiciones de borde
    double rest[n] = {0};

    for (int i=0; i<nstep-1; i++) {     
        time_space[i][0] = T0;
        time_space[i][n-1] = TN;
        F(time_space[i],rest);
        tridag(rest,time_space[i+1],n);
        time_space[i][0] =T0;
        time_space[i][n-1] =TN;
           
    }
    // Guardamos los resultados en un archivo para su analisis posterior, las columnas represetan los tiempos di-
    // derentes y la primera columna es el eje x
    FILE *fp = fopen("output1.txt", "w");
    for (int i=0; i<n; i++) {
        fprintf(fp, "%10.5lf, ",i*h);
        for (int j = 0; j < nstep; ++j)
        {

           fprintf(fp, "%10.5lf, ",time_space[j][i]);
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
    // Imprimimos los valores finales como muestra.
    
    return 0;
}


// funcion que obtiene el vector igualado , es la derivada parcial mas la solucion particular de la edo.
void F(double y[n], double dydt[n]) {

    dydt[0] = T0;
    dydt[n-1] = TN;
    for (int i = 1; i < n-1; ++i)
    {
       dydt[i] = y[i] + kappa*dt/(2*h*h)*(y[i-1]-2*y[i]+y[i+1]) - dt*omega*exp(-pow(h*i - 0.5*L, 2.0)/pow(l, 2.0))/l;
    }
    
    return;
}

// from Numerica recipes in C and modified by michael (a,b,c) are all costants after expanding with crank-michelson
void tridag(double r[], double u[], int ne)
//Solves for a vector u[1..n] the tridiagonal linear set given by equation (2.4.1). a[1..n], b[1..n], c[1..n], and r[1..n] are input vectors and are not modified.
{
double a = -kappa*dt/(2*h*h);
double b = 1+kappa*dt/(h*h);
double c = a;
unsigned long j;
double bet,gam[ne];
u[1]=r[1]/(bet=b);
for (j=2;j<=ne;j++) 
{ //Decomposition and forward substitution.
    gam[j]=c/bet;
    bet=b-a*gam[j];
    if (bet == 0.0) printf("ERROR bet = 0\n");
    u[j]=(r[j]-a*u[j-1])/bet;
    //Algorithm fails; see be- low.
}
for (j=(ne-1);j>=1;j--) u[j] -= gam[j+1]*u[j+1];
}
//Backsubstitution.










