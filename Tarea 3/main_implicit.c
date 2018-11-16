//  main.c
//  p01-explicit_euler
//
//  Resolucion de la ley de Stefan-Boltzmann mediante diferentes
//  integradores.
//
//  Created by Edgardo Dörner on 9/3/18.
//  Copyright © 2018 Edgardo Dörner. All rights reserved.
//  Modified by Michael Haacke

// Parametros del problema

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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

void F(double y[n], double dydt[n]);
void tridag(double r[], double u[], int ne);

int main(int argc, const char * argv[]) {
   // Parametros del sistema.
    
    dt = 0.1;
    
    // A partir de los parametros del sistema calculamos el numero de
    // nodos de la grilla
    

    
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

   
    
    printf("Parametros del problema\n");
    printf("\t a \t\t= \t%10.2f \n", time_space[0][0]);
    printf("\t b \t\t= \t%10.2f \n", time_space[nstep-1][n-1]);
    printf("\t n \t\t= \t%10d\n", n);
    printf("\t h \t\t= \t%10.3f \n", h);
    
    // Resolucion de la ecuacion diferencial mediante el Metodo de eRK-4.
    
    double rest[n] = {0};

    for (int i=0; i<nstep-1; i++) {     
        time_space[i][0] = T0;
        time_space[i][n-1] = TN;
        F(time_space[i],rest);   
        // Evolucion temporal mediante e-RK4.
        tridag(rest,time_space[i+1],n);
        time_space[i][0] =T0;
        time_space[i][n-1] =TN;
           
    }
    // Guardamos los resultados en un archivo para su analisis posterior.
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

void F(double y[n], double dydt[n]) {

    dydt[0] = T0;
    dydt[n-1] = TN;
    for (int i = 1; i < n-1; ++i)
    {
       dydt[i] = y[i] + kappa*dt/(2*h*h)*(y[i-1]-2*y[i]+y[i+1]) - dt*omega*exp(-pow(h*i - 0.5*L, 2.0)/pow(l, 2.0))/l;
    }
    
    return;
}

// from Numerica recipes in C and modified vy michael (a,b,c) are all costants after expanding with crank-michelson
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










