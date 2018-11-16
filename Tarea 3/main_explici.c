//  Resolucion de la la ecuacion de calor de una barra 1-d dependiente del tiempo 
//
//
//  Created by Edgardo Dörner on 9/3/18.
//  Copyright © 2018 Edgardo Dörner. All rights reserved.
//  modified by Michael Haacke Concha and Vicente Salinas


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
const int nstep = 300;      // numero de pasos temporales
const double xi = 0.0;      // frontera inicial
const double xf  = L;       // frontera final
const int n = 21;//1 +  (int)(((xf - xi)/h)); // a partir de los parametros del sistema calculamos el numero de nodos de la grilla
const double dt = 0.1;//h*h/(2*kappa); // se ocupa esta definicion asi el codigo es estable

// Declaracion de las funciones a utilizar.
void F(double y[n], double dydt[n]);
void eRK4(double yn[n], double ynp[n], double dt);

// Funcion principal
int main(int argc, const char * argv[]) {

    printf("*************************************************\n");
    printf(" Resolucion de la ecuacion de calor  \n");
    printf(" mediante el metodo de Runge-Kutta de orden 4    \n");
    printf("*************************************************\n");
    
    // Discretizacion del eje espacial y temporal.
    double time_space[nstep][n];
    for (int i = 1; i < n-1; ++i)
    {
        time_space[0][i] = 0;
    }
    time_space[0][0] = T0;
    time_space[0][n-1] = TN;


    // Resolucion de la ecuacion diferencial mediante el Metodo de eRK-4.
    
    
    for (int i=0; i<nstep-1; i++) {        
        // Evolucion temporal mediante e-RK4.
        eRK4(time_space[i], time_space[i+1], dt);
    }
    // Guardamos los resultados en un archivo para su analisis posterior.
    FILE *fp = fopen("output.txt", "w");
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

// Funcion que contiene la ecuacion diferencial de problema.
// Esta consiste en la ecuacion de calor dependiente del tiempo bajo un sumidero de calor.
void F(double y[n], double dydt[n]) {
    for (int i = 1; i < n-1; ++i)
    {
    // dydt correponde ala derivada respecto al tiempo de la temperatura.
    // Al lado derecho de la diferencia se tiene la segunda derivada respecto a la distancia, utilizando la derivada central.
    // Al lado izquierdo, se tiene la funcion que representa el sumidero gaussiano centrado en L/2.
       dydt[i] = kappa*(y[i-1]-2*y[i]+y[i+1])/(h*h) - omega*exp(-pow(h*i - 0.5*L, 2.0)/pow(l, 2.0))/l;
    }
    dydt[0] = T0;   // se fija el borde izquierdo de la barra a una temperatura T0.
    dydt[n-1] = TN; // se fija el borde derecho de la barra a una temperatura TN.
    
    return;
}

// Funcion que contiene el metodo eRK-4. 
// Se aplica a este metodo a todos los valores del eje x al mismo tiempo.
void eRK4(double yn[n], double ynp[n], double dt) {
        double Y1[n],Y2[n],Y3[n],Y4[n];
        for (int i = 0; i < n; ++i){Y1[i] = yn[i];}
        F(Y1,Y2);
        for (int i = 0; i < n; ++i){Y2[i] = yn[i]+(dt/2)*Y2[i];}
        F(Y2,Y3);
        for (int i = 0; i < n; ++i){Y3[i] = yn[i]+(dt/2)*Y3[i];}
        F(Y3,Y4);
        for (int i = 0; i < n; ++i){Y4[i] = yn[i]+(dt)*Y4[i];ynp[i] = Y4[i];}
        
        double Y1S[n],Y2S[n],Y3S[n],Y4S[n];
        F(Y1,Y1S);
        F(Y2,Y2S);
        F(Y3,Y3S);
        F(Y4,Y4S);
        for (int i = 0; i < n; ++i)
            {
                ynp[i] = yn[i] + (dt/6)*(Y1S[i]+2*Y2S[i]+2*Y3S[i]+Y4S[i]);
            }
        ynp[0] = T0;
        ynp[n-1] = TN;
    
    return;
}
