//
//  main.c
//
//  Resolucion numerica del problema del pendulo doble mediante el metodo
//  Runge-Kutta explicito de orden cuatro (o eRK-4).
//
//  Creado por Michael Haacke Concha y Vicente Salinas Froemel el 11/10/18.
//
//  Requisitos: C o C++
//  Como correr: Compilar el programa (ej: g++ -o main.x main.c)
//               y luego correr en el terminal (ej: ./main.x)

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Parametros del problema

const double m = 1.0;         // masa de las masas del pendulo (Kg)
const double l = 1.0;         // longitud de los brazos del pendulo (m)
const double g = 9.8067;      // aceleracion de gravedad (m/s^2)
const double div_ml = 1.0/(m*pow(l,2));

const double phi10 = 0.0;     // angulo inicial masa 1
const double phi20 = 0.0;     // angulo inicial masa 2
const double p10 = 4.0;       // momento generalizado inicial masa 1
const double p20 = 2.0;       // momento generalizado inicial masa 2

void F(double y[4], double dydt[4]);
void eRK4(double yn[4], double ynp[4], double dt);

int main(int argc, const char * argv[]) {
    // Parametros del sistema.
    const double t0 = 0.0;              // instante inicial (s)
    const double tf  = 60.0;            // instante final (s)
    const double dt = 0.001;             // ancho temporal (s)
    
    // A partir de los parametros del sistema calculamos el numero de pasos
    // necesario.
    const int nsteps = 1 + (tf-t0)/dt;      // numero de pasos temporales
    
    printf("******************************************\n");
    printf(" Resolucion del pendulo doble mediante un \n");
    printf(" algoritmo eRK-4                          \n");
    printf("******************************************\n");
    
    // Discretizacion del eje temporal.
    double *t = (double*) malloc(nsteps*sizeof(double));
    for (int i=0; i<nsteps; i++) {
        t[i] = i*dt;
    }
    
    printf("Parametros del problema\n");
    printf("\t t0 \t\t= \t%10.2f (s)\n", t0);
    printf("\t tf \t\t= \t%10.2f (s)\n", t[nsteps-1]);
    printf("\t nsteps \t= \t%10d\n", nsteps);
    printf("\t dt \t\t= \t%10.2f (s)\n", dt);
    
    double *phi1 = (double*) malloc(sizeof(double)*nsteps);
    double *phi2 = (double*) malloc(sizeof(double)*nsteps);
    double *p1 = (double*) malloc(sizeof(double)*nsteps);
    double *p2 = (double*) malloc(sizeof(double)*nsteps);
    
    phi1[0] = phi10;
    phi2[0] = phi20;
    p1[0] = p10;
    p2[0] = p20;
    
    // Resolucion de la ecuacion diferencial mediante el Metodo de eRK-4.
    
    double yn[4] = {0.0, 0.0, 0.0, 0.0};
    double ynp[4] = {0.0, 0.0, 0.0, 0.0};
    
    for (int i=0; i<nsteps-1; i++) {
        // Variables en el instante actual.
        yn[0] = phi1[i];
        yn[1] = phi2[i];
        yn[2] = p1[i];
        yn[3] = p2[i];
        
        // Evolucion temporal mediante e-RK4.
        eRK4(yn, ynp, dt);
        
        // Guardamos los resultados para el instante siguiente.
        phi1[i+1] = ynp[0];
        phi2[i+1] = ynp[1];
        p1[i+1] = ynp[2];
        p2[i+1] = ynp[3];
    }
    
    // Guardamos los resultados en un archivo para su analisis posterior.
    FILE *fp = fopen("res_erk4.txt", "w");
    for (int i=0; i<nsteps; i++) {
        fprintf(fp, "%10.5lf, %10.5lf, %10.5lf, %10.5lf, %10.5lf\n", t[i], phi1[i], phi2[i], p1[i], p2[i]);
    }
    fclose(fp);

    // Guardamos la posicion de cada masa en coordenadas cartesianas.
    FILE *fp1 = fopen("res_erk4_car.txt", "w");
    for (int i=0; i<nsteps; i++) {
        fprintf(fp, "%10.5lf, %10.5lf, %10.5lf, %10.5lf, %10.5lf\n", t[i], l*sin(phi1[i]), 2*l-l*cos(phi1[i]),
        l*(sin(phi1[i])+sin(phi2[i])), 2*l-l*(cos(phi1[i])+cos(phi2[i])));
    }
    fclose(fp1);
    
    // Imprimimos los valores finales como muestra.
    printf("Valores finales eRK-4\n");
    printf("tf = %10.5f, phi1f = %10.5f, phi2f = %10.5f, p1f = %10.5f, p2f = %10.5f\n", t[nsteps-1], phi1[nsteps-1], phi2[nsteps-1], p1[nsteps-1], p2[nsteps-1]);
    
    // Limpieza
    free(t);
    free(phi1);
    free(phi2);
    free(p1);
    free(p2);
    
    return 0;
}

// Funcion con las ecuaciones diferenciales del problema del pendulo doble.
void F(double y[4], double dydt[4]) {
    dydt[0] = div_ml * (y[2]-y[3]*cos(y[0]-y[1]))/(1+pow(sin(y[0]-y[1]),2));
    dydt[1] = div_ml * (2*y[3]-y[2]*cos(y[0]-y[1]))/(1+pow(sin(y[0]-y[1]),2));
    dydt[2] = div_ml * (1/(1+pow(sin(y[0]-y[1]),2))) * (-y[2]*y[3]*sin(y[0]-y[1]) + (pow(y[2],2) + 
              2*pow(y[3],2) - 2*y[2]*y[3]*cos(y[0]-y[1]))/(1+pow(sin(y[0]-y[1]),2)) * sin(y[0]-y[1])*cos(y[0]-y[1])) - 
              2*m*g*l*sin(y[0]);
    dydt[3] = div_ml * (1/(1+pow(sin(y[0]-y[1]),2))) * (y[2]*y[3]*sin(y[0]-y[1]) - (pow(y[2],2) + 
              2*pow(y[3],2) - 2*y[2]*y[3]*cos(y[0]-y[1]))/(1+pow(sin(y[0]-y[1]),2)) * sin(y[0]-y[1])*cos(y[0]-y[1])) - 
              m*g*l*sin(y[1]);
    return;
}

// Funcion con el metodo de Runge-Kutta de orden cuatro.
void eRK4(double yn[4], double ynp[4], double dt) {
    double Y1[4],Y2[4],Y3[4],Y4[4];
    Y1[0] = yn[0];Y1[1] = yn[1];Y1[2] = yn[2];Y1[3] = yn[3];
    F(Y1,Y2);
    Y2[0] = yn[0]+(dt/2)*Y2[0]; Y2[1] = yn[1]+(dt/2)*Y2[1]; Y2[2] = yn[2]+(dt/2)*Y2[2]; Y2[3] = yn[3]+(dt/2)*Y2[3];
    F(Y2,Y3);
    Y3[0] = yn[0]+(dt/2)*Y3[0]; Y3[1] = yn[1]+(dt/2)*Y3[1]; Y3[2] = yn[2]+(dt/2)*Y3[2]; Y3[3] = yn[3]+(dt/2)*Y3[3];
    F(Y3,Y4);
    Y4[0] = yn[0]+(dt)*Y4[0]; Y4[1] = yn[1]+(dt)*Y4[1]; Y4[2] = yn[2]+(dt)*Y4[2]; Y4[3] = yn[3]+(dt)*Y4[3];
    ynp[0] = Y4[0];ynp[1] = Y4[1];ynp[2] = Y4[2];ynp[3] = Y4[3];
    double Y1S[4],Y2S[4],Y3S[4],Y4S[4];
    F(Y1,Y1S);
    F(Y2,Y2S);
    F(Y3,Y3S);
    F(Y4,Y4S);
    ynp[0] = yn[0] + (dt/6)*(Y1S[0]+2*Y2S[0]+2*Y3S[0]+Y4S[0]);
    ynp[1] = yn[1] + (dt/6)*(Y1S[1]+2*Y2S[1]+2*Y3S[1]+Y4S[1]);
    ynp[2] = yn[2] + (dt/6)*(Y1S[2]+2*Y2S[2]+2*Y3S[2]+Y4S[2]);
    ynp[3] = yn[3] + (dt/6)*(Y1S[3]+2*Y2S[3]+2*Y3S[3]+Y4S[3]);
    
    return;
}

