#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


//-------------------------------------------------------------------------------------------------------
//									Declaracion de Funciones
//-------------------------------------------------------------------------------------------------------
//									Generador random

struct rng {
    // Parametros del pseudo generador de numeros aleatorios y numero pseudo-aleatorio 
    long int a ,c ,m, x;
};
typedef struct rng Rng;

struct particula {
    // Parametros Usados en la reconstruccion de la particula
    double x,y,z,u,v,w;
    int region,nueva_region, escapada, interactua;
};
typedef struct particula Particula;

struct esfera {
    // Parametros Usados para la geometria de la esfera
    double distancias[N_REGION], sigma_t[N_REGION], sigma_a[N_REGION];
    // parametros usado para la estadistica
    double particulas_depositadas[N_REGION + 1], particulas_depositadas_pow2[N_REGION + 1];
};
typedef struct esfera Esfera;


void initRng(long int seed, long int a, long int c, long int m,
             Rng *rng);
long int setRng(Rng *rng);

double getRng(Rng *rng);

//									Resto de las funciones 

void dispercion_isotropica(Rng *rng,double u, double v, double w,double output[3]);

double director_por_pos(Particula *particula);

double pos_pow_2(Particula *particula);




