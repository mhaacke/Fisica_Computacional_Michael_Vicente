//  main.c 
//
//  Modelo transporte de particulas sobre un nucleo fisionable de forma esferica con varias tipo de 
//	materiales mediante el metodo Monte Carlo. En especifico determina el porcentage de particulas 
//	transmitida y absorbidas si estas vienen del nucleo de la esfera y escribe en un archivo todas
//	las historia, se calcula el tiempo de vuelo como la exponecial de -∑_t y luego que tipo de 
//	interraccion es la que sucede, Se ocupa un generador de numeros pseudo-random del estilo 
//	"Generadores lineales congruenciales" con valores obtenidos en Numerical Recipes in C.
//
//	Parte del codigo es obtenido de la clase de fisica computacional del profesor Edgardo Dörner 

#define N_REGION 4				// Numero de regiones es una constante al tiempo de compilacion ¡TIENE QUE ESTAR AQUI!

#include "main.h"

//-------------------------------------------------------------------------------------------------------
//									Seteo de parametros del problema
//-------------------------------------------------------------------------------------------------------

const double d[] = {0.5, 0.1, 0.1,0.1};    // espesor de las regiones del atenuador (cm)
const double sigma_t[] = {0.05, 0.1, 10.0,100}; // seccion eficaz total (cm^-1)
const double sigma_a[] = {0.005,0.01, 0.1, 10.0};//seccion eficaz de absorcion (cm^-1)
const int nhist = 10000000; // numero de historias
const int nbatch = 10;  // grupos estadisticos
const int ptracks = 1;  // guardar posicion particulas
const double velocidad_inicial = 1;

//-------------------------------------------------------------------------------------------------------
//									Definicion del generador pseudo-Random
//-------------------------------------------------------------------------------------------------------


void initRng(long int seed, long int a, long int c, long int m,Rng *rng) {
    /* Inicializacion del generador */
    rng->a = a; rng->c = c; rng->m = m; rng->x = seed;

    printf("Generador pseudo-aleatorio LCG con parametros:\n");
    printf("\t a = %ld\n", rng->a);
    printf("\t c = %ld\n", rng->c);
    printf("\t m = %ld\n", rng->m);
    printf("\t seed = %ld\n", rng->x);
    return;
}

double getRng(Rng *rng)
{
	rng->x = ((rng->a*rng->x + rng->c) % rng->m);
	return (double)rng->x / rng->m; 
}

//-------------------------------------------------------------------------------------------------------
//									Funcion auxiliares
//-------------------------------------------------------------------------------------------------------

void dispercion_isotropica(Rng *rng,double u, double v, double w,double output[3]) 	
{																				
	// esta funcion toma el vector director  anterior y retorna un nuevo vector director output que es una 
	// dispersion isotropica del vector director original, codigo obtenido de lo visto en clases
	double rnno = getRng(rng);
	double theta0 = acos(2.0*rnno - 1.0);
	rnno = getRng(rng);
	double phi0 = 2.0*M_PI*rnno;
	// Calculamos la nueva direccion de la particula. Si el cambio del angulo polar es muy pequeño, w ~ 1 y s ~ 0,
	// por lo tanto agregamos un caso especial para manejar esa situacion.
	if (pow(u, 2.0) + pow(v, 2.0) < 1.0E-20)
	{
		output[0] = sin(theta0)*cos(phi0);
		output[1] = sin(theta0)*sin(phi0);
		output[2] = w*cos(theta0);
	}
	else 
	{
		double s = sqrt(1 - pow(w, 2.0));
		output[0] = (w*u/s)*sin(theta0)*cos(phi0) - (v/s)*sin(theta0)*sin(phi0) + u*cos(theta0);
		output[1] = (w*v/s)*sin(theta0)*cos(phi0) + (u/s)*sin(theta0)*sin(phi0) + v*cos(theta0);
		output[2] = -s*(sin(theta0)*cos(phi0)) + w*cos(theta0);
	}
}

double director_por_pos(Particula *particula)
{// esta funcion retorna el producto punto del vector director con la poscion
	return particula->u*particula->x + particula->v*particula->y + particula->w*particula->z;
}


double pos_pow_2(Particula *particula)
{
	return particula->x*particula->x + particula->y*particula->y + particula->z*particula->z;
}


//-------------------------------------------------------------------------------------------------------
//									Funcion Main
//-------------------------------------------------------------------------------------------------------


int main(int argc, const char * argv[]) 
{
	//------------Definicion y inicilizacion del generador pseudo-random---------------------------------

	Rng rng1;
	initRng(0, 1664525, 1013904223, 4294967296, &rng1);
	
	//------------Definicion y inicilizacion de la struct de esfera y particula--------------------------

	Esfera *esfera = (Esfera *)malloc(sizeof(Esfera));
	Particula *particula = (Particula *)malloc(sizeof(Particula));
	
	// Setea los primos valores a la estructura
	esfera->distancias[0] = d[0];
	esfera->sigma_a[0] = sigma_a[0];
	esfera->sigma_t[0] = sigma_t[0];
	for (int i = 1; i < N_REGION; ++i) // loop que setea los array a la estrucura esfera
	{
		esfera->distancias[i] = d[i] + esfera->distancias[i-1];
		esfera->sigma_a[i] = sigma_a[i];
		esfera->sigma_t[i] = sigma_t[i];
	}
	for (int i=0; i<N_REGION+1; i++) 
		{
			esfera->particulas_depositadas[i] = 0;
			esfera->particulas_depositadas_pow2[i] = 0;
		}
	*particula = (Particula){0}; // Setea la estructura a 0

	//------------Inicializacion del archivo-------------------------------------------------------------
	// Inicializamos el archivo donde guardaremos las trayectorias de cada particula
    FILE *fp;
    fp = fopen("ptracks.lst","w");
    // Guardamos el numero de historias que seran registradas.
    fprintf(fp, "%d\n", nhist);

    //------------Declaracion de variables usadas en el loop---------------------------------------------
    double particulas_depositadas[N_REGION + 1] = {0};
    double new_director[3] = {0};
    double random_number,camino_libre;
    double x_temp,y_temp,z_temp,factor;


    //------------Loop principal-------------------------------------------------------------------------
	long nperbatch = nhist/nbatch;
    
    for (long ibatch=0; ibatch<nbatch; ibatch++)
	{     
		for (long ihist=0; ihist<nperbatch; ihist++)
        {
        	//-------- Marcamos el inicio de cada historia-----------------------------------------------
            fprintf(fp, "%s\n", "&");

			//--------Inicializacion de la particula-----------------------------------------------------
            particula->x = 0;particula->y = 0;particula->z = 0;
            dispercion_isotropica(&rng1,0,0,1,new_director);// se calcula una dispersion isotropica con w
            particula->u = new_director[0];					// 1 y se setea como los valores iniciales, notamos 
            particula->v = new_director[1];					// que el vector direcort es unitario
            particula->w = new_director[2];
            particula->region = 0;particula->nueva_region = 0;
            particula->escapada = 0;particula->interactua = 1;

            // Comienza el recorrido de la particula, esta parte es muy similar a la clase 27
            while (1) 
            {
				//-------------------------TRANSPORTE DE LA PARTICULA----------------------------------------------------------------
            	do
            	{
	            	//obtenemos el camino libre de la exponecial negativa
	            	random_number = getRng(&rng1);
	            	camino_libre = -(1.0/esfera->distancias[particula->region])*log(1.0 - random_number);
	            	// guardamos la region actual como la nueva region
	            	particula->nueva_region = particula->region;

	            	// Calculamos si con el nuevo transporte cambia de region, de ser asi calculamos en que parte del
	            	// cascaron sucede esto, calculamos primero si la particula solo sigue su camino
					x_temp = particula->x + particula->u * camino_libre;
					y_temp = particula->y + particula->v * camino_libre;
					z_temp = particula->z + particula->w * camino_libre;

					if (sqrt(x_temp*x_temp+y_temp*y_temp+z_temp*z_temp) > esfera->distancias[particula->region])
					{
						particula->nueva_region +=1; // marcamos que esta una region mas afuera y chekeamos si se salio de la geometria
						if (particula->nueva_region > N_REGION-1){ particula->escapada = 1;}// La particula abandona la geometria.

						// el factor viene de el punto de interseccion de la recta que genera el camino nuevo con el cascaron esferico es
						// x_r = vector director    x_i = posicion inicial     punto = x_r*factor + x_i
						// factor = (x_r*x_i + sqrt((x_r*x_i)^2 - x_i^2 +r^2)  , esto viene de sacar la ecuacion de la recta, sale una ecuacion
						// cuadratica y se exije que el factor sea positvo asi apunta a la misma direccion que el vector direcor
						factor = -director_por_pos(particula) + sqrt(director_por_pos(particula)*director_por_pos(particula) - pos_pow_2(particula) + esfera->distancias[particula->region]*esfera->distancias[particula->region]);
						x_temp = particula->u*factor + particula->x;
						y_temp = particula->v*factor + particula->y;
						z_temp = particula->w*factor + particula->z;
						//printf("radio observado: %lf limite del radio: %lf\n", sqrt(x_temp*x_temp+y_temp*y_temp+z_temp*z_temp) , esfera->distancias[particula->region]);
					}
					if (particula->escapada == 1) {
                        // La particula abandona la geometria.
                        break;
                    }
					// Guardamos la posicion de la particula antes del transporte.
	                fprintf(fp, "%10.5f %10.5f %10.5f ", particula->x, particula->y, particula->z);
					
					// Transportamos la particula
					particula->x = x_temp; particula->y = y_temp; particula->z = z_temp;

					// Guardamos la posicion de la particula despues del transporte.
					fprintf(fp, "%10.5f %10.5f %10.5f\n",  particula->x, particula->y, particula->z);

					// si no se cambia de region la particula debe interactuar y si cambio volvemos los cambios permanentes
					if (particula->region == particula->nueva_region){particula->interactua = 0;}
					else {particula->region = particula->nueva_region;}// Actualizamos la region de la particula y continuamos el transporte
				}while(particula->interactua);
            	//-------------------------FIN DEL TRANSPORTE DE LA PARTICULA-------------------------------------------------------

				if (particula->escapada == 1)
				{
					// La particula fue descartada. La depositamos y finalizamos
                    // su rastreo.
                    particula->region = particula->nueva_region;
                    particulas_depositadas[particula->region] += 1;
                    break;
				}


            	//-------------------------INTERACCION DE LA PARTICULA-------------------------------------------------------
				// Debemos seleccionar una interaccion
                random_number = getRng(&rng1);
                if (random_number <= sigma_a[particula->region]/sigma_t[particula->region]) 
                {
                    // La particula es absorbida
                    particulas_depositadas[particula->region] += 1;
                    break;
                } else 
                {
                	// La particula es dispersada, debemos muestrear el angulo
                    // de dispersion. Muestreo isotropico.
                    dispercion_isotropica(&rng1,particula->u,particula->v,particula->w,new_director);
                    // actualizamos el vector director
					particula->u = new_director[0];					 
					particula->v = new_director[1];
					particula->w = new_director[2];
					particula->interactua = 1;
                }
			}
		}
		// AL finalizar el grupo estadistico acumulamos los resultados en los
		// arreglos respectivos.
		for (int i=0; i<N_REGION+1; i++) 
		{

			esfera->particulas_depositadas[i] += particulas_depositadas[i];
			esfera->particulas_depositadas_pow2[i] += (particulas_depositadas[i]*particulas_depositadas[i]);
		}
		
		// Anulamos pdep y comenzamos el nuevo grupo estadistico.
        memset(particulas_depositadas, 0.0, (N_REGION+1)*sizeof(double));
        }
    //------------------------Fin del Loop principal------------------------------------------------------------


    //------------------------Analisis de los resultados--------------------------------------------------------
	

    // Analisis de resultados. Obtenemos el promedio y error asociado a las
    // fracciones de reflexion, absorcion y transmision.
    
    double *updep = (double*) malloc((N_REGION+1)*sizeof(double));
    for (int i=0; i<N_REGION+1; i++) {
    	particulas_depositadas[i] = 1.0*esfera->particulas_depositadas[i]/nbatch; 
    	updep[i] = (1.0*esfera->particulas_depositadas_pow2[i] - pow(1.0*particulas_depositadas[i], 2.0)*nbatch);
        updep[i] /= (nbatch*(nbatch-1));
        updep[i] = sqrt(1.0*updep[i]);
          
    }
    
    double pdep_abs = 0.0;
    double updep_abs = 0.0;
    for (int i=0; i<=N_REGION-1; i++) {
        pdep_abs += particulas_depositadas[i];
        updep_abs += pow(updep[i], 2.0);
    }
    
    updep_abs = sqrt(updep_abs);
    printf("Probabilidad de absorcion en nucleo : %6.5f +/- %6.5f \n",
           1.0*particulas_depositadas[0]/nperbatch, 1.0*updep[0]/nperbatch);
    printf("Probabilidad de absorcion : %6.5f +/- %6.5f \n",
           1.0*pdep_abs/nperbatch, updep_abs/nperbatch);
    printf("Probabilidad de transmision : %6.5f +/- %6.5f \n",
           1.0*particulas_depositadas[N_REGION]/nperbatch, 1.0*updep[N_REGION]/nperbatch);
    
    // Limpieza
    free(esfera);
    free(particula);
    free(updep);
    fclose(fp);




return 0;
}






























