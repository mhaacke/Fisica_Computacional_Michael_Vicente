### Proyecto Final Del Curso

Autores:

	Nombre: Michael Haacke Concha

	Nombre: Vicente Salinas Froemel

	Para compilar el codigo se debe correr 

	gcc -o main main.c

Luego para correr el codigo se debe realizar 

	./main

El programa consta de un Modelo transporte de particulas sobre un nucleo fisionable de forma esferica con varios tipo de materiales mediante el metodo Monte Carlo. En especifico determina el porcentage de particulas transmitida y absorbidas si estas vienen del nucleo de la esfera y escribe en un archivo todas las historia, se calcula el tiempo de vuelo como la exponecial de -∑_t y luego que tipo de interraccion es la que sucede, Se ocupa un generador de numeros pseudo-random del estilo "Generadores lineales congruenciales" con valores obtenidos en Numerical Recipes in C. Parte del codigo es obtenido de la clase de fisica computacional del profesor Edgardo Dörner 
