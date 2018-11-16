### La ecuación de calor dependiente del tiempo


Nombre: Michael Haacke Concha

Nombre: Vicente Salinas Froemel

Para compilar el codigo se debe correr 

	gcc -o main_implict main_implicit.c
	gcc -o main_explici main_explici.c

Luego para correr el codigo se debe realizar 

	./main_implict
	./main_explici


Como requerimientos, para graficar los datos se debe tener Python 2 junto con los paquetes numpy y matplotlib.

Luego se debe correr plot.py para realizar las animaciones. Esto se hace escribiendo

	python plot.py

La animacion tiene problemas en mac os.

Por una parte main_explicit utiliza un metodo explicito, el metodo de Runge-Kutta de orden 4, para resolver el problema de la ecuación de calor dependiente del tiempo bajo un sumidero de forma gaussiana. Por otro lado, main_impicit utiliza un metodo implicito, el metodo de Crank-Michelson, en el cual se resuelve un sistema de ecuaiones tridiagonal. El codigo de este último metodo fue obtendio de “Numerical recipes in C.” Cambridge University Press 1 (1988).
	
