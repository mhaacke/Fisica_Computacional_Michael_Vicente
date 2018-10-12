### Tareas Fisica Computacional 

Nombre: Michael Haacke Concha

Nombre: Vicente Salinas Froemel

Aqui se tienen dos archivos, main.c y plot.py. Estos sirven para resolver numericamente el problema de los dos pendulos mediante el metodo Runge-Kutta de orden cuatro, y graficar los resultados, respectivamente. El primero de ellos es un codigo escrito en C, por lo que para correrlo se debe compilar (gcc -o main.x main.c) y luego correr (./main.x) el archivo. Al correrlo entregara dos archivos, res_erk4.txt y res_erk4_car.txt, conteniendo los valores del angulo y del momento generalizado en funcion del tiempo, y las coordanadas cartesianas de las masas en funcion del tiempo, respectivamente. El segundo es un codigo escrito en Python, y usa los paquetes numpy y matplotlib. Para correrlo basta escribir python plot.py en el terminal.

El programa fue comparado con los resultados mostrados en el libro "Basic Concepts in Computational Physics", de Benjamin A. Stickler y Ewald Schachinger, y se noto que se llegan a resultados identicos a simple vista.

