#
#  plot.py
#
#  Creacion de graficos a partir de los resultados entregados por main_implicit.c
#
#  Creado por Michael Haacke Concha y Vicente Salinas Froemel el 15/11/18.
#
#  Requisitos: Python 2 junto con los paquetes numpy y matplotlib.
#  Como correr: Escribir en el terminal "python plot.py"

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation

# Funcion para abrir archivos de datos separados por comas como listas de listas.
def Open(file):
	try: 
		file = open(file,'r')
		content = file.readlines()
		file.close()
		data = []
		for i in range(len(content)):
			content[i] = content[i].replace('\n','')
			line = content[i].split(',')
			data.append(line)
		return data
	except IOError:
		print('error: file not found or disk full')

# Funcion para transponer una lista de listas.
def Transpose(data):
	tdata = []
	rows = len(data)
	cols = len(data[0])
	return [[data[k][j] for k in range(rows)] for j in range(cols)]

# Funcion para actualizar los datos de la animacion.
def animate(i):
    line.set_ydata(data[i+1])
    return line,

# Funcion usada luego para inicializar la animacion.
def init():
    line.set_ydata(np.ma.array(x, mask=True))
    return line,

def Animation(title):
	global line,x,data
	# Graficar los datos
	fig, ax = plt.subplots()
	x = data[0]
	y = data[1]

	# Convertir los datos a punto flotante.
	x = np.array(x).astype(float)
	y = np.array(y).astype(float)
	lasty = np.array(data[300]).astype(float)

	# Determinar las proporciones de los limites del grafico.
	Maxx,Minx,Maxy,Miny = max(x),min(x),max(lasty),min(lasty)
	propx = (Maxx-Minx)/20.
	propy = (Maxy-Miny)/20.
	lim1 = Minx-propx
	lim2 = Maxx+propx
	lim3 = Miny-propy
	lim4 = Maxy+propy

	# Definir propiedades del grafico.
	line, = ax.plot(x,y,color = 'darkred',marker = '',markerfacecolor = 'darkred',markersize = 5,linestyle = '-',linewidth = 1)
	ax.set_xlim((lim1,lim2))
	ax.set_ylim((lim3,lim4))
	ax.set_title(title,fontsize = '34',color = 'black')	
	ax.set_xlabel('Distancia',fontsize = '30',color = 'black')	
	ax.set_ylabel('Temperatura',fontsize = '30',color = 'black')
	ax.grid(color = '0.5', linestyle = ':', linewidth = 1)
	ani = animation.FuncAnimation(fig, animate, np.arange(1, 300), init_func=init, interval=50, blit=True)
	plt.show()

data1 = Transpose(Open('output1.txt'))	# Abrir el archivo con los datos del metodo explicito.
data2 = Transpose(Open('output2.txt'))	# Abrir el archivo con los datos del metodo implicito.
line,x,data = None,None,data1		# Defino algunas variables globales que seran usadas en la funcion Animation
Animation('Metodo explicito eRK-4')	# Animar los datos del metodo explicito.
data = data2
Animation('Metodo implicito')		# Animar los datos del metodo implicito.




