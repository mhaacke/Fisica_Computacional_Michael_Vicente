#
#  plot.py
#
#  Creacion de graficos a partir de los resultados entregados por main.c
#
#  Creado por Michael Haacke Concha y Vicente Salinas Froemel el 11/10/18.
#
#  Requisitos: Python junto con los paquetes numpy y matplotlib.
#  Como correr: Escribir en el terminal "python plot.py"

import numpy as np
import matplotlib.pyplot as plt

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

# Funcion para graficar los datos.
def Plot(x,y,xlabel='x',ylabel='y',title='',color='darkblue'):
	# Convertir los datos a punto flotante.
	x = np.array(x).astype(float)
	y = np.array(y).astype(float)
        
	# Determinar las proporciones de los limites del grafico.
	Maxx,Minx,Maxy,Miny = max(x),min(x),max(y),min(y)
	propx = (Maxx-Minx)/20.
	propy = (Maxy-Miny)/20.
	lim1 = Minx-propx
	lim2 = Maxx+propx
	lim3 = Miny-propy
	lim4 = Maxy+propy

	# Plotear
	plt.plot(x,y,color = color,marker = '',markerfacecolor = color,markersize = 5,linestyle = '-',linewidth = 1)
	plt.xlim((lim1,lim2))
	plt.ylim((lim3,lim4))
	plt.xlabel(xlabel,fontsize = '30',color = 'black')
	plt.ylabel(ylabel,fontsize = '30',color = 'black')
	plt.title(title,fontsize = 34,color = 'black')
	plt.xticks(fontsize = 16,color = 'black')
	plt.yticks(fontsize = 16,color = 'black')
	plt.grid(color = '0.5', linestyle = ':', linewidth = 1)
	plt.figure(1).set_facecolor('white')
	plt.show()
	
data1 = Transpose(Open('res_erk4.txt'))       # Abrir el archivo con los angulos y momentos.
data2 = Transpose(Open('res_erk4_car.txt'))   # Abrir el archivo con las coord. cartesianas.

# Graficar los datos
Plot(data1[1],data1[2],xlabel='phi1',ylabel='phi2')
Plot(data1[3],data1[4],xlabel='p1',ylabel='p2')
Plot(data2[1],data2[2],ylabel='z')
Plot(data2[3],data2[4],ylabel='z')



