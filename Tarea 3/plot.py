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

def animate(i):
    line.set_ydata(data[i+1])  # update the data
    return line,

def init():
    line.set_ydata(np.ma.array(x, mask=True))
    return line,

data = Transpose(Open('output.txt'))

# Graficar los datos
fig, ax = plt.subplots()
#plt.switch_backend("TkAgg")
x = data[0]
y = data[1]
# Convertir los datos a punto flotante.
x = np.array(x).astype(float)
y = np.array(y).astype(float)
lasty = np.array(data[299]).astype(float)
# Determinar las proporciones de los limites del grafico.
Maxx,Minx,Maxy,Miny = max(x),min(x),max(lasty),min(lasty)
propx = (Maxx-Minx)/20.
propy = (Maxy-Miny)/20.
lim1 = Minx-propx
lim2 = Maxx+propx
lim3 = Miny-propy
lim4 = Maxy+propy
line, = ax.plot(x,y,color = 'darkblue',marker = '',markerfacecolor = 'darkblue',markersize = 5,linestyle = '-',linewidth = 1)
ax.set_xlim((lim1,lim2))
ax.set_ylim((lim3,lim4))
ax.set_xlabel('x',fontsize = '30',color = 'black')
ax.set_ylabel('y',fontsize = '30',color = 'black')
ax.grid(color = '0.5', linestyle = ':', linewidth = 1)
ani = animation.FuncAnimation(fig, animate, np.arange(1, 300), init_func=init, interval=50, blit=True)
plt.show()




