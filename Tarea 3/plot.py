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
    line.set_ydata(np.array(data[i+1]).astype(float))  # update the data
    return [line,]

def init():
    line.set_ydata(np.array(data[1]).astype(float))
    return [line,]

data = Transpose(Open('output.txt'))

# Graficar los datos
fig, ax = plt.subplots()
#plt.switch_backend("TkAgg")
x = data[0]
y = data[2]
# Convertir los datos a punto flotante.
x1 = np.array(x).astype(float)
y1 = np.array(y).astype(float)
lasty = np.array(data[299]).astype(float)
# Determinar las proporciones de los limites del grafico.
Maxx,Minx,Maxy,Miny = max(x1),min(x1),max(lasty),min(lasty)
propx = (Maxx-Minx)/20.
propy = (Maxy-Miny)/20.
lim1 = 0
lim2 =20
lim3 =0
lim4 =4
line, = ax.plot(x,y,color = 'darkblue',marker = '',markerfacecolor = 'darkblue',markersize = 5,linestyle = '-',linewidth = 1)
ax.set_xlim((lim1,lim2))
ax.set_ylim((lim3,lim4))
ax.set_xlabel('x',fontsize = '30',color = 'black')
ax.set_ylabel('y',fontsize = '30',color = 'black')
# ax.grid(color = '0.5', linestyle = ':', linewidth = 1)
# ani = animation.FuncAnimation(fig, animate, np.arange(1, 299), init_func=init, interval=50, blit=True)
plt.show()
print(fig,ax)

# import numpy as np
# import matplotlib.pyplot as plt
# import matplotlib.animation as animation

# fig, ax = plt.subplots()

# #x = np.arange(0, 2*np.pi, 0.01)
# # line, = ax.plot(x, data[1])


# def init():  # only required for blitting to give a clean slate.
#     line.set_ydata(data[1])
#     return line,


# def animate(i):
#     line.set_ydata(data[i+1])  # update the data.
#     return line,


# ani = animation.FuncAnimation(
#     fig, animate,np.arange(1, 299), init_func=init, interval=2, blit=True, save_count=50)

# # To save the animation, use e.g.
# #
# # ani.save("movie.mp4")
# #
# # or
# #
# # from matplotlib.animation import FFMpegWriter
# # writer = FFMpegWriter(fps=15, metadata=dict(artist='Me'), bitrate=1800)
# # ani.save("movie.mp4", writer=writer)

# plt.show()




