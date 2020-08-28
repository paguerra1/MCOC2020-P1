# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

import numpy as np
import scipy as sp
from scipy.integrate import odeint
import matplotlib.pylab as plt

# masa de la tierra
mt=(5.972*(10**24)) #kg

# Distancia de la tierra al satélite:
d= 700000 #m

radio_tierra=6371   # km   
radio_m= radio_tierra*1000 #m
#  = 6,67·10–11 N m^2/kg^2
# 1N= 1kg*m/s^2
# Constante gravitacional:
G=(6.67*(10**-11))# m3/s2kg
# Velocidad angular de la rotación de la tierra
om =(7.27*(10**-5))# rad/s


def matrices(om,t):
    
    R=np.array([[np.cos(om*t),-np.sin(om*t),0],[np.sin(om*t), np.cos(om*t),0], [ 0, 0, 1]])
    
    #Primera derivada de R(teta):
    Rp=om*np.array([[-np.sin(om*t),-np.cos(om*t),0], [np.cos(om*t), -np.sin(om*t),0],[ 0, 0, 0]])
    
    # Segunda derivada de R(teta):
    R2p=(om**2)*np.array([[-np.cos(om*t),  np.sin(om*t),0], [-np.sin(om*t), -np.cos(om*t),0],[0, 0, 0]])
    if om==0:
        V=[R,np.zeros((3,3)),np.zeros((3,3))]
        return V
    V=(R,Rp,R2p)
    return V

# Vector de estado:
    # z[0] -> x
    # z[1] -> y
    # z[2] -> z   
    # z[3] -> vx
    # z[4] -> vy
    # z[5] -> vz

def satelite(z,t): 
    # Se crea un vector de 6x1 para trabajar:
    zp = sp.zeros(6)

    zp[0] = z[3]
    zp[1] = z[4]
    zp[2] = z[5]
   
# Se toma en cuenta la ecuación presentada en el video:
# Rp= R punto= 1ra derivada de R
# R2p= R 2 puntos= 2da derivada de R
# Z2p= ((-G*mt)/r^3)Z1 - R^t(R2p*Z1+2Rp*Z2)
# Z2=z[3:6]
# Z1=z[0:3]
    
    # r= distancia de la tierra al satélite:   
    # r_3=(((np.dot((z[0:3]),(z[0:3])))(1/2))**3
    r_3=(radio_m+d)**3
    # El lado derecho de la ecuacion queda:
    lado_derecho1=np.dot(((-G*mt)/r_3),z[0:3])
    
    V=matrices(om, t)
    R=V[0]
    Rp=V[1]
    R2p=V[2]
    lado_derecho2= np.dot(-np.transpose(R) , np.dot(R2p,z[0:3])+2*np.dot(Rp,z[3:6]))
    
    Z2p=lado_derecho1+lado_derecho2
    zp[3:6]=Z2p   
    return zp

hora=3600  #seg
t = sp.linspace(0,hora*8 ,(hora*8)+1)

# Velocidad:

vi = (23000*1000)/hora

# Coordenadas iniciales:[xi=radio_m+d,yi=0,z0=0,vx=0,vy=vi,vz=0].
z0 = sp.array([radio_m+d ,0,0,0,vi,0])
# Se resuelve la EDO:
sol = odeint(satelite,z0,t)
x = sol[:,0]
y = sol[:,1]
z = sol[:,2]

plt.figure(1)

# Se grafica el eje x(t):
plt.subplot(3,1,1)
plt.plot(t,x)
yticks1=[10000000,0,-10000000]
yticks_text1=["-10.000 km","0 km","10.000 km"]
xticks_1=[0,2500,5000,10000,15000,20000,25000,30000]
plt.yticks(yticks1,yticks_text1)
plt.xticks(xticks_1,[])
plt.title("Historias de tiempo\n x(t)")
plt.ylabel("Distancia (km)")

# Se grafica el eje y(t):
plt.subplot(3,1,2)
plt.plot(t,y,color="red")
yticks2=[10000000,0,-10000000]
yticks_text2=["-10000 km","0 km","10000 km"]
xticks_2=[0,2500,5000,10000,15000,20000,25000,30000]
plt.yticks(yticks2,yticks_text2)
plt.xticks(xticks_2,[])
plt.title("y(t)")
plt.ylabel("Distancia (km")

# Se grafica el eje z(t):
plt.subplot(3,1,3)
plt.plot(t,z,color="cyan")
yticks3=[0.05,0,-0.05]
yticks_text3=["1 km","0 km","-1 km"]
xticks_3=[0,2500,5000,10000,15000,20000,25000,30000]
plt.yticks(yticks3,yticks_text3)
plt.xticks(xticks_3,xticks_3,rotation=45)
plt.title("z(t)")
plt.ylabel("Distancia (km)")
plt.xlabel("Tiempo transcurrido" )
plt.tight_layout()

plt.show()

# ---------------------------------------------------------
# Gráfico de la distancia desde el satélite a la tierra
# en el periodo de tiempo transcurrido
# ----------------------------------------------------------
plt.figure(1)
norma=((x**2+y**2)**(0.5))
plt.plot(t,norma)
tierra=radio_m+80000

yticks4=[8000000,5000000]
yticks_text4=["80000 km", "5000 km"]
xticks_4=[0,2500,5000,10000,15000,20000,25000,30000]

plt.yticks(yticks4, yticks_text4)
plt.xticks(xticks_4,xticks_4,rotation=45)
plt.ylabel("Distancia(Km)")
plt.xlabel("Tiempo transcurrido (s)")
# Tierra:
plt.hlines(radio_m,0,hora*8,color="green")
# Atmosfera:
plt.hlines((radio_m+80000),0,hora*8,color="cyan")
plt.title("Distancia del satélite V/S tiempo")
plt.tight_layout()
plt.show()


plt.figure(3)
P= np.linspace(0, 2*np.pi , 1001)
xx1= tierra *np.cos(P)
yy1= tierra *np.sin(P)
plt.plot(xx1,yy1, color="cyan")
plt.xlabel("X")
plt.ylabel("Y")
plt.plot(x,y)
plt.gca().set_aspect("equal")
plt.title("Trayetoria del satélite")
plt.tight_layout()
plt.show()















