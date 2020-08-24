# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 10:25:10 2020

@author: pauli
"""
# ENTREGA 1- INTEGRACION DE ECUACIONES DIFERENCIALES 

import scipy as sp
from scipy.integrate import odeint

# Parametros base:
g= 9.81  # m/s**2
cm=0.01
inch=2.54*cm
 
# coeficiente de arrastre:
ro=1.225  # kg/m**3 
cd=0.47
D= 8.5 *inch
r=D/2

A= sp.pi * r**2
CD= (1/2)*ro*cd*A

# masa= kilos
m=15.


#  Funcion a integrar:
    # z es el vector de estado
    # z=[x,y,vx,vy]
    # dz/dt= bala (z,t)
     #  dz1/dt =z2
     #        [ z2       ]
     # dz/dt= [          ]
     #        [ FD/m   -g]
    
    
# Vector de estado:
    # z[0] -> x
    # z[1] -> y
    # z[2] -> vx
    # z[3] -> vy
    
def bala(z,t):
    zp= sp.zeros(4) 
    zp[0]= z[2]
    zp[1]= z[3]
    # vector velocidad
    v = z[2:4]  # saca los ltimos dos componentes
    v[0]=v[0]-V

    v2= sp.dot(v,v)
    vnorm= sp.sqrt(v2)
    
    #  fuerza de arrastre:
    FD= - CD*v2* (v / vnorm)
    zp[2]= FD[0]/ m
    zp[3]= FD [1] /m - g
    
    return zp

V=[0,10,20]

for i in V:
    V= i
    # vector de tiempo
    t= sp.linspace (0, 30, 1001)
    
    # velocidad inicial en km/h:
    vi=100*(1000/3600)
    
    # parte en el origen, tiene vx=vy=2m/s
    z0= sp.array([ 0, 0, vi, vi])
    
    sol= odeint(bala, z0, t)
    
    import matplotlib.pylab as plt
    x=sol[:,0]
    y=sol[:,1]
    
    plt.figure(1)
    plt.plot(x, y, label = f"V ={i} m/s")
    plt.axis([0,150,0,50])



plt.title("Trayectoria para distintos vientos")
plt.ylabel("Y (m)")
plt.xlabel("X (m)")
plt.grid()
plt.legend()
plt.tight_layout() 
plt.savefig("Trayectoria para distintos vientos.png")



