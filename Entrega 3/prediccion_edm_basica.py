# -*- coding: utf-8 -*-
"""
Created on Fri Aug 28 15:39:37 2020

@author: pauli
"""

import math
import numpy as np
import scipy as sp
from scipy.integrate import odeint

# -----------------------------------------------
# Entrega 3: Predicción de vector de estado 
# de la EDM básica
# ------------------------------------------


#  DATOS INCIALES:
hora= 3600 # seg
km= 1000 #m
# masa de la tierra
mt=(5.972*(10**24)) #kg

# Distancia de la tierra al satélite:
d= 700* km
radio= 6371*km #km

# Constante gravitacional:
G=(6.67408*(10**-11))# m3/s2kg
# Velocidad angular de la rotación de la tierra
om =(7.2921*(10**-5))# rad/s

FdMax= (G*mt)/radio**2

# Vector de estado:
    # z[0] -> x
    # z[1] -> y
    # z[2] -> z   
    # z[3] -> vx
    # z[4] -> vy
    # z[5] -> vz

def satelite(z,t): 
    zp = sp.zeros(6)
    # Se crea un vector de 6x1 para trabajar:
    c= np.cos(om*t)
    s= np.sin (om*t)
    R=np.array([[c,-s,0],[s, c,0], [ 0, 0, 1]])
    
    #Primera derivada de R(teta):
    Rp=om*np.array([[-s,-c,0], [c, -s,0],[ 0, 0, 0]])
    
    # Segunda derivada de R(teta):
    R2p=(om**2)*np.array([[-c,  s,0], [-s, -c,0],[0, 0, 0]])
    
    z2=z[3:6]
    z1=z[0:3]

# Se toma en cuenta la ecuación presentada en el video:
# Rp= R punto= 1ra derivada de R
# R2p= R 2 puntos= 2da derivada de R
# Z2p= ((-G*mt)/r^3)Z1 - R^t(R2p*Z1+2Rp*Z2)

    # r= distancia de la tierra al satélite:   
    r_3= (math.sqrt (np.dot(z1,z1)))**3
    
    # El lado derecho de la ecuacion queda:
    lado_derecho1=np.dot(((-G*mt)/r_3),z1)
    lado_derecho2= np.dot(- R.T , np.dot(R2p,z1)+2*np.dot(Rp,z2))
    
    Z2p=lado_derecho1+lado_derecho2
    zp[0:3]=z2
    zp[3:6]=Z2p  
    return zp



import datetime as dt
utc_EOF_format=   "%Y-%m-%dT%H:%M:%S.%f"
t1= dt.datetime.strptime ("2020-07-31T22:59:42.000000",utc_EOF_format)
t2= dt.datetime.strptime ("2020-08-02T00:59:42.000000",utc_EOF_format)
# Tomando en cuenta las dos fechas :
intervalo= t2 - t1 
intervalo_en_segundos= intervalo.total_seconds()
print (f"Iintervalo = {intervalo} s")
print (f"Intervalo_en_segundos = {intervalo_en_segundos} s")


# Parámetros iniciales:
x_i= -294271.795997
y_i= -2511619.758656
z_i= 6598339.200344

vx_i= -2401.659901
vy_i= 6760.253877
vz_i= 2460.740181

# Parámetros finales:
x_f= -1753088.623736
y_f= -6852318.248084
z_f= -241271.374219

vx_f= -1594.856888
vy_f= 155.635645
vz_f= 7426.026842


# Vector de tiempo:
t= np.linspace (0, intervalo_en_segundos, 9361)
# Posición inicial:
x0= radio + d
vt= 6820
# Vector de los parámetros iniciales del satélite:
z0= np.array ([x_i, y_i, z_i, vx_i, vy_i, vz_i])
sol = odeint (satelite,z0,t)

x= sol [:,0:3]
# Vector de los parámetros finales del satélite:
final= sp.array ([x_f, y_f, z_f, vx_f,vy_f,vz_f]) - sol [-1]

# Diferencia entre los parámetros finales e iniciales:
deriva= math.sqrt(final[0]**2 + final[1]**2 + final [2]**2)

print (" La diferencia entre la ubicación inicial y final del satélite :",deriva, "m")








