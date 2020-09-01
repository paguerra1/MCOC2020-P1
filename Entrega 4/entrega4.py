# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 19:25:24 2020

@author: pauli
"""

from scipy.integrate import odeint
import numpy as np
import matplotlib.pylab as plt

#DATOS DEL PROBLEMA:
m = 1       #kg
f = 1       #Hz
chi = 0.2
om = 2*(np.pi)*f
k = m*(om**2)
c = 2*chi*om*m
p=(1-chi**2)**0.5
# z = [x, xp]
# zp =[xp, xpp]


def eulerint(zp, z0, t, Nsubdivisiones=1):
    Nt = len(t)
    Ndim = len(np.array(z0))
    z = np.zeros((Nt,Ndim))
    z[0,:] = z0
    
    for i in range(1, Nt):
    
        t_anterior = t[i-1]
        
        dt = (t[i] - t[i-1])/Nsubdivisiones
        
        z_temp = z[i-1, :].copy()
        
        for k in range(Nsubdivisiones):
            z_temp += dt*zp(z_temp , t_anterior + k*dt)
            
        z[i,:] = z_temp
        
    return z

# z' = a*z   
# oscilador armónico: m * xpp + c*xp + k*x = 0
#  xpp= (-c*xp - k*x)/m

# cond. inicial: x(0)=1  , xp(0)=1
def z_punto(z,t):
    xp= z[1]
    x=z[0]  
    zp = np.zeros(2)
    zp[0] = z[1]
    zp[1] = -(c*xp+k*x)/m
   
    return zp

z0 = [1,1]
t = np.linspace(0,4.,100)

sol= odeint(z_punto,z0,t)
z_odein= sol[:,0]

# Para condiciones iniciales del problema se tiene: 
# z_real(t)= e^(-c/2m) * (m*cos(om*t* ) + 81sin(om*t)

z_real = np.exp((-c/(2*m))*t)*(m*np.cos(om*p*t) + ((1 + om*chi*m)/(om*(p)))*np.sin(om*p*t))

sol_1 = eulerint(z_punto,z0,t,1)
sol_10 = eulerint(z_punto,z0,t,10)
sol_100 = eulerint(z_punto,z0,t,100)
z_euler1=  sol_1[:,0]
z_euler10=  sol_10[:,0]
z_euler100=  sol_100[:,0]

# Se grafica:
fig = plt.figure() 
plt.plot(t,z_odein, color="b", label="Odeint")
plt.plot(t,z_real, color="k",linewidth=2,label="Real")
plt.plot(t,z_euler1, color="green", linestyle="--", label="Euler N=1")
plt.plot(t,z_euler10, color="red", linestyle="--",  label="Euler N=10")
plt.plot(t,z_euler100, color="orange", linestyle="--",  label="Euler N=100")


plt.title("Oscilador Armónico")
plt.xlabel("Tiempo [s]")
plt.ylabel("u[m]")

plt.legend()
plt.show()