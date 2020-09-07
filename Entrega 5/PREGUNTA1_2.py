# -*- coding: utf-8 -*-
"""
Created on Sun Sep  6 19:22:19 2020

@author: pauli
"""

#  PREGUNTA 1 Y 2:
from scipy.integrate import odeint
from time import perf_counter
from leer_eof import leer_eof
from matplotlib import pyplot as plt
import math
import numpy as np

archivo = leer_eof("S1A_OPER_AUX_POEORB_OPOD_20200821T121202_V20200731T225942_20200802T005942.EOF")
tiempo = archivo[0]
x = archivo[1]
y = archivo[2]
z = archivo[3]
vx = archivo[4]
vy = archivo[5]
vz = archivo[6]
z0 = [x[0],y[0],z[0],vx[0],vy[0],vz[0]]
zf = [x[-1],y[-1],z[-1],vx[-1],vy[-1],vz[-1]]

#  DATOS INCIALES:
mt=(5.972*(10**24)) #masa de la tierra: kg
# Constante gravitacional:
G=(6.67408*(10**-11))# m3/s2kg
# Velocidad angular de la rotación de la tierra
om =(7.2921*(10**-5))# rad/s

# Se utiliza la mismas funciones "satelite" y "eulerint" utilizada en entregas pasadas:

def satelite(z,t): 
    zp = np.zeros(6)
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

    # r= distancia de la tierra al satélite:   
    r_3= (math.sqrt (np.dot(z1,z1)))**3
    
    # El lado derecho de la ecuacion queda:
    lado_derecho1=np.dot(((-G*mt)/r_3),z1)
    lado_derecho2= np.dot(- R.T , np.dot(R2p,z1)+2*np.dot(Rp,z2))
    
    Z2p=lado_derecho1+lado_derecho2
    zp[0:3]=z2
    zp[3:6]=Z2p  
    return zp

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

delta = tiempo[-1]
t = np.linspace(0, delta, len(tiempo))
t1=perf_counter()

sol_odeint= odeint(satelite,z0,t)
t2=perf_counter()

sol_eulerint = eulerint(satelite,z0,t,1)
t3=perf_counter()

tiempo_odeint = t2-t1
tiempo_euler = t3-t2

x_od = sol_odeint[:,0]
y_od = sol_odeint[:,1]
z_od = sol_odeint[:,2]
vx_od = sol_odeint[:,3]
vy_od = sol_odeint[:,4]
vz_od = sol_odeint[:,5]


print (f"tiempo odeint = {tiempo_odeint} s")
print (f"tiempo eulerint = {tiempo_euler} s")

# Se grafica la posición del satélite
plt.figure(1)
plt.subplot(3,1,1)
plt.title("Posición del Sentinel")
plt.ylabel("X(t) [KM]")
plt.yticks((-5000000, 0 ,5000000 ),("-5000","0","5000"))
plt.xticks ((0,18000,36000,54000,72000,90000),("0","5","10","15","20","25"))
plt.plot(tiempo,x,c="blue")
plt.plot(tiempo,x_od, c="orange")
plt.tight_layout()

plt.subplot(3,1,2)
plt.ylabel("Y(t) [KM]")
plt.yticks((-5000000, 0 ,5000000 ),("-5000","0","5000"))
plt.xticks ((0,18000,36000,54000,72000,90000),("0","5","10","15","20","25"))
plt.plot(tiempo,y,c="blue")
plt.plot(tiempo,y_od, c="orange")
                                   
plt.subplot(3,1,3)
plt.ylabel("Z(t)[KM]")
plt.xlabel("Tiempo, t [horas]")
plt.yticks((-5000000, 0 ,5000000 ),("-5000","0","5000"))
plt.xticks ((0,18000,36000,54000,72000,90000),("0","5","10","15","20","25"))
plt.plot(tiempo,z,c="blue")
plt.plot(tiempo,z_od, c="orange")
plt.tight_layout()
plt.savefig("Preg1_Posición")

x_eu=sol_eulerint[:,0]
y_eu=sol_eulerint[:,1]
z_eu=sol_eulerint[:,2]

grad_x = np.gradient(vx_od,tiempo)
grad_y = np.gradient(vy_od,tiempo)
grad_z = np.gradient(vz_od,tiempo)

grad_realx = np.gradient(vx,tiempo)
grad_realy = np.gradient(vy,tiempo)
grad_realz = np.gradient(vz,tiempo)

delta_odeint = np.sqrt((x_od-x)**2+(y_od-y)**2+(z_od-z)**2)
delta_euler = np.sqrt((x_eu-x)**2+(y_eu-y)**2+(z_eu-z)**2)
print (f"La deriva de eulerint vs odeint: {delta_euler[-1]/1000-delta_odeint[-1]/1000} Km")
hora= t/3600
nueva = np.zeros(len(tiempo))
for i in range(len(tiempo)):
    nueva[i] = math.sqrt(np.dot((sol_odeint[i,:3] - sol_eulerint[i,:3]), (sol_odeint[i,:3] - sol_eulerint[i,:3])))
cuadrado=sol_odeint[-1,:3]
final = math.sqrt(np.dot(cuadrado,cuadrado))
error = np.round_(nueva[-1]/final,1)
print (f"Error = {error*100} %")


plt.figure(2)
plt.plot(hora,delta_odeint/1000,label="Odeint", c="blue")
plt.plot(hora,delta_euler/1000,label="Eulerint",c="green")
plt.title(f" Distancias: $\\delta_{{odeint}} = {delta_odeint[-1]/1000:.1f}$ [Km] , $\\delta_{{eulerint}} = {delta_euler[-1]/1000:.1f}$ [Km]")
plt.ylabel(" $\\delta$ [KM]")
plt.xlabel("Tiempo[hora]")
plt.tight_layout()
plt.legend()
plt.savefig("Preg2_Deriva_EulerintVSOdeint")


plt.figure(3)
plt.plot(hora,delta_euler/1000,label="Eulerint", c="green")
plt.title(f" Distancia entre posición real y predicha $\\delta_{{max}} = {delta_euler[-1]/1000:.1f}$ [Km]")
plt.ylabel(" $\\delta$ [KM]")
plt.xlabel("Tiempo[hora]")
plt.tight_layout()
plt.legend()
plt.savefig("Preg2_Deriva_Eulerint")

plt.figure(4)
plt.plot(hora,delta_odeint/1000,label="Odeint",c="blue")
plt.title(f" Distancia entre posición real y predicha $\\delta_{{max}} = {delta_odeint[-1]/1000:.1f}$ [Km]")
plt.ylabel(" $\\delta$ [KM]")
plt.xlabel("Tiempo[hora]")
plt.tight_layout()
plt.legend()
plt.savefig("Preg2_Deriva_Odeint")
plt.show()


