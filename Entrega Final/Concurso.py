# -*- coding: utf-8 -*-
"""
Created on Fri Sep 11 19:19:21 2020

@author: pauli
"""


from scipy.integrate import odeint
from leer_eof import leer_eof
from sys import argv
from datetime import datetime,timedelta
import math
import numpy as np


fname= argv[1]
t, x, y, z, vx, vy, vz =  leer_eof(fname)
fname_out = fname.replace(".EOF", ".PRED")
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
J2 = 1.75553e10*(1000**5) # km5 -s2 -> m5
J3 = -2.61913e11*(1000**6) # km6 s-2 -> m6
Rt = 6378136.3   # radio de la tierra
mu = 398600.440*(1000.)**3 #G*Mtierra



# Se aplicará el mismo código utilizado por el profesor 
#Que permite obtener los coeficientes hasta  n = m = 8 de la expansión en armónicos esféricos del geopotencial terrestre
# Solo se aplicaron esta cantidad de términos


#Terminos zonales  ("Jn"), son los C[n,0], notar que S[n,0] = 0

Nt = 8
C = np.zeros((Nt+1,Nt+1))
S = np.zeros((Nt+1,Nt+1))

# n,m     n m       Jn
C[2,0], S[2,0] =   -0.10826360229840E-02   ,     0.0
C[3,0], S[3,0] =    0.25324353457544E-05   ,     0.0
C[4,0], S[4,0] =    0.16193312050719E-05   ,     0.0
C[5,0], S[5,0] =    0.22771610163688E-06   ,     0.0
C[6,0], S[6,0] =   -0.53964849049834E-06   ,     0.0
C[7,0], S[7,0] =    0.35136844210318E-06   ,     0.0
C[8,0], S[8,0] =    0.20251871520885E-06   ,     0.0

#Terminos teserales Cmn y Smn
# n,m     n m       Cmn                          Smn
C[2,1], S[2,1] =   -0.24140000522221E-09   ,     0.15430999737844E-08
C[3,1], S[3,1] =    0.21927988018965E-05   ,     0.26801189379726E-06
C[4,1], S[4,1] =   -0.50872530365024E-06   ,    -0.44945993508117E-06
C[5,1], S[5,1] =   -0.53716510187662E-07   ,    -0.80663463828530E-07
C[6,1], S[6,1] =   -0.59877976856303E-07   ,     0.21164664354382E-07
C[7,1], S[7,1] =    0.20514872797672E-06   ,     0.69369893525908E-07
C[8,1], S[8,1] =    0.16034587141379E-07   ,     0.40199781599510E-07 
C[2,2], S[2,2] =    0.15745360427672E-05   ,    -0.90386807301869E-06
C[3,2], S[3,2] =    0.30901604455583E-06   ,    -0.21140239785975E-06
C[4,2], S[4,2] =    0.78412230752366E-07   ,     0.14815545694714E-06
C[5,2], S[5,2] =    0.10559053538674E-06   ,    -0.52326723987632E-07
C[6,2], S[6,2] =    0.60120988437373E-08   ,    -0.46503948132217E-07
C[7,2], S[7,2] =    0.32844904836492E-07   ,     0.92823143885084E-08
C[8,2], S[8,2] =    0.65765423316743E-08   ,     0.53813164055056E-08
C[3,3], S[3,3] =    0.10055885741455E-06   ,     0.19720132389889E-06
C[4,3], S[4,3] =    0.59215743214072E-07   ,    -0.12011291831397E-07
C[5,3], S[5,3] =   -0.14926153867389E-07   ,    -0.71008771406986E-08
C[6,3], S[6,3] =    0.11822664115915E-08   ,     0.18431336880625E-09
C[7,3], S[7,3] =    0.35285405191512E-08   ,    -0.30611502382788E-08
C[8,3], S[8,3] =   -0.19463581555399E-09   ,    -0.87235195047605E-09
C[4,4], S[4,4] =   -0.39823957404129E-08   ,     0.65256058113396E-08
C[5,4], S[5,4] =   -0.22979123502681E-08   ,     0.38730050770804E-09
C[6,4], S[6,4] =   -0.32641389117891E-09   ,    -0.17844913348882E-08
C[7,4], S[7,4] =   -0.58511949148624E-09   ,    -0.26361822157867E-09
C[8,4], S[8,4] =   -0.31893580211856E-09   ,     0.91177355887255E-10
C[5,5], S[5,5] =    0.43047675045029E-09   ,    -0.16482039468636E-08
C[6,5], S[6,5] =   -0.21557711513900E-09   ,    -0.43291816989540E-09
C[7,5], S[7,5] =    0.58184856030873E-12   ,     0.63972526639235E-11
C[8,5], S[8,5] =   -0.46151734306628E-11   ,     0.16125208346784E-10
C[6,6], S[6,6] =    0.22136925556741E-11   ,    -0.55277122205966E-10
C[7,6], S[7,6] =   -0.24907176820596E-10   ,     0.10534878629266E-10
C[8,6], S[8,6] =   -0.18393642697634E-11   ,     0.86277431674150E-11
C[7,7], S[7,7] =    0.25590780149873E-13   ,     0.44759834144751E-12 
C[8,7], S[8,7] =    0.34297618184624E-12   ,     0.38147656686685E-12 
C[8,8], S[8,8] =   -0.15803322891725E-12   ,     0.15353381397148E-12

# Se define cada uno de los "J" con los parámetros mencionados anteriormente:

J2 = -C[2,0]*mu*Rt**2
J3 = -C[3,0]*mu*Rt**3
J4 = -C[4,0]*mu*Rt**4
J5 = -C[5,0]*mu*Rt**5
J6 = -C[6,0]*mu*Rt**6
J7 = -C[7,0]*mu*Rt**7
J8 = -C[8,0]*mu*Rt**8

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
    r = math.sqrt(np.dot(z1,z1))
    
       
    # El lado derecho de la ecuacion queda:
    lado_derecho1=np.dot(((-mu)/r**2),z1)
    lado=-(mu)/r**2 * R@z1/r
    lado_derecho2= np.dot(- R.T , np.dot(R2p,z1)+2*np.dot(Rp,z2))
    
    xx=R@z1
    zz2 = xx[2]**2
    zz1= xx[1]**2
    zz0= xx[0]**2
    
  # Se incluye J2 y J3 usados en entregas anteriores:
    FJ2 = np.zeros(3)
    J22=J2*(xx)/r**7
    FJ2[0] = J22[0]*(6*zz2 -1.5*(zz0 + zz1))
    FJ2[1] = J22[1]*(6*zz2 -1.5*(zz0 + zz1))
    FJ2[2] = J22[2]*(3*zz2 -4.5*(zz0 + zz1))
    
    FJ3 = np.zeros(3)
    FJ3[0] = J3*xx[0]*xx[2]/r**9 * (10*zz2 - 7.5*(zz0 + zz1))
    FJ3[1] = J3*xx[1]*xx[2]/r**9 * (10*zz2 - 7.5*(zz0 + zz1))
    FJ3[2] = J3/r**9 * (4*zz2 *(zz2 - 3*(zz0 + zz1)) +1.5*(zz0 + zz1)**2) 
    
   # Se incluyen correcciones desde el J4 al J7:
    
    FJ4 = np.zeros(3)
    FJ4[0] = -5*xx[0]*(35*xx[2]**4/(8*(zz0 + zz1 + zz2)**2) - 15*zz2/(4*(zz0+ zz1 + zz2)) + 0.375)/(zz0 + zz1 + zz2)**(7/2) + (-35*xx[0]*xx[2]**4/(2*(zz0 + zz1 + zz2)**3) + 15*xx[0]*zz2/(2*(zz0 + zz1 + zz2)**2))/(zz0 + zz1+ zz2)**(5/2)
    FJ4[1] = -5*xx[1]*(35*xx[2]**4/(8*(zz0 + zz1 + zz2)**2) - 15*zz2/(4*(zz0 + zz1 + zz2)) + 0.375)/(zz0 + zz1 + zz2)**(7/2) + (-35*xx[1]*xx[2]**4/(2*(zz0 + zz1 + zz2)**3) + 15*xx[1]*zz2/(2*(zz0 + zz1 + zz2)**2))/(zz0 + zz1 + zz2)**(5/2)
    FJ4[2] = -5*xx[2]*(35*xx[2]**4/(8*(zz0 + zz1 + zz2)**2) - 15*zz2/(4*(zz0 + zz1 + zz2)) + 0.375)/(zz0 + zz1 + zz2)**(7/2) + (-35*xx[2]**5/(2*(zz0 + zz1 + zz2)**3) + 25*xx[2]**3/(zz0 + zz1 + xx[2]**2)**2 - 15*xx[2]/(2*(zz0 + zz1+ zz2)))/(zz0 + zz1 + zz2)**(5/2)
    
    FJ5 = np.zeros(3)
    FJ5[0] = -6*xx[0]*(63*xx[2]**5/(8*(zz0 + zz1 + zz2)**(5/2)) - 35*xx[2]**3/(4*(zz0 + zz1 + zz2)**(3/2)) + 15*xx[2]/(8*math.sqrt(zz0 +zz1 + zz2)))/(zz0 + zz1 + zz2)**4 + (-315*xx[0]*xx[2]**5/(8*(zz0 + zz1 + zz2)**(7/2)) + 105*xx[0]*xx[2]**3/(4*(zz0 + zz1+ zz2)**(5/2)) - 15*xx[0]*xx[2]/(8*(zz0 + zz1 + zz2)**(3/2)))/(zz0 + zz1 + zz2)**3
    FJ5[1] = -6*xx[1]*(63*xx[2]**5/(8*(zz0 + zz1 + zz2)**(5/2)) - 35*xx[2]**3/(4*(zz0 + zz1 + zz2)**(3/2)) + 15*xx[2]/(8*math.sqrt(zz0 + zz1 + zz2)))/(zz0 +zz1 + zz2)**4 + (-315*xx[1]*xx[2]**5/(8*(zz0 + zz1 + zz2)**(7/2)) + 105*xx[1]*xx[2]**3/(4*(zz0 + zz1 + zz2)**(5/2)) - 15*xx[1]*xx[2]/(8*(zz0 + zz1 + zz2)**(3/2)))/(zz0 + zz1 + zz2)**3
    FJ5[2] = -6*xx[2]*(63*xx[2]**5/(8*(zz0 + zz1 + zz2)**(5/2)) - 35*xx[2]**3/(4*(zz0+ zz1 + zz2)**(3/2)) + 15*xx[2]/(8*math.sqrt(zz0 + zz1 + zz2)))/(zz0 +zz1 + zz2)**4 + (-315*xx[2]**6/(8*(zz0 + zz1 + zz2)**(7/2)) + 525*xx[2]**4/(8*(zz0 + zz1 + zz2)**(5/2)) - 225*zz2/(8*(zz0 + zz1 + zz2)**(3/2)) + 15/(8*math.sqrt(zz0 + zz1 + zz2)))/(zz0 + zz1 + zz2)**3
    
    FJ6 = np.zeros(3)
    FJ6[0] = -7*xx[0]*(231*xx[2]**6/(16*(zz0 + zz1 + zz2)**3) - 315*xx[2]**4/(16*(zz0 + zz1 + zz2)**2) + 105*zz2/(16*(zz0 + zz1 + zz2)) - 0.3125)/(zz0 + zz1 + zz2)**(9/2) + (-693*xx[0]*xx[2]**6/(8*(zz0 + zz1 + zz2)**4) + 315*xx[0]*xx[2]**4/(4*(zz0+ zz1 + zz2)**3) - 105*xx[0]*zz2/(8*(zz0 + zz1 + zz2)**2))/(zz0 + zz1 + zz2)**(7/2)
    FJ6[1] = -7*xx[1]*(231*xx[2]**6/(16*(zz0 + zz1 + zz2)**3) - 315*xx[2]**4/(16*(zz0 + zz1 + zz2)**2) + 105*zz2/(16*(zz0 + zz1 + zz2)) - 0.3125)/(zz0 + zz1 + zz2)**(9/2) + (-693*xx[1]*xx[2]**6/(8*(zz0 +zz1  + zz2)**4) + 315*xx[1]*xx[2]**4/(4*(zz0 + zz1 + zz2)**3) - 105*xx[1]*zz2/(8*(zz0 + zz1 + zz2)**2))/(zz0 + zz1 + zz2)**(7/2)
    FJ6[2] = -7*xx[2]*(231*xx[2]**6/(16*(zz0 + zz1 + zz2)**3) - 315*xx[2]**4/(16*(zz0 + zz1 + zz2)**2) + 105*zz2/(16*(zz0 + zz1 + zz2)) - 0.3125)/(zz0 + zz1 + zz2)**(9/2) + (-693*xx[2]**7/(8*(zz0 + zz1 + zz2)**4) + 1323*xx[2]**5/(8*(zz0 + zz1 + zz2)**3) - 735*xx[2]**3/(8*(zz0 + zz1 + zz2)**2) + 105*xx[2]/(8*(zz0 + zz1 + zz2)))/(zz0 + zz1 + zz2)**(7/2)
    

    FJ7 = np.zeros(3)
    FJ7[0] = -8*xx[0]*(429*xx[2]**7/(16*(zz0 + zz1 + zz2)**(7/2)) - 693*xx[2]**5/(16*(zz0 + zz1 + zz2)**(5/2)) + 315*xx[2]**3/(16*(zz0 + zz1 + zz2)**(3/2)) - 35*xx[2]/(16*math.sqrt(zz0 + zz1 + zz2)))/(zz0 + zz1 + zz2)**5 + (-3003*xx[0]*xx[2]**7/(16*(zz0 + zz1 + zz2)**(9/2)) + 3465*xx[0]*xx[2]**5/(16*(zz0 + zz1 + zz2)**(7/2)) - 945*xx[0]*xx[2]**3/(16*(zz0 + zz1 + zz2)**(5/2)) + 35*xx[0]*xx[2]/(16*(zz0 + zz1 + zz2)**(3/2)))/(zz0 + zz1 + zz2)**4
    FJ7[1] = -8*xx[1]*(429*xx[2]**7/(16*(zz0 + zz1 + zz2)**(7/2)) - 693*xx[2]**5/(16*(zz0 + zz1 + zz2)**(5/2)) + 315*xx[2]**3/(16*(zz0 + zz1+ zz2)**(3/2)) - 35*xx[2]/(16*math.sqrt(zz0 + zz1 + zz2)))/(zz0 + zz1 + zz2)**5 + (-3003*xx[1]*xx[2]**7/(16*(zz0 + zz1 + zz2)**(9/2)) + 3465*xx[1]*xx[2]**5/(16*(zz0 + zz1 + zz2)**(7/2)) - 945*xx[1]*xx[2]**3/(16*(zz0 + zz1 + zz2)**(5/2)) + 35*xx[1]*xx[2]/(16*(zz0 + zz1 + zz2)**(3/2)))/(zz0 + zz1 + zz2)**4
    FJ7[2] = -8*xx[2]*(429*xx[2]**7/(16*(zz0 + zz1 + zz2)**(7/2)) - 693*xx[2]**5/(16*(zz0 + zz1 + zz2)**(5/2)) + 315*xx[2]**3/(16*(zz0 + zz1 + zz2)**(3/2)) - 35*xx[2]/(16*math.sqrt(zz0 +zz1  + zz2)))/(zz0 + zz1 + zz2)**5 + (-3003*xx[2]**8/(16*(zz0 + zz1 + zz2)**(9/2)) + 1617*xx[2]**6/(4*(zz0 + zz1 + zz2)**(7/2)) - 2205*xx[2]**4/(8*(zz0 + zz1 + zz2)**(5/2)) + 245*zz2/(4*(zz0 + zz1 + zz2)**(3/2)) - 35/(16*math.sqrt(zz0 + zz1 + zz2)))/(zz0 + zz1 + zz2)**4
     
    zp = np.zeros(6)
    zp[0:3]=z2
        
    zp[3:6] = R.T@(lado+ FJ2 + FJ3 + J4*FJ4 + J5*FJ5 + J6*FJ6 + J7*FJ7 - (2*Rp@z2+R2p@z1)) 
    
    return zp


tiempo= t
delta = tiempo[-1]
t = np.linspace(0, delta, len(tiempo))

sol_odeint= odeint(satelite,z0,t)

x_od = sol_odeint[:,0]
y_od = sol_odeint[:,1]
z_od = sol_odeint[:,2]
vx_od = sol_odeint[:,3]
vy_od = sol_odeint[:,4]
vz_od = sol_odeint[:,5]



# Se escribe datos del Sentinel 1A/B en el archivo EOF
with open(fname_out,"w") as fout:
    cantidad = len(t)
    fout.write("<?xml version=\"1.0\"?>\n")
    fout.write("<Earth_Explorer_File>\n")
    fout.write("  <Earth_Explorer_Header>\n")
    fout.write("    <Fixed_Header>\n")
    fout.write("      <File_Name>S1A_OPER_AUX_POEORB_OPOD_20200821T121202_V20200731T225942_20200802T005942</File_Name>\n")
    fout.write("      <File_Description>Precise Orbit Ephemerides (POE) Orbit File</File_Description>\n")
    fout.write("      <Notes></Notes>\n")
    fout.write("      <Mission>Sentinel-1B</Mission>\n")
    fout.write("      <File_Class>OPER</File_Class>\n")
    fout.write("      <File_Type>AUX_POEORB</File_Type>\n")
    fout.write("      <Validity_Period>\n")
    fout.write("        <Validity_Start>UTC=2020-07-31T22:59:42</Validity_Start>\n")
    fout.write("        <Validity_Stop>UTC=2020-08-02T00:59:42</Validity_Stop>\n")
    fout.write("      </Validity_Period>\n")
    fout.write("      <File_Version>0001</File_Version>\n")
    fout.write("      <Source>\n")
    fout.write("        <System>OPOD</System>\n")
    fout.write("        <Creator>OPOD</Creator>\n")
    fout.write("        <Creator_Version>0.0</Creator_Version>\n")
    fout.write("        <Creation_Date>UTC=2020-08-21T12:12:02</Creation_Date>\n")
    fout.write("      </Source>\n")
    fout.write("    </Fixed_Header>\n")
    fout.write("    <Variable_Header>\n")
    fout.write("      <Ref_Frame>EARTH_FIXED</Ref_Frame>\n")
    fout.write("      <Time_Reference>UTC</Time_Reference>\n")
    fout.write("    </Variable_Header>\n")
    fout.write("  </Earth_Explorer_Header>\n")
    fout.write("<Data_Block type=\"xml\">\n ")
    fout.write ("<List_of_OSVs count=\"9361\">\n")
    # fecha y datos vector
    for i in range(cantidad):
        fecha1 = datetime(2020,7,31,22,59,42,000000)
        fecha2 = (fecha1 + timedelta(seconds=t[i])).strftime("%Y-%m-%dT%H:%M:%S.%f")
        fout.write(" <OSV>\n")
        fout.write(" <UTC>UTC={fecha}</UTC>\n")
        fout.write(" <Absolute_Orbit>+22823</Absolute_Orbit>\n")
        fout.write(" <X unit=\"m\">{x[i]}</X>\n")
        fout.write(" <Y unit=\"m\">{y[i]}</Y>\n")
        fout.write(" <Z unit=\"m\">{z[i]}</Z>\n")
        fout.write(" <VX unit=\"m/s\">{vx[i]}</VX>\n")
        fout.write(" <VY unit=\"m/s\">{vy[i]}</VY>\n")
        fout.write(" <VZ unit=\"m/s\">{vz[i]}</VZ>\n")
        fout.write("   <Quality>NOMINAL</Quality>\n")
        fout.write("   </OSV>\n")
    fout.write("  </List_of_OSVs>\n        </Data_Block>\n       </Earth_Explorer_File>")
