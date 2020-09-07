# MCOC2020-P1: Predicción de Órbitas
# Entrega 1:

![Trayectoria para distintos vientos](https://user-images.githubusercontent.com/69210578/91094449-c9fb6e80-e628-11ea-976b-0585eb8ecd6a.png)
# Entrega 2: 
![Historias_tiempo](https://user-images.githubusercontent.com/69210578/91517269-650b7700-e8bb-11ea-9058-042baab5d821.png)
![Distancia_satéliteVS_Tiempo](https://user-images.githubusercontent.com/69210578/91517319-7b193780-e8bb-11ea-8dc3-34ba77961500.png)
* ¿Cuanto debe valer vt de modo que el satélite efectivamente orbite sin caer de vuelta dentro de la atmosfera (asuma que esta comienza a una altura de 80km)?
  * Asumiendo las dos vueltas que da el satélite, se obtuvo que alrrededor de 24.000 m/s orbitaba sin caer a la atmosfera de tierra

* ¿Como encontró vt?
  * Se observó que para velocidad mayores como 100.000 m/s o menores como 1.000 m/s el satélite caía hacia la atmósfera.
De este modo, se probaron para diversas velocidades, llegando a la conclucion que cercano a los 24.000 m/s, o  un poco inferior a este, el satélite no caía, hasta 23.000 m/s.

* ![Trayectoria del_satelite](https://user-images.githubusercontent.com/69210578/91518124-748bbf80-e8bd-11ea-91fb-d8f4561cd1fe.png)



# Entrega 4:
* El gráfico presenta la solución de la ecuación diferencial del oscilador armónico. Fueron presentadas la solucion analítica, odeint, y eulerint para 1, 10 y 100 subdivisiones.   
* Para la solución real se utilizó la fórmula:  ``` z_real = np.exp((-c/(2*m))*t)*(m*np.cos(om*p*t) + ((1 + om*chi*m)/(om*(p)))*np.sin(om*p*t)) ``` 
* donde p= (1-chi^2)^0.5
* Se observa que para las soluciones de Euler, con subdivisiones de 10 y 100 tienen un comportamiento mas preciso y similar a la solución real. El caso de subdvision=1 pierde precisión, como se muestra en el primer gráfico para ```linspace(0,4.,100) ``` (Gráfico pedido para entrega 4).

![Grafico_entrega4](https://user-images.githubusercontent.com/69210578/91870194-bf526200-ec44-11ea-9457-a38d91a37073.png)

* Luego, a modo de comparación se grafica el caso de 1000 números en vez de solo 100, con ```linspace(0,4.,1000) ```  y se observa que para la solucion de Euler con 1 subdivisión reflejaría mayor precisión que para el caso pedido.

![linspace(0,4 ,1000)](https://user-images.githubusercontent.com/69210578/91870228-c8433380-ec44-11ea-95de-624141c1ef9f.png)


# Entrega 5:
* 1. Se graficó la posición (x,y,z) en el tiempo del vector estado del Sentinel 1A/B correspondiente.
  *
* 2. Se realizó una comparación de las coluciones Odeint y Eurelint para Nsubdivisiones=1:
  *
  * ¿Cuánto deriva Eulerint de odeint en este caso al final del tiempo?
  *
  * ¿ Cuánto se demora Odeint y Eulerint respectivamente en producir los resultados?

* 3. ¿ Cuantas subdivisiones hay que usar para que la predicción con eulerint al final del tiempo esté en menos de un 1% de error? (Comentar tiempo de ejecución de eulerint)
  *
* 4. Se implementaron las correciones J2 Y J3. 
  *
  * ¿Cuanta deriva incurre al agregar las correciones J2 y J3? ¿Cuánto se demora su código en correr?




