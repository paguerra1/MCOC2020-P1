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
* Se observa que para las soluciones de Euler, con subdivisiones de 10 y 100 tienen un comportamiento mas preciso y similar a la solución real. El caso de subdvision=1 pierde presición, como se muestra en el primer gráfico para ```linspace(0,4.,100) ``` (Gráfico pedido para entrega 4).

![linspace(0,4 ,100)](https://user-images.githubusercontent.com/69210578/91866195-3f29fd80-ec40-11ea-81dd-c3301a46ca65.png)

* Luego, a modo de comparación se grafica el caso de 1000 números en vez de solo 100, con ```linspace(0,4.,1000) ```  y se observa que para la solucion de Euler con 1 subdivisión reflejaría mayor precisión que para el caso pedido.

![grafico2linspace(0,4 ,1000)](https://user-images.githubusercontent.com/69210578/91866198-3fc29400-ec40-11ea-88e7-692978a0eeb7.png)

