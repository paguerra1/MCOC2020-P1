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
1. Se graficó la posición (x,y,z) en el tiempo del vector estado del Sentinel 1A/B correspondiente.
  ![Preg1_Posición](https://user-images.githubusercontent.com/69210578/92346643-bbd73480-f0a3-11ea-8369-d08ef5d1a811.png)
2. Se realizó una comparación de las coluciones Odeint y Eurelint para Nsubdivisiones=1:
![Preg2_Deriva_EulerintVSOdeint](https://user-images.githubusercontent.com/69210578/92346654-c2fe4280-f0a3-11ea-83b1-17aaa70540ea.png)
  ![Preg2_Deriva_Eulerint](https://user-images.githubusercontent.com/69210578/92346650-c1cd1580-f0a3-11ea-8815-ccab743cdefd.png)
![Preg2_Deriva_Odeint](https://user-images.githubusercontent.com/69210578/92346657-c5609c80-f0a3-11ea-9d2a-6098c987679e.png) 

  * ¿Cuánto deriva Eulerint de odeint en este caso al final del tiempo?
    *  Se notan las diferencias, Eulerint tarda mas tiempo en realizar el trabajo, derivan en 17mil Km aprox.
  * ¿ Cuánto se demora Odeint y Eulerint respectivamente en producir los resultados?
    * Odeint tarda 0.108 segundos, mientras que Eulerint 0.336 segundos.

3. ¿ Cuantas subdivisiones hay que usar para que la predicción con eulerint al final del tiempo esté en menos de un 1% de error? (Comentar tiempo de ejecución de eulerint)

 ![PREGUNTA3_Nsub5000](https://user-images.githubusercontent.com/69210578/92347230-ceeb0400-f0a5-11ea-912f-cf0e8626445a.png)
 ![Preg3_Deriva_Eulerint_Nsub=5000](https://user-images.githubusercontent.com/69210578/92347896-ba0f7000-f0a7-11ea-8acb-4f849db667e8.png)


   * Se realizaron Nsubdivisiones= 5000 Con esto el error fue menor al 1% No se continuo probando para otras subdivisiones ya que el tiempo de ejecución fue bastante. El tiempo de ejecución de Eulerint para este caso fue de: 2374.32 seg  y de Odeint 0.11 seg.
   * La deriva de Eulerint vs Odeint fue -226 km como se logra apreciar en el gráfico.
 

4. Se implementaron las correciones J2 Y J3. 

  ![Preg4_Posición_correción J2 y J3](https://user-images.githubusercontent.com/69210578/92346782-356f2280-f0a4-11ea-9c2e-c9433fdc357a.png)
  * Correción J2:
  
![Preg4J2_Deriva_Odeint](https://user-images.githubusercontent.com/69210578/92347067-2a68c200-f0a5-11ea-9888-83744243e1fd.png)
* Correción J2 y J3:

![Preg4_Deriva_Odeint](https://user-images.githubusercontent.com/69210578/92347072-2c328580-f0a5-11ea-9e04-3725c0c58108.png)

  * ¿Cuanta deriva incurre al agregar las correciones J2 y J3? ¿Cuánto se demora su código en correr?
    * Las mejoras al implementar las correciones J2 y J3 son notorias, resaltando el gráfico de posición donde la orbita real (línea azul) es semejante a la redicha (naranja). Además del gráfico de la distancia entre la posición real y predicha (odeint)  esta es de: 3km
    * El código se demora menos de 1 segundo en correr. Tiempor Odeint: 0.39 seg, Tiempo Eulerint: 0.73seg



# Entrega Final
* Para esta entrega se agregaron 8 términos utilizando el código presentado por el profesor.
* El resultado se puede ver a continuación para los gráficos de la posición,  y la distancia entre la predicha y odeint: 
![Posición_CorrecionFINAL](https://user-images.githubusercontent.com/69210578/92982367-37eec500-f474-11ea-9a4f-cf7348c85b28.png)
![Distancia_odVSeu_CorreciónFinal](https://user-images.githubusercontent.com/69210578/92982369-3c1ae280-f474-11ea-801a-a6f5cf640656.png)
![Final_Deriva_Odeint](https://user-images.githubusercontent.com/69210578/92982377-3fae6980-f474-11ea-92f4-866cfb4e7daa.png)

* Cuando se aplica la correción J2 se llega a una distancia odeint de 5.8 [Km] como se observa en las entregas anteriores. Para esta entrega se efectuó una correción hasta el término J8 (Presentado en el código: EntregaFinal) , obteniendo una distancia de 4.0 [km] La clara diferencia demuestra una ventaja notoria al utilizar estas correcciones. Es un camino tedioso producto de las ecuacionestulizadas, pero finalmente presenta mejoras.
* Además el gráfico de posición no se ve perjudicado de ninguna forma con estas correcciones, al contrario, presenta mayor presición.
* El código demora menos de 3segundos en correr. Pero  ha aumentando los tiempos de las funciones, explicado por el aumento de correciones y ecuaciones utilizadas. Tiempo Odeint: 0.99 seg , Tiempo Eulerint: 3.53

* Finalmente, para realizar correcciones a esto mismo se decide utilizar solo la función "Odeint" (presentada en archivo: Concurso)
