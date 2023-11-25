# Codigos-transient-hysteresis-TFM
- main.m : Este es el programa que hace las gráficas, lo que hago al principio es cargar directamente las matrices con los resultados de las simulaciones. Los nombres para las matrices que he usado son de la forma 
ci7_p1_i1_500 
ci es la condición inicial (en mi caso tenia 2 IPTG = 10^-7 ci7 o 10^-2 ci2)
p1 o p2 es la proteína 
i1, i2 i3 es el intervalo (tuve que hacer la simulación por intervalos porque a mi ordenador llegó un punto que hacer las simulaciones le llevaba la vida
y sino se quedaba sin memoria, igual a ti no te hace falta hacer las simulaciones así)
500 es el tiempo hasta el que ejecutaba la simulación

- diagrama.m : Es el código que utilizaba para ejecutar el Selansi y dónde almacenaba las matrices en cada simulación
