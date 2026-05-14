# Algoritmo-Genetico
Proyecto Final 


***Mejoras del algoritmo***

**1. Selección por torneo**

¿Cómo funciona?
Se seleccionan varios individuos aleatorios.
Se comparan sus fitness.
Gana el mejor.
Ese individuo pasa a reproducirse.

Ejemplo

| Individuo | Fitness |
| --------- | ------- |
| A         | 120     |
| B         | 90      |
| C         | 150     |

Si el torneo toma:

A
B
C

Gana C porque tiene mayor fitness.

¿Por qué mejora?

-Mantiene competencia
-Conserva diversidad
-Evita que todos los individuos sean iguales muy rápido

**2. Elitismo**

¿Qué significa?

Los mejores individuos:

-NO se eliminan
-pasan directamente a la siguiente generación

¿Por qué ayuda?

Sin elitismo:

-Una mutación mala puede destruir la mejor solución

Con elitismo:

-Siempre se conservan los mejores individuos encontrados

**3. Mutación adaptativa**

¿Qué hace la mutación adaptativa?
Si el fitness no mejora durante varias generaciones:
if historial[-1] == historial[-2] == historial[-3]:

el algoritmo:
aumenta temporalmente la mutación

¿Qué logra?

-Introduce nuevas soluciones
-Escapa de mínimos locales
-Encuentra mejores alineamientos


