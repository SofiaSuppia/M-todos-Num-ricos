# Problema 2 — Examen 2025 (Sistemas de ecuaciones lineales 6×6)

Resolución **completa a mano**, paso a paso, con las fórmulas y los datos que
hay que escribir en cada parte.

## Archivos de esta carpeta

| Archivo | Qué es |
|---|---|
| `problema2.c++` | Código: eliminación gaussiana (a) y Gauss-Seidel (b), con el enunciado en los comentarios |
| `problema2` | Ejecutable ya compilado |
| `README - Problema 2.md` | Este documento |

Compilar (si modificás el `.c++`) y ejecutar:

```bash
g++ -O2 -o problema2 problema2.c++ -lm
./problema2
```

> El código no pide datos por teclado: el sistema, el valor inicial y la
> tolerancia están fijos. Se ejecuta y lista.

---

## Enunciado

Dado el sistema:

```
 4·x1 +  x2                     −  x6 =  0
 2·x1 + 10·x2 −  x3        − 2·x5    =  9
        3·x2 +  8·x3        −  x5 − 3·x6 =  7
        − x2         + 4·x4       +  x6 = 20
  x1   − 2·x2                 + 5·x5    = 22
               3·x3  − 2·x4       + 6·x6 = 37
```

**a)** Obtener la solución con **eliminación gaussiana**. ¿Es necesario emplear
pivoteo? Justificar.

**b)** ¿Es posible resolverlo con **Gauss-Seidel**? Justificar. En caso
afirmativo, calcular con error absoluto `10⁻⁴`, con valor inicial `x0j = 1`
(`j = 1..6`). Escribir la solución (10 cifras decimales), el error absoluto
estimado y el número de iteraciones.

---

## El sistema y su solución exacta

La solución del sistema (que saldrá en la parte a) es:

```
x = (1, 2, 3, 4, 5, 6)
```

Se puede comprobar en la primera ecuación: `4·1 + 2 − 6 = 4 + 2 − 6 = 0` ✓.
Sirve para calcular el **error exacto** de Gauss-Seidel después.

---

## Parte a) Eliminación gaussiana (a mano)

### ¿Hace falta pivoteo? NO. Justificación

Se necesita pivoteo solamente si algún elemento de la diagonal es cero (se
dividiría por cero). La diagonal del sistema es:

```
a[1][1]=4   a[2][2]=10   a[3][3]=8   a[4][4]=4   a[5][5]=5   a[6][6]=6
```

Todos son **distintos de cero** → no hay división por cero → la eliminación
simple avanza sin problema.

Además, la matriz es **estrictamente diagonalmente dominante por filas**
(|diagonal| > suma del resto de la fila), por lo que en cada columna el máximo
valor absoluto **ya está en la diagonal**: ni siquiera el **pivoteo parcial**
haría un intercambio.

### Paso 1 — eliminar la columna 1 (pivote `a[1][1] = 4`)

Multiplicadores (solo donde el coeficiente no es cero):

```
m[2][1] = 2/4 = 0.5
m[5][1] = 1/4 = 0.25
```

**Fila 2** ← Fila 2 − 0.5·Fila 1:

```
a[2][2]: 10 − 0.5·1      = 9.5
a[2][3]: −1 − 0.5·0      = −1
a[2][5]: −2 − 0.5·0      = −2
a[2][6]:  0 − 0.5·(−1)   = 0.5
b[2]:     9 − 0.5·0      = 9
```

**Fila 5** ← Fila 5 − 0.25·Fila 1:

```
a[5][2]: −2 − 0.25·1      = −2.25
a[5][5]:  5 − 0.25·0      = 5
a[5][6]:  0 − 0.25·(−1)   = 0.25
b[5]:    22 − 0.25·0      = 22
```

Sistema luego de la columna 1:

```
|  4    1    0    0    0    −1    |   0
|  0    9.5  −1   0    −2    0.5   |   9
|  0    3    8    0    −1    −3    |   7
|  0    −1   0    4    0     1     |  20
|  0    −2.25 0   0    5     0.25  |  22
|  0    0    3    −2   0     6     |  37
```

### Paso 2 — eliminar la columna 2 (pivote `a[2][2] = 9.5`)

Multiplicadores:

```
m[3][2] = 3/9.5     = 0.315789
m[4][2] = −1/9.5    = −0.105263
m[5][2] = −2.25/9.5 = −0.236842
```

Fila 3 ← Fila 3 − 0.315789·Fila 2:

```
a[3][3]: 8 − 0.315789·(−1)                        = 8.3158
a[3][5]: −1 − 0.315789·(−2)                       = −0.3684
a[3][6]: −3 − 0.315789·0.5                        = −3.1579
b[3]:     7 − 0.315789·9                          = 4.1579
```

Fila 4 ← Fila 4 − (−0.105263)·Fila 2:

```
a[4][5]: 0 − (−0.105263)·(−2) = −0.2105
a[4][6]: 1 − (−0.105263)·0.5  = 1.0526
b[4]:    20 − (−0.105263)·9   = 20.9474
```

Fila 5 ← Fila 5 − (−0.236842)·Fila 2:

```
a[5][5]: 5 − (−0.236842)·(−2) = 4.5263
a[5][6]: 0.25 − (−0.236842)·0.5 = 0.3684
b[5]:    22 − (−0.236842)·9   = 24.1316
```

### Pasos 3, 4 y 5 — eliminar columnas 3, 4 y 5

Multiplicadores que usa el código (la matriz ya va quedando triangular):

```
col. 3:  m[4][3] = −0.105263/8.315789 = −0.012658
         m[5][3] = −0.236842/8.315789 = −0.028481
         m[6][3] =  3.000000/8.315789 =  0.360759

col. 4:  m[6][4] = −2/4 = −0.5

col. 5:  m[6][5] = 0.025316/4.515823 = 0.005606
```

### Matriz triangular superior final

```
|  4       1       0       0       0      −1       |   0
|  0       9.5    −1       0      −2       0.5     |   9
|  0       0       8.3158  0      −0.3684  −3.1579  |   4.1579
|  0       0       0       4      −0.2152   1.0127  |  21
|  0       0       0       0       4.5158   0.2785  |  24.25
|  0       0       0       0       0        7.6440  |  45.8641
```

**Determinante = 4·9.5·8.3158·4·4.5158·7.6440 = 43632**
(no hubo intercambio de filas → no se corrige signo). Como es ≠ 0 → solución
única.

### Retrosustitución (de abajo hacia arriba)

```
x6 = 45.8641 / 7.6440 = 6

x5 = (24.25 − 0.2785·x6) / 4.5158
   = (24.25 − 1.6710) / 4.5158 = 22.579 / 4.5158 = 5

x4 = (21 − 1.0127·x6 − (−0.2152)·x5) / 4
   = (21 − 6.0762 + 1.0760) / 4 = 16.000 / 4 = 4

x3 = (4.1579 − (−0.3684)·x5 − (−3.1579)·x6) / 8.3158
   = (4.1579 + 1.8420 + 18.9474) / 8.3158 = 24.947 / 8.3158 = 3

x2 = (9 − (−1)·x3 − (−2)·x5 − 0.5·x6) / 9.5
   = (9 + 3 + 10 − 3) / 9.5 = 19 / 9.5 = 2

x1 = (0 − 1·x2 − (−1)·x6) / 4
   = (0 − 2 + 6) / 4 = 4 / 4 = 1
```

### Respuesta de la parte a

```
x = (1, 2, 3, 4, 5, 6)      →  NO fue necesario pivoteo.
```

En la consola lo confirmás: todos los pasos, la matriz triangular final, el
determinante 43632, la retrosustitución y la verificación `‖A·x − b‖ = 0`.

---

## Parte b) Gauss-Seidel

### ¿Es posible? SÍ. Justificación

**Criterio suficiente de convergencia**: si la matriz es **estrictamente
diagonalmente dominante por filas**, Gauss-Seidel converge para **cualquier**
valor inicial. Se verifica fila por fila (`|a[i][i]| >` suma del resto):

| Fila | `|a[i][i]|` | Σ resto | ¿domina? |
|---|---:|---:|---|
| 1 | \|4\| = 4 | \|1\|+\|-1\| = 2 | sí |
| 2 | \|10\| = 10 | \|2\|+\|-1\|+\|-2\| = 5 | sí |
| 3 | \|8\| = 8 | \|3\|+\|-1\|+\|-3\| = 7 | sí |
| 4 | \|4\| = 4 | \|-1\|+\|1\| = 2 | sí |
| 5 | \|5\| = 5 | \|1\|+\|-2\| = 3 | sí |
| 6 | \|6\| = 6 | \|3\|+\|-2\| = 5 | sí |

→ **SÍ es posible**, converge para cualquier inicial (con `x0j = 1` también).

### Despejes de Gauss-Seidel

De cada fila se despeja la incógnita de la diagonal (el elemento dominante):

```
x1 = ( 0  − x2                  + x6 ) / 4
x2 = ( 9  − 2·x1  + x3        + 2·x5 ) / 10
x3 = ( 7  − 3·x2              + x5 + 3·x6 ) / 8
x4 = ( 20 + x2                 − x6 ) / 4
x5 = ( 22 − x1  + 2·x2                 ) / 5
x6 = ( 37 − 3·x3 + 2·x4                ) / 6
```

**Regla de oro de Gauss-Seidel:** al calcular cada `x[i]` se usan los valores
**recién calculados** en esa misma iteración (los de arriba) y los **viejos**
para los de abajo.

### Iteración 1 — con `x⁽⁰⁾ = (1, 1, 1, 1, 1, 1)`

```
x1 = ( 0  − 1 + 1 ) / 4              = 0.00
x2 = ( 9  − 2·0 + 1 + 2·1 ) / 10     = 12/10 = 1.20
x3 = ( 7  − 3·1.2 + 1 + 3·1 ) / 8    = 7.4/8 = 0.925
x4 = ( 20 + 1.2 − 1 ) / 4            = 20.2/4 = 5.05
x5 = ( 22 − 0 + 2·1.2 ) / 5          = 24.4/5 = 4.88
x6 = ( 37 − 3·0.925 + 2·5.05 ) / 6   = 44.325/6 = 7.3875
```

Nota bien el encadenamiento: `x2` usa el `x1=0` recién calculado; `x3` usa el
`x2=1.2` recién calculado; y así. Resultado:

```
x⁽¹⁾ = ( 0.00, 1.20, 0.925, 5.05, 4.88, 7.3875 )
```

Error absoluto estimado (norma euclídea entre `x⁽¹⁾` y `x⁽⁰⁾`):

```
‖x⁽¹⁾−x⁽⁰⁾‖ = √( (0−1)² + (1.2−1)² + (0.925−1)² + (5.05−1)² +
                  (4.88−1)² + (7.3875−1)² )
            = √( 1 + 0.04 + 0.005625 + 16.4025 + 15.0544 + 40.8002 )
            = √(73.3027)
            = 8.5617
```

### Iteración 2 — se nota cómo está convergiendo

```
x1 = ( 0  − 1.2 + 7.3875 ) / 4               = 1.546875
x2 = ( 9  − 2·1.546875 + 0.925 + 2·4.88 )/10 = 1.659125
x3 = ( 7  − 3·1.659125 + 4.88 + 3·7.3875 )/8 = 3.633141
x4 = ( 20 + 1.659125 − 7.3875 ) / 4          = 3.567906
x5 = ( 22 − 1.546875 + 2·1.659125 ) / 5      = 4.754275
x6 = ( 37 − 3·3.633141 + 2·3.567906 ) / 6    = 5.539398
```

```
x⁽²⁾ = ( 1.546875, 1.659125, 3.633141, 3.567906, 4.754275, 5.539398 )
‖x⁽²⁾−x⁽¹⁾‖ = 3.9453
```

El error **bajó de 8.56 a 3.95** en una sola iteración; los x ya se están
acercando a (1, 2, 3, 4, 5, 6).

### Convergencia (criterio: error absoluto < 1e-4)

El programa corta cuando `‖x⁽ᵏ⁾ − x⁽ᵏ⁻¹⁾‖ < 10⁻⁴`. Evolución del error:

| iter | error |
|---|---:|
| 1 | 8.56170 |
| 2 | 3.94532 |
| 3 | 1.38404 |
| 6 | 0.03994 |
| 9 | 0.00133 |
| 10 | 0.00042 |
| 11 | 0.00014 |
| **12** | **0.000043 < 1e-4** → corta |

### Respuesta de la parte b

```
x1 = 1.0000024733      x4 = 3.9999955257
x2 = 1.9999980382      x5 = 4.9999987206
x3 = 3.0000072039      x6 = 5.9999949066

Error absoluto estimado: 4.3141×10⁻⁵
Iteraciones empleadas:   12
```

La solución exacta es `(1, 2, 3, 4, 5, 6)`: la aproximación difiere en menos de
`7×10⁻⁶` (error exacto máximo en `x3`). El código además imprime la verificación
`‖A·x − b‖ = 7.5×10⁻⁵` ≈ 0 → correcto.

---

## Qué datos va tomando cada parte

| Parte | Datos del enunciado | Qué se escribe como respuesta |
|---|---|---|
| **a** | la matriz `[A\|b]` 6×6 | Solución `(1, 2, 3, 4, 5, 6)`; justificación (diagonal ≠ 0, NO pivoteo) + captura de consola |
| **b** | valor inicial `x0j = 1`, error absoluto `1e-4` | Sí es posible (dominancia); solución `(1.0000024733, 1.9999980382, 3.0000072039, 3.9999955257, 4.9999987206, 5.9999949066)`; error `4.3141e-05`; `12` iteraciones + captura de consola |

En el código están fijos: el sistema `A`, `b`, `x0 = (1,1,1,1,1,1)` y la
tolerancia `1e-4`.

---

## Resumen

| | Eliminación gaussiana | Gauss-Seidel |
|---|---|---|
| Pivoteo necesario | **No** (diagonal 4,10,8,4,5,6 ≠ 0) | — |
| Convergencia | — | Sí (estrict. diagonalmente dominante) |
| Solución | `(1, 2, 3, 4, 5, 6)` exacta | `(1.0000024733, 1.9999980382, 3.0000072039, 3.9999955257, 4.9999987206, 5.9999949066)` |
| Error | 0 (residuo = 0) | estimado `4.3141×10⁻⁵` |
| Iteraciones | 1 (directo) | 12 |

La eliminación gaussiana da la solución **exacta** (det = 43632 ≠ 0). Gauss-
Seidel necesita 12 iteraciones para llegar a un error absoluto menor a `10⁻⁴`,
y converge tan rápido precisamente porque la matriz es estrictamente diagonal
dominante.