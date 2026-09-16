# Problema 3 — Examen 2025 (Mínimos cuadrados: f(x) = a·sen x + b·cos x)

Resolución **completa a mano**, paso a paso, con las fórmulas y los datos que
hay que escribir en cada parte.

## Archivos de esta carpeta

| Archivo | Qué es |
|---|---|
| `problema3.c++` | Código: ecuaciones normales 2×2, Gauss (`gauss.h`), St, Sr y r. Enunciado en los comentarios |
| `problema3` | Ejecutable ya compilado |
| `gauss.h` | Eliminación gaussiana con pivoteo parcial (igual a `4. Ajuste de curvas/`) |
| `data.txt` | Datos: primera línea el nº de puntos (7) y luego pares `xi yi` |
| `README - Problema 3.md` | Este documento |

Compilar (si modificás el `.c++`) y ejecutar:

```bash
g++ -O2 -o problema3 problema3.c++ -lm
./problema3
```

> El programa lee `data.txt` (de la carpeta donde estás parado) y no pide nada
> por teclado.

---

## Enunciado

Los datos de la tabla pueden ajustarse por una función de la forma

```
f(x) = a·sen(x) + b·cos(x)
```

El método de mínimos cuadrados minimiza

```
ε = Σ(k=1..N) ( yk − f(xk) )²
```

**a)** Obtener el sistema de ecuaciones lineales que deben satisfacer `a` y
`b` para que `ε` sea mínimo (derivando respecto de `a` y `b`). Escribir la
**matriz de coeficientes** y el **vector independiente**.

**b)** Resolver el sistema por **eliminación gaussiana**. Escribir el **vector
solución**, `St`, `Sr` y el **coeficiente de correlación** `r`.

Datos de la tabla:

| x | 1.0 | 1.2 | 1.5 | 2.0 | 2.6 | 2.9 | 3.0 |
|---|---|---|---|---|---|---|---|
| y | −0.236 | 0.209 | 0.853 | 1.746 | 2.231 | 2.163 | 2.110 |

---

## Parte a) Sistema de ecuaciones normales (derivando a mano)

Se quiere minimizar la suma de residuos al cuadrado:

```
ε(a, b) = Σ ( yk − a·sen(xk) − b·cos(xk) )²
```

**Regla:** para minimizar, se deriva respecto de cada parámetro y se iguala a
cero.

### Derivada respecto de `a`

Por la regla de la cadena, la derivada de `( yk − a·sen − b·cos )²` respecto
de `a` es `2·( yk − a·sen − b·cos )·(−sen)`:

```
∂ε/∂a = Σ 2·( yk − a·sen(xk) − b·cos(xk) )·(−sen(xk)) = 0

Divido 2 y distribuyo:
 Σ yk·sen(xk)  −  a·Σ sen²(xk)  −  b·Σ sen(xk)·cos(xk)  = 0

→  a·Σsen²  +  b·Σsen·cos  =  Σ y·sen          (ecuación 1)
```

### Derivada respecto de `b`

Idéntica, pero multiplicando por `(−cos)`:

```
∂ε/∂b = Σ 2·( yk − a·sen(xk) − b·cos(xk) )·(−cos(xk)) = 0

→  a·Σsen·cos  +  b·Σcos²  =  Σ y·cos          (ecuación 2)
```

### El sistema 2×2

```
┌  Σsen²      Σsen·cos  ┐ ┌a┐   ┌  Σy·sen  ┐
│                       │ │ │ = │          │
└  Σsen·cos   Σcos²     ┘ └b┘   └  Σy·cos  ┘
```

Solo falta calcular las cinco sumatorias.

---

## Tabla de sumatorias (cálculo a mano)

Se arma una columna por cada término que aparece en el sistema. Por punto se
calcula `sen(xi)`, `cos(xi)` y los productos:

| i | xi | yi | sen(xi) | cos(xi) | sen²(xi) | sen·cos | cos²(xi) | y·sen | y·cos |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1.0 | −0.236 | 0.841471 | 0.540302 | 0.708073 | 0.454649 | 0.291927 | −0.198587 | −0.127511 |
| 2 | 1.2 | 0.209 | 0.932039 | 0.362358 | 0.868697 | 0.337732 | 0.131303 | 0.194796 | 0.075733 |
| 3 | 1.5 | 0.853 | 0.997495 | 0.070737 | 0.994996 | 0.070560 | 0.005004 | 0.850863 | 0.060339 |
| 4 | 2.0 | 1.746 | 0.909297 | −0.416147 | 0.826822 | −0.378401 | 0.173178 | 1.587633 | −0.726592 |
| 5 | 2.6 | 2.231 | 0.515501 | −0.856889 | 0.265742 | −0.441727 | 0.734258 | 1.150084 | −1.911719 |
| 6 | 2.9 | 2.163 | 0.239249 | −0.970958 | 0.057240 | −0.232301 | 0.942760 | 0.517496 | −2.100183 |
| 7 | 3.0 | 2.110 | 0.141120 | −0.989992 | 0.019915 | −0.139708 | 0.980085 | 0.297763 | −2.088884 |
| **Σ** | | | | | **3.741485** | **−0.329197** | **3.258515** | **4.400049** | **−6.818818** |

Ejemplo (punto 1, x = 1.0):

```
sen(1.0) = 0.841471      cos(1.0) = 0.540302
sen²  = (0.841471)²         = 0.708073
sen·cos = 0.841471·0.540302 = 0.454649
cos²  = (0.540302)²         = 0.291927
y·sen = (−0.236)·0.841471   = −0.198587
y·cos = (−0.236)·0.540302   = −0.127511
```

> Se puede verificar que `sen² + cos² = 1` en cada fila (identidad
> pitagórica). El ángulo está en **radianes** (la calculadora/el código en
> radianes: `sin(1.0 rad) = 0.841471`, no 0.01745).

### Respuesta de la parte a (matriz y vector)

```
┌  3.741485   −0.329197 ┐ ┌a┐   ┌  4.400049 ┐
│                       │ │ │ = │           │
└ −0.329197    3.258515 ┘ └b┘   └ −6.818818 ┘
```

---

## Parte b) Resolución por eliminación gaussiana (a mano)

### Matriz aumentada

```
|  3.741485   −0.329197 |  4.400049 |
| −0.329197    3.258515 | −6.818818 |
```

### Paso 1 — eliminar debajo del pivote `a[1][1] = 3.741485`

```
m[2][1] = −0.329197 / 3.741485 = −0.087985
```

Fila 2 ← Fila 2 − (−0.087985)·Fila 1:

```
a[2][2]: 3.258515 − (−0.087985)·(−0.329197) = 3.258515 − 0.028965 = 3.229550
b[2]:    −6.818818 − (−0.087985)·4.400049  = −6.818818 + 0.387142 = −6.431676
```

```
|  3.741485   −0.329197 |  4.400049 |
|  0.000000    3.229550 | −6.431676 |
```

### Determinante

```
det = 3.741485·3.229550 = 12.083314   (≠ 0 → solución única)
```

### Retrosustitución

```
b = b[2] / a[2][2] = −6.431676 / 3.229550 = −1.991508

a = ( b[1] − a[1][2]·b_coef ) / a[1][1]
  = ( 4.400049 − (−0.329197)·(−1.991508) ) / 3.741485
  = ( 4.400049 − 0.655607 ) / 3.741485
  = 3.744442 / 3.741485
  = 1.000792
```

### Respuesta de la parte b (vector solución)

```
a = 1.0007923964
b = −1.9915083520

f(x) = 1.0007923964·sen(x) − 1.9915083520·cos(x)
```

> Da casi los coeficientes "enteros agradables" `a ≈ 1` y `b ≈ −2`, o sea
> `f(x) ≈ sen(x) − 2·cos(x)`.

---

## St, Sr y coeficiente de correlación (a mano)

El modelo ya está ajustado. Ahora se predice en cada punto, se calcula el
residuo `ei = yi − f(xi)` y se arman las sumas:

| i | xi | yi | f(xi) | ei = yi − f(xi) | ei² |
|---|---:|---:|---:|---:|---:|
| 1 | 1.0 | −0.236 | −0.233879 | −0.002121 | 0.0000044995 |
| 2 | 1.2 | 0.209 | 0.211139 | −0.002139 | 0.0000045759 |
| 3 | 1.5 | 0.853 | 0.857412 | −0.004412 | 0.0000194628 |
| 4 | 2.0 | 1.746 | 1.738778 | +0.007222 | 0.0000521594 |
| 5 | 2.6 | 2.231 | 2.222411 | +0.008589 | 0.0000737716 |
| 6 | 2.9 | 2.163 | 2.173110 | −0.010110 | 0.0001022162 |
| 7 | 3.0 | 2.110 | 2.112810 | −0.002810 | 0.0000078970 |

Suma de los residuos al cuadrado:

```
Sr = 0.0000044995 + 0.0000045759 + 0.0000194628 + 0.0000521594 +
     0.0000737716 + 0.0001022162 + 0.0000078970
   = 0.0002645825
```

Promedio de las `y`:

```
ȳ = ( −0.236 + 0.209 + 0.853 + 1.746 + 2.231 + 2.163 + 2.110 ) / 7
  = 9.076 / 7
  = 1.2965714286
```

Suma total de cuadrados `St = Σ (yi − ȳ)²`:

```
St = (−1.532571)² + (−1.087571)² + (−0.443571)² + (0.449429)² +
     (0.934429)² + (0.866429)² + (0.813429)²
   = 6.2158497143
```

Coeficiente de correlación:

```
r = √( (St − Sr) / St )
  = √( (6.2158497143 − 0.0002645825) / 6.2158497143 )
  = √( 0.9999574 )
  = 0.9999787169
```

### Respuesta final de la parte b

```
a = 1.0007923964     b = −1.9915083520

St = 6.2158497143
Sr = 0.0002645825
r  = 0.9999787169  → el ajuste es EXCELENTE (r muy cercano a 1)
```

`r` tan cercano a 1 me ayuda a ver cuán bien el modelo seno–coseno explica los
datos: solo queda un `Sr` del orden de 10⁻⁴ (los residuos son menores a 0.02).

---

## Qué datos va tomando cada parte

| Parte | Datos que necesita | Qué se escribe como respuesta |
|---|---|---|
| **a** | la tabla `(xi, yi)` y el modelo `f = a·sen + b·cos` | Derivación (∂ε/∂a = 0, ∂ε/∂b = 0) y el sistema 2×2 con `A = [[3.741485, −0.329197], [−0.329197, 3.258515]]`, `b = [4.400049, −6.818818]` |
| **b** | el sistema de ecuaciones normales | `a = 1.0007923964`, `b = −1.9915083520`; `St = 6.2158497143`, `Sr = 0.0002645825`, `r = 0.9999787169` + captura de consola |

En el código están fijos: `data.txt` (7 puntos) y la construcción de las
ecuaciones normales. La eliminación gaussiana la hace `gauss.h` (con pivoteo
parcial; acá no haría falta porque la matriz es bien condicionada y de
diagonal dominada por `Σsen²`, `Σcos²`).

---

## Resumen

| Concepto | Valor |
|---|---|
| Ecuaciones normales | `a·Σsen² + b·Σsen·cos = Σy·sen` ; `a·Σsen·cos + b·Σcos² = Σy·cos` |
| Matriz | `[[3.741485, −0.329197], [−0.329197, 3.258515]]` |
| Vector independiente | `[4.400049, −6.818818]` |
| Solución | `a = 1.0007923964`, `b = −1.9915083520` |
| Modelo | `f(x) = 1.0007923964·sen(x) − 1.9915083520·cos(x)` |
| St | `6.2158497143` |
| Sr | `0.0002645825` |
| r | `0.9999787169` (excelente) |