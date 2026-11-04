# Problema 1 — Examen 2025 (Raíces de ecuaciones no lineales)

Resolución **completa a mano** del problema, paso a paso, con las fórmulas y
los datos que hay que escribir en cada parte.

## Archivos de esta carpeta

| Archivo | Qué es |
|---|---|
| `problema1.c++` | Código: raíz analítica (a), bisección (b), Newton-Raphson 1ra y 4ta iteración (c), con el enunciado en los comentarios |
| `problema1` | Ejecutable ya compilado |
| `README - Problema 1.md` | Este documento |

Para compilar (si modificás el `.c++`):

```bash
g++ -O2 -o problema1 problema1.c++ -lm
```

Para ejecutar:

```bash
./problema1
```

> El código no pide datos por teclado: el intervalo, la tolerancia y `x0`
> están fijos, se ejecuta y lista.

---

## Enunciado

Dada la función

```
f(x) = (x + 1)/(x + 4) − 0.25·x
```

**a)** Determinar la raíz positiva de `f` analíticamente.

**b)** Obtener la raíz por **bisección** con intervalo inicial `[1.75, 3.2]` y
error porcentual de `0.01%`, usando el criterio

```
ε = |x_nuevo − x_viejo| / ( 1/2·|x_nuevo + x_viejo| ) · 100      (#)
```

donde `x_viejo` en la primera iteración es `a = 3.2`. Escribir la raíz (10
cifras decimales), el error estimado, el error porcentual exacto y el número de
iteraciones.

**c1)** A mano: primera iteración de Newton-Raphson con `x0 = 1.0`. Escribir la
raíz obtenida, el error exacto y el error estimado (criterio #).

**c2)** Con el código: raíz en la 4º iteración (10 cifras decimales) y error
estimado (#).

---

## Parte a) Raíz analítica (a mano)

Se resuelve `f(x) = 0`. La incógnita aparece en un denominador, así que
primero se despeja multiplicando todo por `(x + 4)` (recordando que debe ser
`x ≠ −4`):

```
(x + 1)/(x + 4) − 0.25·x = 0

x + 1 − 0.25·x·(x + 4) = 0        ← multipliqué todo por (x + 4)
x + 1 − 0.25·x² − x    = 0        ← apliqué distributiva 0.25·x·(x+4)
1 − 0.25·x²            = 0        ← x − x = 0
0.25·x²                = 1
x²                     = 4
x                      = ±√4 = ±2
```

**La raíz positiva es `x = 2`** (la negativa es `−2`, que no se pide).

> Este valor `x = 2` es el que se usa después para calcular el **error
> porcentual exacto**: `|x_k − 2|/2·100`.

---

## Parte b) Bisección

### Datos de entrada

| Dato | Valor |
|---|---|
| Intervalo inicial | `[a, b] = [1.75, 3.2]` |
| Tolerancia (error porcentual) | `0.01%` |
| `x_viejo` de la 1ª iteración | `3.2` *(como dice el enunciado)* |
| Criterio de error | `ε = |x_nuevo − x_viejo| / (½·|x_nuevo + x_viejo|) · 100` |
| Raíz exacta (para el error exacto) | `2` |

### Fórmula del método

En cada paso se toma el **punto medio** del intervalo:

```
c = (a + b) / 2
```

- Si `f(a)·f(c) > 0`: la raíz está a la derecha de `c` → `a = c`.
- Si `f(a)·f(c) < 0`: la raíz está a la izquierda de `c` → `b = c`.

Se sigue hasta que el error (#) sea `ε < 0.01`.

### Verificación previa (hay raíz en el intervalo)

Como `f(a)·f(b) < 0`, hay al menos una raíz en `[1.75, 3.2]`:

```
f(1.75) = (2.75)/(5.75) − 0.25·1.75 = 0.478261 − 0.4375   = 0.040761  > 0
f(3.20) = (4.20)/(7.20) − 0.25·3.20 = 0.583333 − 0.8      = −0.216667 < 0
f(1.75)·f(3.2) = (+)·(−) = (−)   →   hay raíz en [1.75, 3.2]   ✓
```

### Iteración 1

```
c₁ = (a + b)/2 = (1.75 + 3.2)/2 = 2.475
f(2.475) = (3.475)/(6.475) − 0.25·2.475 = 0.536680 − 0.618750 = −0.082070

Como f(a)·f(c₁) = (+)(−) < 0  →  la raíz está en [1.75, 2.475]  →  b = 2.475
```

Error estimado (#), con `x_viejo = 3.2`:

```
ε₁ = |2.475 − 3.2| / ( ½·|2.475 + 3.2| ) · 100
   = 0.725 / ( ½·5.675 ) · 100
   = 0.725 / 2.8375 · 100
   = 25.5507%
```

### Iteración 2

```
c₂ = (1.75 + 2.475)/2 = 2.1125
f(2.1125) = (3.1125)/(6.1125) − 0.25·2.1125 = 0.509202 − 0.528125 = −0.018923

f(a)·f(c₂) < 0  →  b = 2.1125   (intervalo [1.75, 2.1125])
ε₂ = |2.1125 − 2.475| / ( ½·|2.1125 + 2.475| ) · 100
   = 0.3625 / 2.29375 · 100
   = 15.8038%
```

### Iteración 3

```
c₃ = (1.75 + 2.1125)/2 = 1.93125
f(1.93125) = (2.93125)/(5.93125) − 0.25·1.93125 = 0.494203 − 0.482813 = 0.011392

f(a)·f(c₃) = (+)(+) > 0  →  a = 1.93125   (intervalo [1.93125, 2.1125])
ε₃ = |1.93125 − 2.1125| / ( ½·|1.93125 + 2.1125| ) · 100
   = 0.18125 / 2.021875 · 100
   = 8.9644%
```

### Iteración 4

```
c₄ = (1.93125 + 2.1125)/2 = 2.021875
f(2.021875) = −0.003653 < 0   →  b = 2.021875
ε₄ = |2.021875 − 1.93125| / ( ½·|2.021875 + 1.93125| ) · 100
   = 0.090625 / 1.9765625 · 100
   = 4.5850%
```

### ...y así hasta ε < 0.01%

Cada paso **divide el intervalo a la mitad** y el error (#) también baja a la
mitad. Después de la iteración 12 el error es `0.0177%` (todavía > 0.01) y en
la **iteración 13** ya es `0.0088% < 0.01%`. Ese es el corte.

### Resultados de la parte b

| Resultado | Valor |
|---|---|
| Raíz (10 cifras decimales) | `2.0001037598` |
| Error estimado (criterio #) | `0.008849 %` |
| Error porcentual exacto (`\|2.0001037598 − 2\|/2·100`) | `0.005188 %` |
| Iteraciones empleadas | `13` |

La tabla completa de las 13 iteraciones (a, b, f(a), c, f(c), error abs, error
# y error exacto) la imprime el programa en la salida de la consola.

---

## Parte c1) Newton-Raphson a mano (1ª iteración)

### La derivada, paso a paso

La función es una **diferencia** de dos términos:

```
f(x) = (x + 1)/(x + 4)     −     0.25·x
       └─────┬─────┘              └───┬───┘
        cociente u/v                polinomio
```

**Término 1 — regla del cociente** `(u/v)' = (u'·v − u·v')/v²`:

```
u = x + 1   →  u' = 1
v = x + 4   →  v' = 1

(u/v)' = ( 1·(x + 4) − (x + 1)·1 ) / (x + 4)²
       = ( x + 4 − x − 1 ) / (x + 4)²
       = 3 / (x + 4)²
```

**Término 2 — regla de la potencia**: la derivada de `0.25·x` es `0.25`, así
que la de `−0.25·x` es `−0.25`.

**Derivada total:**

```
f'(x) = 3 / (x + 4)²  −  0.25
```

> Es la derivada que se usa "manualmente". El código `problema1.c++` en cambio
> usa la derivada **numérica** `(f(x+h)−f(x))/h` (igual que tu
> `newtonRapshon.cpp`), que sirve para cualquier `f`.

### Evaluación en x₀ = 1

```
f(1)  = (1 + 1)/(1 + 4) − 0.25·1 = 2/5 − 0.25 = 0.40 − 0.25 = 0.15
f'(1) = 3/(1 + 4)² − 0.25 = 3/25 − 0.25 = 0.12 − 0.25 = −0.13
```

### Fórmula de Newton-Raphson y primera iteración

```
x_{n+1} = x_n − f(x_n) / f'(x_n)

x₁ = x₀ − f(x₀)/f'(x₀)
   = 1 − 0.15/(−0.13)
   = 1 + 1.153846
   = 2.153846
```

### Errores de la primera iteración

```
Error estimado (criterio #):  ε = |x₁ − x₀| / ( ½·|x₁ + x₀| ) · 100
                                = |2.153846 − 1| / ( ½·|2.153846 + 1| ) · 100
                                = 1.153846 / 1.576923 · 100
                                = 73.1707%

Error exacto (absoluto):        |x₁ − 2| = |2.153846 − 2| = 0.153846
Error porcentual exacto:        0.153846 / 2 · 100 = 7.6923%
```

**Respuesta c1: `x₁ = 2.1538461538`, error exacto `0.1538461538`, error
estimado (#) `73.170732%`.**

> Con la derivada numérica del código (h = 0.001) el resultado difiere apenas en
> las últimas cifras (`x₁ = 2.1536332180`), porque `f'(1) ≈ −0.130024` en vez de
> exactamente −0.13. En el papel va la cuenta con la derivada analítica.

---

## Parte c2) Newton-Raphson con el código (4ª iteración)

Iterando la fórmula desde `x₀ = 1` (el código usa la derivada numérica):

```
 iter |   x_nuevo    |   f(x_nuevo)    f'(x_nuevo)  |  error est. (#)
------+--------------+------------------+----------------------------
  0   |  1.0000000000 |   0.150000   −0.130024   |       ----
  1   |  2.1536332180 |  −0.025925   −0.170789   |   73.16216809
  2   |  2.0018363910 |  −0.000306   −0.166732   |    7.30588074
  3   |  2.0000004336 |   0.000000   −0.166681   |    0.09175574
  4   |  2.0000000000 |   0.000000   −0.166681   |    0.00002168
```

**Respuesta c2: raíz `2.0000000000` (4ª iteración), error estimado (#)
`0.000022 %`.**

Se observa la **convergencia cuadrática** de Newton: en 4 pasos el error (#)
pasa de `73%` a `0.00002%`.

> Nota: si se partiera de la `x₁` exacta del papel (2.153846) la convergencia
> sería casi idéntica (la diferencia entre ambas `x₁` es del orden de 10⁻⁴ y no
> cambia el resultado de la 4ª iteración).

---

## Qué datos va tomando cada parte

| Parte | Datos que necesita | Qué se entrega/escribe |
|---|---|---|
| **a** | Ninguno (álgebra a mano) | La deducción y `x = 2` |
| **b** | `a = 1.75`, `b = 3.2`, tolerancia `0.01%`, `x_viejo` inicial `= 3.2` | Raíz `2.0001037598`, `ε(#) = 0.008849%`, error exacto `0.005188%`, `13` iteraciones + captura de consola |
| **c1** | `x₀ = 1.0`, derivada analítica `f'(x) = 3/(x+4)² − 0.25` | `x₁ = 2.1538461538`, error exacto `0.1538461538`, `ε(#) = 73.170732%` |
| **c2** | `x₀ = 1.0` (derivada numérica en el código) | Raíz `2.0000000000`, `ε(#) = 0.000022%` en la 4ª iteración + captura de consola |

---

## Resumen de resultados

| | Bisección | Newton (1ª iter.) | Newton (4ª iter.) |
|---|---|---|---|
| Raíz | `2.0001037598` | `2.1538461538` | `2.0000000000` |
| Error estimado (#) | `0.008849 %` | `73.170732 %` | `0.000022 %` |
| Error exacto | `0.005188 %` | `0.1538461538` | `≈ 0` |
| Iteraciones | `13` | 1 | 4 |

Raíz exacta de `f`: **`x = 2`**. Bisección: seguro pero lento (necesita 13
iteraciones). Newton: mucho más rápido (converge en 4 soles, si el valor
inicial está cerca).