# Sistemas de Ecuaciones Lineales

Resolución de sistemas lineales `A·x = b` con tres métodos: **eliminación
gaussiana con pivoteo parcial** (directo), **Jacobi** e **Gauss-Seidel**
(iterativos). Todo con el **paso a paso** impreso en pantalla para después
explicarlo/copiarlo en la carpeta del parcial.

## Archivos de esta carpeta

| Archivo | Qué es |
|---|---|
| `EliminGaussiana.c++` | Eliminación gaussiana con pivoteo parcial, determinante y retrosustitución, todo paso a paso |
| `EliminGaussiana` | Ejecutable ya compilado |
| `JacobiYGaussSeidel.c++` | Métodos iterativos Jacobi y Gauss-Seidel con tabla de iteraciones y verificación por residuo |
| `JacobiYGaussSeidel` | Ejecutable ya compilado |
| `matriz.dat` | Archivo de datos: la matriz aumentada `[ A \| b ]` del sistema |
| `README.md` | Este documento |

---

## 1) Qué hace cada código

### EliminGaussiana.c++ (`./EliminGaussiana`)

Resuelve `A·x = b` **directamente** en 6 pasos:

1. Lee la matriz aumentada de `matriz.dat`.
2. **Pivoteo parcial**: si el pivote es casi cero, intercambia con la fila de
   mayor valor absoluto (más abajo) — por eso la matriz puede **tener un 0 en
   la diagonal** sin drama.
3. **Triangulación**: elimina hacia abajo restando `m[j][i]` veces la fila
   pivote, ¡mostrando cada factor `m[j][i] = a[j][i] / a[i][i]`!
4. **Determinante**: producto de la diagonal. Avisa y corrige si hubo
   intercambios de filas (cada intercambio cambia el signo del det).
5. **Retrosustitución**: resuelve de abajo hacia arriba, mostrando cada cuenta.
6. Imprime la solución.

### JacobiYGaussSeidel.c++ (`./JacobiYGaussSeidel`)

Los dos métodos **iterativos** más conocidos. Primero **verifica la diagonal
dominante** y después:

- **Jacobi**: en cada iteración calcula TODOS los `x[i]` nuevos usando **solo
  los valores de la iteración anterior**.
- **Gauss-Seidel**: al calcular `x[i]` ya usa los valores **nuevos** recién
  calculados (`x[0]..x[i-1]`), por eso suele converger **en menos
  iteraciones**.

Ambos imprimen una **tabla** con el vector completo en cada iteración y el
**error** (norma euclídea entre iteración e iteración), y al final calculan el
**residuo** `‖A·x − b‖` para verificar el resultado.

> **IMPORTANTE:** los métodos iterativos necesitan que en la diagonal no haya
> ceros y conviene que la matriz sea **diagonalmente dominante**. Si el sistema
> original no lo es, **reordená las filas** en `matriz.dat` (ver §3).

---

## 2) Cómo se ejecuta

### Compilar (una vez, si modificaste el .c++)

```bash
g++ -O2 -o EliminGaussiana EliminGaussiana.c++ -lm
g++ -O2 -o JacobiYGaussSeidel JacobiYGaussSeidel.c++ -lm
```

### Ejecutar

```bash
./EliminGaussiana
./JacobiYGaussSeidel
```

`EliminGaussiana` no pide nada: lee `matriz.dat` y resuelve. `JacobiYGaussSeidel`
muestra un menú:

```
--- MENU ---
1- Jacobi
2- Gauss-Seidel
3- Salir
------------
OPCION:      <- elegís 1 o 2
```

Después te pide la **tolerancia** (ej. `0.0001` o `0.000001`). El programa corta
cuando el error entre iteraciones es menor a esa tolerancia.

> Ambos leen `matriz.dat` desde la carpeta donde estás parado. Si se ejecutan
> desde otra carpeta, copiá `matriz.dat` al lado.

---

## 3) Formato de `matriz.dat`

Una fila por ecuación, **sin encabezado**. Cada fila: los `n` coeficientes y el
término independiente:

```
a11 a12 ... a1n  b1
a21 a22 ... a2n  b2
  ...
an1 an2 ... ann  bn
```

Sistema actual (el del parcial 3):

```
0 2 5 -14
1 -4 2 -19
5 1 -2 16
```

Es decir:

```
0·x1 + 2·x2 + 5·x3 = −14
1·x1 − 4·x2 + 2·x3 = −19
5·x1 + 1·x2 − 2·x3 = 16
```

**OJO con el orden de las filas:**

- `EliminGaussiana` **tolera el 0** en la diagonal (hace pivoteo solo).
- `JacobiYGaussSeidel` **corta si hay un 0** en la diagonal. Para este sistema
  hay que **reordenar** las filas dejando una matriz diagonal dominante:

```
5 1 -2 16
1 -4 2 -19
0 2 5 -14
```

(o sea: fila 3, fila 2, fila 1 del sistema original). Con ese orden la
convergencia está asegurada.

---

## 4) Teoría en dos párrafos

**Eliminación gaussiana.** Se convierte la matriz en **triangular superior**
haciendo ceros debajo de la diagonal con `Fila_j ← Fila_j − m·Fila_i`, donde
`m = a[j][i] / a[i][i]`. Si el pivote es 0 (o muy chico), se **pivotea**:
intercambiar filas. Después, la solución sale de atrás hacia adelante:

```
x[n-1] = b[n-1] / a[n-1][n-1]
x[i]   = ( b[i] − Σ a[i][j]·x[j] ) / a[i][i]      (j = i+1 .. n-1)
```

El **determinante** de la matriz triangular es el producto de la diagonal; cada
intercambio de filas le cambia el signo.

**Métodos iterativos.** Partiendo de un vector inicial (el código usa
`x = (0,0,0)`), en cada iteración se despeja cada incógnita de su ecuación:

```
Jacobi:       x[k+1][i] = ( b[i] − Σ a[i][j]·x[k][j]       ) / a[i][i]   (j ≠ i)
Gauss-Seidel: x[k+1][i] = ( b[i] − Σ a[i][j]·x[k+1][j]     ) / a[i][i]
                                   − Σ a[i][j]·x[k][j]                    (j>i)
```

Convergen si la matriz es **estrictamente diagonalmente dominante** (por fila,
`|a[i][i]| > Σ|a[i][j]|`) o si el **radio espectral** de la matriz de iteración
es `< 1`. Se corta cuando `‖x[k+1] − x[k]‖ < tolerancia`.

---

## 5) Cálculo a mano — Eliminación gaussiana (parcial 3)

Sistema original (en el orden del enunciado):

```
[ 0   2   5 | −14 ]
[ 1  −4   2 | −19 ]
[ 5   1  −2 |  16 ]
```

### Paso 1 — Columna 1: pivote `a[1][1] = 0` → hay que pivotear

Se busca el mayor valor absoluto en la columna 1: `|5|` (fila 3) `> |1|` (fila 2)
`> |0|` (fila 1). Se **intercambian fila 1 ↔ fila 3**:

```
[ 5   1  −2 |  16 ]
[ 1  −4   2 | −19 ]
[ 0   2   5 | −14 ]
```

### Paso 2 — Eliminar debajo del pivote `5`

```
m[2][1] = 1 / 5 = 0.2
m[3][1] = 0 / 5 = 0
```

Fila 2 ← Fila 2 − 0.2·Fila 1:
`1−0.2·5 = 0  |  −4−0.2·1 = −4.2  |  2−0.2·(−2) = 2.4  |  −19−0.2·16 = −22.2`

Fila 3 ← Fila 3 − 0·Fila 1: no cambia.

```
[ 5    1     −2    |  16    ]
[ 0   −4.2    2.4  | −22.2  ]
[ 0    2      5    | −14    ]
```

### Paso 3 — Columna 2: pivote `a[2][2] = −4.2` → no hace falta pivoteo

```
m[3][2] = 2 / (−4.2) = −0.476190…
```

Fila 3 ← Fila 3 − (−0.476190)·Fila 2:
`0 | 2−(−0.476190)·(−4.2) = 0  |  5−(−0.476190)·2.4 = 6.142857  |  −14−(−0.476190)·(−22.2) = −24.571429`

```
[ 5    1     −2     |  16      ]
[ 0   −4.2    2.4   | −22.2    ]
[ 0    0      6.142857 | −24.571429 ]
```

### Paso 4 — Determinante

Producto de la diagonal: `5 · (−4.2) · 6.142857 = −129`. Hubo **1 intercambio**
de filas (número impar) → se le cambia el signo: **det(A) = 129**. Como es
distinto de cero, el sistema tiene **solución única**.

### Paso 5 — Retrosustitución (de abajo hacia arriba)

```
x3 = −24.571429 / 6.142857 = −4

x2 = ( b2 − a23·x3 ) / a22
   = ( −22.2 − 2.4·(−4) ) / (−4.2)
   = ( −22.2 + 9.6 ) / (−4.2) = ( −12.6 ) / (−4.2) = 3

x1 = ( b1 − a12·x2 − a13·x3 ) / a11
   = ( 16 − 1·3 − (−2)·(−4) ) / 5
   = ( 16 − 3 − 8 ) / 5 = 5 / 5 = 1
```

**Solución:** `x = (1, 3, −4)`.

---

## 6) Cálculo a mano — Gauss-Seidel (parcial 3)

Como `a[1][1] = 0`, el método iterativo **no se puede aplicar directamente**:
primero se **reordenan las filas** (E3, E2, E1) dejando una matriz estrictamente
diagonalmente dominante:

```
[ 5   1  −2 |  16 ]
[ 1  −4   2 | −19 ]
[ 0   2   5 | −14 ]
```

Verificamos dominancia por fila, `|a[i][i]| > Σ|a[i][j]|`:

| Fila | |diag| | Σ|fuera de diag| | ¿domina? |
|---|---:|---:|---|
| 1 | 5 | 1+2 = 3 | sí |
| 2 | 4 | 1+2 = 3 | sí |
| 3 | 5 | 2 | sí |

→ convergencia asegurada. (El radio espectral de la matriz de iteración es
`ρ ≈ 0.29 < 1`, por eso converge rápido.)

### Despejes (formula de Gauss-Seidel)

Fórmulas de iteración (leyendo la fila `i` y despejando `x[i]`):

```
x1 = ( 16 − 1·x2 + 2·x3 ) / 5
x2 = ( −19 − 1·x1 − 2·x3 ) / (−4)
x3 = ( −14 − 0·x1 − 2·x2 ) / 5
```

Partimos de `x⁽⁰⁾ = (0, 0, 0)`.

### Iteración 1 — se usa lo que ya se calculó en esta misma iteración

```
x1⁽¹⁾ = (16 − 1·0 + 2·0) / 5    = 16/5   = 3.2
x2⁽¹⁾ = (−19 − 1·3.2 − 2·0) / (−4) = (−22.2)/(−4) = 5.55
x3⁽¹⁾ = (−14 − 0·3.2 − 2·5.55) / 5 = (−25.1)/5 = −5.02
```

Nota: el `x1 = 3.2` recién calculado se **usa en la misma iteración** para
`x2` (y `x1`, `x2` para `x3`) — eso es lo propio de Gauss-Seidel.

**x⁽¹⁾ = (3.2, 5.55, −5.02)**

Error (norma euclídea entre iteraciones):

```
‖x⁽¹⁾−x⁽⁰⁾‖ = √( 3.2² + 5.55² + (−5.02)² ) = 8.138974…
```

(Como la solución exacta es (1,3,−4), el error exacto de la primera iteración
es `‖x⁽¹⁾−x*‖ = 3.518934…`.)

### Iteración 2

```
x1⁽²⁾ = (16 − 1·5.55 + 2·(−5.02)) / 5 = (16 − 5.55 − 10.04)/5 = 0.41/5 = 0.082
x2⁽²⁾ = (−19 − 1·0.082 − 2·(−5.02)) / (−4) = (−9.042)/(−4) = 2.2605
x3⁽²⁾ = (−14 − 0·0.082 − 2·2.2605) / 5 = (−18.521)/5 = −3.7042
```

**Nota: aquí se ve cómo el error baja rápido** — los valores ya se acercan a
(1, 3, −4):

```
‖x⁽²⁾−x⁽¹⁾‖ = √( (0.082−3.2)² + (2.2605−5.55)² + (−3.7042+5.02)² ) = 4.719541…
```

### Convergencia (tol = 10⁻⁴)

El programa corta en la **iteración 11** (error estimado = `7.84×10⁻⁵ < tol`):

```
  iter |     x[0]          x[1]          x[2]      |      error
  -----------------------------------------------------
    1  |   3.20000000   5.55000000  -5.02000000 | 8.138974136831
    2  |   0.08200000   2.26050000  -3.70420000 | 4.719540643961
   ...
   11  |   1.00001332   3.00001073  -4.00000429 | 0.000078428874

--- SOLUCION (Gauss-Seidel) ---
x[0] = 1.000013317560
x[1] = 3.000010728034
x[2] = -4.000004291214
Error: 0.000078428874
Iteraciones: 11
Residuo ||A*x - b|| = 0.000093999974
```

El residuo `‖A·x − b‖ ≈ 9.4×10⁻⁵` confirma que el resultado es correcto (debería
dar ~0 si `x` fuera exacto).

---

## 7) Verificación del resultado

Dos maneras independientes de checkear:

1. **Reemplazar en el sistema:** con `x = (1, 3, −4)`:
   ```
   E1:  2·3 + 5·(−4) = 6 − 20 = −14  ✓
   E2:  1 − 4·3 + 2·(−4) = 1 − 12 − 8 = −19  ✓
   E3:  5·1 + 3 − 2·(−4) = 5 + 3 + 8 = 16  ✓
   ```
2. **Residuo:** `‖A·x − b‖` chico (el código lo calcula al final de cada método
   iterativo).

---

## 8) Relación con el parcial

El sistema de §3, §5 y §6 es el del **Problema 3 del Parcial 1 (2024)**.
Además de estos códigos, tenés resuelto todo el problema (a, b y c) en:

```
Practica/Examen 2024/Parcial 1/problema3.c++
```

que hace **eliminación con pivoteo** (parte a), el **análisis de convergencia**
(parte b) y Gauss-Seidel con error estimado **y** error exacto (parte c).

---

## 9) Limitaciones

- **EliminGaussiana:** si la matriz es singular (det = 0) avisa y corta. Para
  matrices casi singulares el pivoteo reduce (no elimina) los errores de
  redondeo.
- **Jacobi/Gauss-Seidel:** necesitan diagonal sin ceros (el código corta si
  encuentra un 0). Si la matriz **no** es diagonalmente dominante solo **avisa**;
  el método puede converger o no — lo seguro es **reordenar filas** hasta tener
  dominancia.
- **Máximo 50 ecuaciones** (`MAX_SIZE`).
- El criterio de corte es `‖x[k+1] − x[k]‖ < tolerancia` (norma euclídea).