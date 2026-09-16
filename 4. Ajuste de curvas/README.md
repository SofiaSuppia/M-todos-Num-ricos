# Ajuste de curvas por Mínimos Cuadrados (Regresión)

## Archivos de esta carpeta

| Archivo | Qué es |
|---|---|
| `regresionMinimosCuadrados.cpp` | Código fuente: ajuste polinomial por mínimos cuadrados (grado configurable) |
| `regresionMinimosCuadrados` | Ejecutable ya compilado |
| `gauss.h` | Eliminación gaussiana con pivoteo parcial (la usa el código) |
| `data.txt` | Archivo de datos: una pareja `xi yi` por línea |
| `README.md` | Este documento |

---

## 1) Qué hace el código

Dados n puntos `(xi, yi)`, encuentra el **polinomio de grado `k`** que mejor se
ajusta en el sentido de mínimos cuadrados:

```
P(x) = a0 + a1·x + a2·x² + ... + ak·x^k
```

- **Grado 1** → regresión lineal: `P(x) = a0 + a1·x`
- **Grado 2** → regresión cuadrática: `P(x) = a0 + a1·x + a2·x²`
- Grado mayor → polinomio de mayor grado.

Además calcula `Sr` (suma de cuadrados de los residuos), `St` (suma total de
cuadrados) y el **coeficiente de correlación `r`**.

---

## 2) Cómo se ejecuta

### Compilar (una vez, si modificaste el .c++)

```bash
g++ -O2 -o regresionMinimosCuadrados regresionMinimosCuadrados.cpp -lm
```

### Ejecutar

```bash
./regresionMinimosCuadrados
```

Apenas arranca te pide:

```
Ingrese el grado del polinomio de ajuste (grado >= 1):
```

Escribís `1` (lineal), `2` (cuadrática) o el grado que quieras.

> El programa lee `data.txt` desde la carpeta donde estás parado (ruta relativa).
> Si se ejecuta desde otra carpeta, copiá `data.txt` al lado.

### Formato de `data.txt`

Una línea por punto, separados por espacio, **sin encabezados** (igual que en
interpolación). Ejemplo (los datos del parcial):

```
1.0 -0.148
1.2 -0.040
1.5 0.181
1.75 0.419
2.0 0.700
```

- Columna 1: `xi`, columna 2: `yi`.
- Máximo 50 puntos (`MAX_PUNTOS`).
- Usar **punto** decimal (`.`), no coma.

---

## 3) Qué ingresa y qué devuelve

| Entrada | Descripción |
|---|---|
| `grado` | Grado del polinomio de ajuste (1 = lineal, 2 = cuadrática, …) |
| `data.txt` | Los pares (xi, yi) |

Salida:

1. Los puntos leídos (para verificar).
2. La matriz y el vector del sistema (ya triangulados por Gauss).
3. El determinante de la matriz.
4. Los **coeficientes** `a[0], a[1], ..., a[k]`.
5. El polinomio de ajuste: `P(x) = a0 + a1·x + a2·x² + ...`.
6. El **coeficiente de correlación `r`** y una interpretación del ajuste.

Con los datos del parcial:

| Grado | Coeficientes | r | Interpretación |
|---|---|---|---|
| 1 (lineal) | a0 = −1.038775, a1 = 0.846426 | 0.992172 | excelente |
| 2 (cuadrática) | a0 = −0.242237, a1 = −0.284133, a2 = 0.377790 | 0.999996 | excelente |

> La cuadrática ajusta mejor: `r` más cerca de 1 y `Sr` muchísimo menor.

---

## 4) Teoría: qué minimizamos

Se busca minimizar la suma de los cuadrados de los **residuos**:

```
Sr = Σ (yi − P(xi))²
```

Derivando respecto de cada `aj` e igualando a cero se obtienen las
**ecuaciones normales** (sistema lineal con matriz simétrica):

```
para grado k:
  Σ (yi − a0 − a1·xi − ...) · xi^j = 0        para j = 0, 1, ..., k
```

que equivale a `A·a = b` con:

```
A[j][m] = Σ xi^(j+m)          (j, m = 0..k)
b[j]    = Σ yi·xi^j
```

Para grado 1 (lineal):

```
┌ n          Σx   ┐ ┌a0┐   ┌ Σy   ┐
│                 │ │  │ = │      │
└ Σx       Σx²   ┘ └a1┘   └ Σxy  ┘
```

Para grado 2 (cuadrática):

```
┌ n       Σx     Σx²  ┐ ┌a0┐   ┌ Σy   ┐
│ Σx     Σx²    Σx³   │ │a1│ = │ Σxy  │
└ Σx²    Σx³    Σx⁴   ┘ └a2┘   └ Σx²y ┘
```

Luego, con los coeficientes ya calculados:

```
St = Σ (yi − ȳ)²          con ȳ = promedio de las y
Sr = Σ (yi − P(xi))²
r  = √( (St − Sr) / St )
```

- `r` cercano a 1 → el modelo explica muy bien los datos.
- `Sr` chico → los residuos son pequeños.

---

## 5) Cálculo a mano — Regresión lineal (grado 1)

### Datos del parcial

| i | xi | yi |
|---|---|---|
| 1 | 1.00 | −0.148 |
| 2 | 1.20 | −0.040 |
| 3 | 1.50 | 0.181 |
| 4 | 1.75 | 0.419 |
| 5 | 2.00 | 0.700 |

### Paso 1: tabla de sumatorias

Cada sumatoria se arma **columna por columna**: un término por punto. Se
calcula una columna por cada potencia de `xi` que necesite el grado y una por
cada `yi·xi^k`:

| i | xi | yi | xi² | xi³ | xi⁴ | xi·yi | xi²·yi |
|---|---|---|---|---|---|---|---|
| 1 | 1.00 | −0.148 | 1.00 | 1.000 | 1.0000 | −0.14800 | −0.148000 |
| 2 | 1.20 | −0.040 | 1.44 | 1.728 | 2.0736 | −0.04800 | −0.057600 |
| 3 | 1.50 | 0.181 | 2.25 | 3.375 | 5.0625 | 0.27150 | 0.407250 |
| 4 | 1.75 | 0.419 | 3.0625 | 5.359375 | 9.37890625 | 0.73325 | 1.2831875 |
| 5 | 2.00 | 0.700 | 4.00 | 8.000 | 16.0000 | 1.40000 | 2.800000 |
| **Σx = 7.45** | | **Σy = 1.112** | **Σx² = 11.7525** | **Σx³ = 19.462375** | **Σx⁴ = 33.51500625** | **Σxy = 2.20875** | **Σx²y = 4.2848375** |

Cada total es la suma de los términos de su columna:

```
Σx   = 1.00 + 1.20 + 1.50 + 1.75 + 2.00                    = 7.45
Σy   = -0.148 - 0.040 + 0.181 + 0.419 + 0.700              = 1.112
Σx²  = 1.00 + 1.44 + 2.25 + 3.0625 + 4.00                  = 11.7525
Σx³  = 1.000 + 1.728 + 3.375 + 5.359375 + 8.000            = 19.462375
Σx⁴  = 1.0000 + 2.0736 + 5.0625 + 9.37890625 + 16.0000     = 33.51500625
Σxy  = -0.14800 - 0.04800 + 0.27150 + 0.73325 + 1.40000    = 2.20875
Σx²y = -0.148000 - 0.057600 + 0.407250 + 1.2831875 + 2.8   = 4.2848375
```

Para la regresión **lineal** solo se usan: `n`, `Σx`, `Σx²`, `Σy`, `Σxy`.
Para la **cuadrática** se agregan: `Σx³`, `Σx⁴`, `Σx²y`.

n = 5, x̄ = 7.45 / 5 = 1.49, ȳ = 1.112 / 5 = 0.2224

### Paso 2: armar el sistema de ecuaciones normales

```
5·a0  + 7.45·a1  = 1.112
7.45·a0 + 11.7525·a1 = 2.20875
```

### Paso 3: resolver el sistema

Con eliminación: fila 2 − 1.49·fila 1 (porque 7.45/5 = 1.49):

```
(7.45 − 1.49·5)·a0 + (11.7525 − 1.49·7.45)·a1 = 2.20875 − 1.49·1.112
       0·a0        +        0.652·a1          = 0.55187
```

```
a1 = 0.55187 / 0.652 = 0.846426
```

Despejando de la primera ecuación:

```
a0 = (1.112 − 7.45·0.846426) / 5
   = (1.112 − 6.305874) / 5
   = −1.038775
```

> Fórmula rápida equivalente:
> a1 = (n·Σxy − Σx·Σy) / (n·Σx² − (Σx)²)
>    = (5·2.20875 − 7.45·1.112) / (5·11.7525 − 7.45²)
>    = 2.75935 / 3.26 = 0.846427
> a0 = ȳ − a1·x̄ = 0.2224 − 0.846427·1.49 = −1.038776

### Paso 4: el polinomio de ajuste

```
y = −1.038775 + 0.846426·x
```

### Paso 5: predecir, calcular residuos, Sr, St y r

ŷi = −1.038775 + 0.846426·xi

| xi | yi | ŷi | ei = yi − ŷi | ei² |
|---|---|---|---|---|
| 1.00 | −0.148 | −0.192349 | +0.044349 | 0.0019668 |
| 1.20 | −0.040 | −0.023064 | −0.016936 | 0.0002868 |
| 1.50 | 0.181 | 0.230864 | −0.049864 | 0.0024864 |
| 1.75 | 0.419 | 0.442471 | −0.023471 | 0.0005509 |
| 2.00 | 0.700 | 0.654077 | +0.045923 | 0.0021089 |

```
Sr = 0.0019668 + 0.0002868 + 0.0024864 + 0.0005509 + 0.0021089 = 0.0073999
```

St con ȳ = 0.2224:

```
St = (−0.3704)² + (−0.2624)² + (−0.0414)² + (0.1966)² + (0.4776)²
   = 0.1371962 + 0.0688538 + 0.0017140 + 0.0386516 + 0.2281018
   = 0.4745172
```

```
r = √( (0.4745172 − 0.0073999) / 0.4745172 )
  = √(0.98441)
  = 0.992172
```

---

## 6) Cálculo a mano — Regresión cuadrática (grado 2)

### Paso 1: sumatorias (agregando potencias mayores)

De la tabla anterior, además:

```
Σx³ = 1³ + 1.2³ + 1.5³ + 1.75³ + 2³
    = 1 + 1.728 + 3.375 + 5.359375 + 8 = 19.462375

Σx⁴ = 1 + 2.0736 + 5.0625 + 9.37890625 + 16 = 33.51500625

Σx²y = 1²·(−0.148) + 1.2²·(−0.040) + 1.5²·0.181 + 1.75²·0.419 + 2²·0.700
     = −0.148 − 0.0576 + 0.40725 + 1.2831875 + 2.8 = 4.2848375
```

### Paso 2: sistema de ecuaciones normales (3×3)

```
┌  5        7.45       11.7525    ┐ ┌a0┐   ┌  1.112    ┐
│  7.45     11.7525    19.462375  │ │a1│ = │  2.20875  │
└  11.7525  19.462375  33.51500625┘ └a2┘   └  4.2848375┘
```

### Paso 3: resolver con Gauss (a mano)

**a) Eliminar debajo del pivote a00 = 5 (multiplicadores 7.45/5 = 1.49 y 11.7525/5 = 2.3505):**

Fila 2 ← Fila 2 − 1.49·Fila 1:

```
0·a0 + (11.7525 − 1.49·7.45)·a1 + (19.462375 − 1.49·11.7525)·a2 = 2.20875 − 1.49·1.112
  0·a0 +          0.652·a1       +         1.95115·a2           = 0.55187
```

Fila 3 ← Fila 3 − 2.3505·Fila 1:

```
0·a0 + (19.462375 − 2.3505·7.45)·a1 + (33.51500625 − 2.3505·11.7525)·a2
     = 4.2848375 − 2.3505·1.112

0·a0 +        1.95115·a1            +        5.890755·a2        = 1.6710815
```

**b) Eliminar debajo del pivote 0.652 (multiplicador 1.95115/0.652 ≈ 2.9926):**

Fila 3 ← Fila 3 − 2.9926·Fila 2:

```
0·a0 + 0·a1 + (5.890755 − 2.9926·1.95115)·a2 = 1.6710815 − 2.9926·0.55187
       0·a0 + 0·a1 +        0.051819·a2       = 0.019577
```

**c) Sustitución hacia atrás:**

```
a2 = 0.019577 / 0.051819 = 0.377790

a1 = (0.55187 − 1.95115·0.377790) / 0.652
   = (0.55187 − 0.737125) / 0.652
   = −0.284133

a0 = (1.112 − 7.45·(−0.284133) − 11.7525·0.377790) / 5
   = (1.112 + 2.116791 − 4.439997) / 5
   = −0.242237
```

### Paso 4: el polinomio de ajuste

```
y = −0.242237 − 0.284133·x + 0.377790·x²
```

### Paso 5: residuos, Sr y r

ŷi = −0.242237 − 0.284133·xi + 0.377790·xi²

| xi | yi | ŷi | ei = yi − ŷi |
|---|---|---|---|
| 1.00 | −0.148 | −0.148580 | +0.000580 |
| 1.20 | −0.040 | −0.039179 | −0.000821 |
| 1.50 | 0.181 | 0.181591 | −0.000591 |
| 1.75 | 0.419 | 0.417512 | +0.001488 |
| 2.00 | 0.700 | 0.700657 | −0.000657 |

```
Sr = (0.000580)² + (0.000821)² + (0.000591)² + (0.001488)² + (0.000657)²
   ≈ 0.00000401

r = √( (0.4745172 − 0.00000401) / 0.4745172 )
  = √(0.99999156)
  = 0.999996
```

**Conclusión:** r cuadrático (0.999996) > r lineal (0.992172) y Sr cuadrático
(~4·10⁻⁶) mucho menor que el lineal (~7.4·10⁻³) → **la parábola ajusta mejor**.

---

## 7) Si te dan el polinomio de ajuste directamente

Cuando el enunciado ya te da la curva (`y = a0 + a1·x` o `y = a0 + a1·x + a2·x²`)
y los puntos, lo que suele pedirse es:

1. **Evaluar el polinomio** en cada xi: `ŷi = P(xi)`.
2. **Calcular los residuos**: `ei = yi − ŷi`.
3. **Calcular Sr = Σ ei²** y **St = Σ (yi − ȳ)²**.
4. **Calcular r = √( (St − Sr) / St )** y comparar modelos: el de mayor **r**
   (más cercano a 1) es el mejor ajuste.

Si además te dan el polinomio **para comparar con el que vos obtuviste**, los
coeficientes deben coincidir; pequeñas diferencias en las últimas cifras son
normales por redondeo.

---

## 8) Limitaciones del código

- **Solo polinomios**: ajusta `a0 + a1·x + ... + ak·x^k`. No sirve para modelos
  exponenciales, potenciales o logarítmicos (ahí habría que linealizar los datos
  o cambiar el modelo).
- **Se necesitan al menos `grado + 1` puntos**, si no, el programa avisa y corta.
- Si `n = grado + 1` el ajuste **interpola exactamente** (pasa por todos los
  nodos, Sr = 0, r = 1).
- **Grados altos (≥ 6–7)** pueden dar coeficientes con error por mal
  condicionamiento de las ecuaciones normales. Para parciales (grados 1–2)
  trabaja exacto.