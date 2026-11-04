# Problema 4 — Examen 2025 (Interpolación de Lagrange — Radar Doppler)

Resolución **completa a mano**, paso a paso, con todos los datos y cálculos
que hay que escribir en cada parte.

## Archivos de esta carpeta

| Archivo | Qué es |
|---|---|
| `problema4.c++` | Código: Lagrange (valor numérico + polinomio expandido) + interpolación cúbica. Enunciado en comentarios |
| `problema4` | Ejecutable ya compilado |
| `data.txt` | Datos: 5 pares `ti vi` (sin encabezado, una línea por punto) |
| `README - Problema 4.md` | Este documento |

Compilar y ejecutar:

```bash
g++ -Wall -O2 -o problema4 problema4.c++ -lm
./problema4
```

---

## Enunciado

Un radar Doppler mide la velocidad de caída de un paracaidista. Los datos:

| Tiempo (s) | Velocidad (cm/s) |
|---:|---:|
| 1 | 800 |
| 3 | 2310 |
| 5 | 3090 |
| 7 | 3940 |
| 13 | 4755 |

**a)** Empleando interpolación de Lagrange, estime las velocidades para
`t = 4 s` y `t = 6 s`. (Adjunte captura de pantalla de la consola).

**b)** Obtenga el polinomio de Lagrange que interpola los datos de la tabla.

**c)** Use interpolación cúbica (normal) para estimar la velocidad en
`t = 10 s`, y escriba el polinomio usado. (Escriba ambas cosas en papel
y adjunte captura de pantalla de la consola).

---

## Parte a) Interpolación de Lagrange: v(4) y v(6)

### Fórmula de Lagrange (repaso)

Con `n+1` nodos `(x₀,y₀), ..., (xₙ,yₙ)`:

```
Pₙ(x) = Σₖ₌₀ⁿ  yₖ · Lₖ(x)

              Πᵢ₌₀ (x - xᵢ)       (omitir i = k)
Lₖ(x) =  ──────────────────
              Πᵢ₌₀ (xₖ - xᵢ)      (omitir i = k)
```

Propiedades: `Lₖ(xₖ) = 1`, `Lₖ(xⱼ) = 0` si `j ≠ k`, y `Σ Lₖ(x) = 1`.

Nodos: `t₀=1`, `t₁=3`, `t₂=5`, `t₃=7`, `t₄=13` → `n = 4` (grado 4).

---

### Estimación de v(4) (grado 4, todos los nodos)

#### Cálculo de cada Lₖ(4)

**L₀(4)** — omite `t₀ = 1`:

```
Numerador   = (4-3)(4-5)(4-7)(4-13) = (1)(-1)(-3)(-9) = -27
Denominador = (1-3)(1-5)(1-7)(1-13) = (-2)(-4)(-6)(-12) = 576

L₀(4) = -27 / 576 = -0.046875
```

**L₁(4)** — omite `t₁ = 3`:

```
Numerador   = (4-1)(4-5)(4-7)(4-13) = (3)(-1)(-3)(-9) = -81
Denominador = (3-1)(3-5)(3-7)(3-13) = (2)(-2)(-4)(-10) = -160

L₁(4) = -81 / (-160) = 0.506250
```

**L₂(4)** — omite `t₂ = 5`:

```
Numerador   = (4-1)(4-3)(4-7)(4-13) = (3)(1)(-3)(-9) = 81
Denominador = (5-1)(5-3)(5-7)(5-13) = (4)(2)(-2)(-8) = 128

L₂(4) = 81 / 128 = 0.632812
```

**L₃(4)** — omite `t₃ = 7`:

```
Numerador   = (4-1)(4-3)(4-5)(4-13) = (3)(1)(-1)(-9) = 27
Denominador = (7-1)(7-3)(7-5)(7-13) = (6)(4)(2)(-6) = -288

L₃(4) = 27 / (-288) = -0.093750
```

**L₄(4)** — omite `t₄ = 13`:

```
Numerador   = (4-1)(4-3)(4-5)(4-7) = (3)(1)(-1)(-3) = 9
Denominador = (13-1)(13-3)(13-5)(13-7) = (12)(10)(8)(6) = 5760

L₄(4) = 9 / 5760 = 0.001563
```

**Verificación:** `Σ Lₖ(4) = -0.046875 + 0.50625 + 0.632812 - 0.09375 + 0.001563 = 1.000000` ✓

#### Valor interpolado

```
P₄(4) = 800·(-0.046875) + 2310·(0.50625) + 3090·(0.632812)
       + 3940·(-0.09375) + 4755·(0.001563)

     = (-37.5) + 1169.4375 + 1955.390625 + (-369.375) + 7.4296875

     = 2725.382812 cm/s
```

> **v(4) ≈ 2725.4 cm/s** (entre v(3)=2310 y v(5)=3090, coherente).

---

### Estimación de v(6) (grado 4, todos los nodos)

#### Cálculo de cada Lₖ(6)

**L₀(6)** — omite `t₀ = 1`:

```
Numerador   = (6-3)(6-5)(6-7)(6-13) = (3)(1)(-1)(-7) = 21
Denominador = 576   (mismo que antes)

L₀(6) = 21 / 576 = 0.036458
```

**L₁(6)** — omite `t₁ = 3`:

```
Numerador   = (6-1)(6-5)(6-7)(6-13) = (5)(1)(-1)(-7) = 35
Denominador = -160   (mismo)

L₁(6) = 35 / (-160) = -0.218750
```

**L₂(6)** — omite `t₂ = 5`:

```
Numerador   = (6-1)(6-3)(6-7)(6-13) = (5)(3)(-1)(-7) = 105
Denominador = 128   (mismo)

L₂(6) = 105 / 128 = 0.820312
```

**L₃(6)** — omite `t₃ = 7`:

```
Numerador   = (6-1)(6-3)(6-5)(6-13) = (5)(3)(1)(-7) = -105
Denominador = -288   (mismo)

L₃(6) = -105 / (-288) = 0.364583
```

**L₄(6)** — omite `t₄ = 13`:

```
Numerador   = (6-1)(6-3)(6-5)(6-7) = (5)(3)(1)(-1) = -15
Denominador = 5760   (mismo)

L₄(6) = -15 / 5760 = -0.002604
```

**Verificación:** `Σ Lₖ(6) = 0.036458 - 0.21875 + 0.820312 + 0.364583 - 0.002604 = 1.000000` ✓

#### Valor interpolado

```
P₄(6) = 800·(0.036458) + 2310·(-0.21875) + 3090·(0.820312)
       + 3940·(0.364583) + 4755·(-0.002604)

     = 29.1667 - 505.3125 + 2534.7656 + 1436.4583 - 12.3828

     = 3482.695312 cm/s
```

> **v(6) ≈ 3482.7 cm/s** (entre v(5)=3090 y v(7)=3940, coherente).

---

## Parte b) Polinomio de Lagrange que interpola la tabla (grado 4)

### Polinomio en forma de Lagrange (suma de productos)

```
P₄(t) = 800·L₀(t) + 2310·L₁(t) + 3090·L₂(t) + 3940·L₃(t) + 4755·L₄(t)
```

con cada `Lₖ(t)` como se calculó arriba (Producto de factores, omitiendo el
término k).

### Polinomio expandido (forma clásica)

Al expandir todos los productos y simplificar (lo hace el código con
`lagrangeCoeficientes`, o a mano con el esquema de Horner):

```
P₄(t) = -663.867188 + 1813.625000·t - 392.869792·t² + 44.875000·t³ - 1.763021·t⁴
```

### Verificación: el polinomio reproduce los datos

| t | P₄(t) | v dato | ¿Correcto? |
|---:|---:|---:|:---:|
| 1 | 800.000 | 800 | ✓ |
| 3 | 2310.000 | 2310 | ✓ |
| 5 | 3090.000 | 3090 | ✓ |
| 7 | 3940.000 | 3940 | ✓ |
| 13 | 4755.000 | 4755 | ✓ |

> El polinomio pasa **exactamente** por los 5 puntos.

---

## Parte c) Interpolación cúbica en t = 10 s

### ¿Por qué estos nodos?

`t = 10` s está entre `7` y `13` s (los datos no son equiespaciados: hay un
salto de 6 s entre el penúltimo y el último punto). Para una interpolación
**cúbica** (grado 3) necesitamos exactamente **4 nodos**.

El criterio es elegir los **4 nodos más cercanos** a `t = 10`:

```
|1 - 10| = 9    |3 - 10| = 7    |5 - 10| = 5
|7 - 10| = 3    |13 - 10| = 3
```

Los 4 más cercanos: **3, 5, 7 y 13** (distancias 7, 5, 3, 3). `t = 10`
queda interpolado entre `7` y `13`, y los otros dos nodos aportan la
curvatura necesaria para el polinomio cúbico.

Nodos seleccionados:

| Índice j | tⱼ (s) | vⱼ (cm/s) |
|---:|---:|---:|
| 0 | 3 | 2310 |
| 1 | 5 | 3090 |
| 2 | 7 | 3940 |
| 3 | 13 | 4755 |

### Cálculo de cada Lₖ(10) — polinomio de Lagrange de grado 3

Denominadores (cada uno se calcula una sola vez con los nodos 3, 5, 7, 13):

```
Den₀ = (3-5)(3-7)(3-13)  = (-2)(-4)(-10) = -80
Den₁ = (5-3)(5-7)(5-13)  = (2)(-2)(-8)   =  32
Den₂ = (7-3)(7-5)(7-13)  = (4)(2)(-6)    = -48
Den₃ = (13-3)(13-5)(13-7) = (10)(8)(6)   = 480
```

**L₀(10)** — omite nodo `t₀ = 3`:

```
Num = (10-5)(10-7)(10-13) = (5)(3)(-3) = -45

L₀(10) = -45 / (-80) = 0.562500
```

**L₁(10)** — omite nodo `t₁ = 5`:

```
Num = (10-3)(10-7)(10-13) = (7)(3)(-3) = -63

L₁(10) = -63 / 32 = -1.968750
```

**L₂(10)** — omite nodo `t₂ = 7`:

```
Num = (10-3)(10-5)(10-13) = (7)(5)(-3) = -105

L₂(10) = -105 / (-48) = 2.187500
```

**L₃(10)** — omite nodo `t₃ = 13`:

```
Num = (10-3)(10-5)(10-7) = (7)(5)(3) = 105

L₃(10) = 105 / 480 = 0.218750
```

**Verificación:** `Σ Lₖ(10) = 0.5625 - 1.96875 + 2.1875 + 0.21875 = 1.000000` ✓

### Valor interpolado

```
P₃(10) = 2310·(0.5625) + 3090·(-1.96875) + 3940·(2.1875) + 4755·(0.21875)

      = 1299.375 - 6083.4375 + 8618.75 + 1040.15625

      = 4874.84375 cm/s
```

> **v(10) ≈ 4874.84 cm/s**

### Polinomio cúbico usado (forma clásica expandida)

Al expandir el polinomio de Lagrange de grado 3 con esos 4 nodos:

```
P₃(t) = 1742.656250 + 1.239583·t + 76.093750·t² - 4.489583·t³
```

Se puede evaluar directamente `P₃(10)` con esta forma para verificar:

```
P₃(10) = 1742.656250 + 1.239583·10 + 76.093750·100 - 4.489583·1000
       = 1742.656250 + 12.395830 + 7609.375000 - 4489.583000
       = 4874.844080  ≈  4874.84375  (diferencia por redondeo) ✓
```

---

## Resumen

| Inciso | Qué se hace | Resultado |
|---|---|---|
| **a)** | Lagrange grado 4, v(4) y v(6) con todos los nodos | v(4) = **2725.38 cm/s**, v(6) = **3482.70 cm/s** |
| **b)** | Polinomio de Lagrange interpolante (grado 4) | P₄(t) = -663.867 + 1813.625·t - 392.870·t² + 44.875·t³ - 1.763·t⁴ |
| **c)** | Interpolación cúbica en t = 10 (nodos 3, 5, 7, 13) | v(10) = **4874.84 cm/s**, P₃(t) = 1742.656 + 1.240·t + 76.094·t² - 4.490·t³ |