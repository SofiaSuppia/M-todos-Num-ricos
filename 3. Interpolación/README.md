# Interpolación de Lagrange

## Archivos de esta carpeta

| Archivo | Qué es |
|---|---|
| `LagrangePolinomio.c++` | Código fuente: interpolación de Lagrange + polinomio interpolante por Gauss |
| `LagrangePolinomio` | Ejecutable ya compilado |
| `data.txt` | Archivo de datos de entrada: una pareja `xi yi` por línea |
| `GraficarLagrange.py` | Script Python que grafica los coeficientes C~N,k~(x), Pn(x) contra f(x) y el error |
| `README.md` | Este documento |

---

## 1) Cómo se ejecuta

### Compilar (una sola vez, si modificaste el .c++)

```bash
g++ -O2 -o LagrangePolinomio LagrangePolinomio.c++ -lm
```

### Ejecutar

```bash
./LagrangePolinomio
```

> El programa **lee `data.txt` de la carpeta donde estás parado** (ruta relativa).
> Si se ejecuta desde otra carpeta, hay que copiar `data.txt` allí o ejecutarlo
> desde la carpeta donde está el archivo.

### Formato de `data.txt`

Una línea por punto, separados por espacio (o tab), **sin encabezados**:

```
1.0 -0.148
1.2 -0.040
1.5 0.181
1.75 0.419
2.0 0.700
```

- Primera columna: `xi`, segunda columna: `yi`.
- Máximo de puntos: `MAX_PUNTOS = 50`.
- Usar **punto** (`.`) como separador decimal, no coma.

---

## 2) Valores que se ingresan

Después de leer los datos, aparece un menú:

```
--- MENU DE INTERPOLACION ---
1. Lagrange (valor interpolado)
2. Polinomio interpolante (con Gauss y pivoteo)
3. Lagrange expandido (coeficientes del polinomio)
0. Salir
```

| Opción | Qué hace | Qué ingresás |
|---|---|---|
| **1** | Calcula Pn(x̂) con la fórmula de Lagrange (numérico) | el valor `x` (xsombrero) |
| **2** | Arma la matriz de Vandermonde `A[i][j]=xi^j`, la resuelve con Gauss + pivoteo y da los coeficientes a0..a4 | primero el `x` a evaluar |
| **3** | Muestra el polinomio expandido (forma clásica a0 + a1·x + ... ) | el valor `x` a interpolar |

- Si `x̂` está **dentro del rango** de los datos → interpolación.
- Si `x̂` está **fuera del rango** → extrapolación (menos confiable).
- En la opción 1 y 3 el programa imprime cada coeficiente de Lagrange
  `L(n-1)k(x̂)`: la suma de todos siempre da **1** (sirve para verificar).

Ejemplo con el parcial (5 puntos, x̂ = 1.6):

```
L40(1.600) =  0.032000
L41(1.600) = -0.136364
L42(1.600) =  0.768000
L43(1.600) =  0.372364
L44(1.600) = -0.036000
suma = 1

P(1.6) = -0.148·0.032 - 0.040·0.136364 + 0.181·0.768 + 0.419·0.372364 + 0.700·(-0.036)
P(1.6) = 0.270547        (el parcial da 0.2705469091 → coincide)
```

---

## 3) Cálculo a mano (fórmula de Lagrange)

Con n+1 puntos `(x0,y0), ..., (xn,yn)`:

```
Pn(x) = Σ(k=0..n)  yk · Lk(x)

             (x - x0)(x - x1)...(x - x[k-1])(x - x[k+1])...(x - xn)
Lk(x) =  ------------------------------------------------------------------
             (xk - x0)(xk - x1)...(xk - x[k-1])(xk - x[k+1])...(xk - xn)
```

Es decir, en cada `Lk` se **omite el punto k** (el numerador y el denominador
tienen todos los nodos salvo el `k`).

**Propiedades útiles:**
- `Lk(xk) = 1` y `Lk(xj) = 0` para `j ≠ k`.
- `Σ Lk(x) = 1` para cualquier x.

### Ejemplo a mano (parcial, calcular P4(1.6))

Nodos: x0=1.0, x1=1.2, x2=1.5, x3=1.75, x4=2.0
y: −0.148, −0.040, 0.181, 0.419, 0.700

Cada Lk se arma **omitiendo el punto k**. Ejemplo completo para k=3 (se omite x3=1.75):

```
L3(1.6) = (1.6-1)(1.6-1.2)(1.6-1.5)(1.6-2)
          -----------------------------------
          (1.75-1)(1.75-1.2)(1.75-1.5)(1.75-2)

Numerador   = 0.6 · 0.4 · 0.1 · (-0.4) = -0.0096
Denominador = 0.75 · 0.55 · 0.25 · (-0.25) = -0.02578125

L3(1.6) = -0.0096 / -0.02578125 = 0.372364
```

Se repite para k = 0, 1, 2, 4 (los cinco términos de la **sumatoria**):

| k | Lk(1.6) = numerador / denominador | resultado |
|---|---|---|
| 0 | 0.4·0.1·(−0.15)·(−0.4) / (−0.2)·(−0.5)·(−0.75)·(−1.0) = 0.0024 / 0.075 | **0.032000** |
| 1 | 0.6·0.1·(−0.15)·(−0.4) / 0.2·(−0.3)·(−0.55)·(−0.8) = 0.0036 / (−0.0264) | **−0.136364** |
| 2 | 0.6·0.4·(−0.15)·(−0.4) / 0.5·0.3·(−0.25)·(−0.5) = 0.0144 / 0.01875 | **0.768000** |
| 3 | 0.6·0.4·0.1·(−0.4) / 0.75·0.55·0.25·(−0.25) = −0.0096 / (−0.02578125) | **0.372364** |
| 4 | 0.6·0.4·0.1·(−0.15) / 1.0·0.8·0.5·0.25 = −0.0036 / 0.1 | **−0.036000** |

Verificación: **Σ Lk = 0.032 − 0.136364 + 0.768 + 0.372364 − 0.036 = 1.0** ✓

Ahora se resuelve la sumatoria `P4(1.6) = Σ yk · Lk(1.6)` término a término:

```
P4(1.6) = y0·L0 + y1·L1 + y2·L2 + y3·L3 + y4·L4
        = (-0.148)·(0.032) + (-0.040)·(-0.136364) + (0.181)·(0.768)
          + (0.419)·(0.372364) + (0.700)·(-0.036)
        = -0.004736 + 0.005455 + 0.139008 + 0.156021 - 0.025200
        = 0.270547
```

### Expandir el polinomio

Haciendo las cuentas (o con la opción 3 del programa) el polinomio queda en
forma clásica (parcial, x0=1.0 ... x4=2.0):

```
P4(x) = -0.39109091 + 0.21493939·x - 0.22293939·x² + 0.30806061·x³ - 0.05696970·x⁴
```

Con ella también se puede evaluar: P4(1.6) = 0.2705469091.

---

## 4) Si te dan el polinomio directamente

Cuando el enunciado **te da Pn(x) ya armado** (por ejemplo la forma expandida),
lo que tenés que hacer normalmente es:

1. **Verificar que interpola**: comprobar que `P(xk) = yk` en todos los nodos
   (si falla en algún nodo, el polinomio no corresponde a esos datos).
2. **Evaluar** en el xsombrero pedido: `P(x̂) = Σ ak·x̂^k` (o con el esquema de Horner).
3. **Comparar contra f(x̂)** si conocés la f real para estimar el error
   `e(x̂) = |f(x̂) − P(x̂)|`.

Si el ejercicio pide "hallar el polinomio interpolante" y además te dan Pn(x),
lo que se compara es el resultado del código (opción 2 o 3) contra el del
enunciado: deben coincidir coeficiente a coeficiente.

> Nota de precisión: la opción 2 (Vandermonde + Gauss) puede dar una diferencia
> mínima en las últimas cifras por el condicionamiento de la matriz; el valor
> interpolado P(x̂) sigue coincidiendo.

---

## 5) Script de gráficos

`GraficarLagrange.py` (sin argumentos) grafica, para los nodos (1, 2, 2.5):

1. Los coeficientes C~2,k~(x) y verifica que valen 1 en su nodo y 0 en los otros.
2. «f(x) vs P₂(x)»: el polinomio pasa exactamente por los nodos.
3. El error |f(x) − P₂(x)|: vale 0 en los nodos.

Para graficar **otra tabla**, usar `GrafDatosTabla.py` (en la raíz del
proyecto), que lee `data.txt` y luego, en la consola del script:

```python
graficar_tabla(x, y, xmin=0.5, xmax=2.5, ymin=-0.3, ymax=0.9)
```

Ajustá los límites según tus datos para que queden visibles.