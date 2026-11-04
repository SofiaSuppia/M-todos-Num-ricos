# Problema 2 — Segundo Parcial 2025 (Diferenciación numérica e integración)

Resolución completa del problema, con el desarrollo manual paso a paso y la
verificación con el programa `problema2.c++`.

## Archivos de esta carpeta

| Archivo | Qué es |
|---|---|
| `problema2.c++` | Programa: derivada exacta, diferencias centradas para \(h=0.5\) y \(h=0.01\), tabla de derivadas con \(h=0.1\) y trapecio compuesto. |
| `problema2` | Ejecutable ya compilado. |
| `salida_problema2.txt` | Salida completa del programa, útil para la captura de consola. |
| `captura_problema2.png` | Captura tipo consola de la salida del programa. |
| `README - Problema 2.md` | Este documento. |

Compilar y ejecutar:

```bash
g++ -Wall -Wextra -O2 -o problema2 problema2.c++ -lm
./problema2
```

---

## Enunciado

Sea

\[
f(x)=\sqrt{1+x^2}.
\]

Se pide calcular \(f'(1)\).

**a)** Calcular la derivada exacta y evaluarla en \(x=1\).

**b)** Aproximar \(f'(1)\) mediante una **diferencia finita centrada de orden
\(h^2\)** con \(h=0.5\). Escribir el resultado con 6 cifras decimales, calcular
el error absoluto exacto y analizar la consistencia del resultado.

**c)** Repetir el inciso b) con \(h=0.01\).

**d)** Construir una tabla de \(f'(x)\) para

\[
x=0.0,0.1,0.2,\ldots,1.0
\]

con \(h=0.1\). Usar:

- diferencia hacia adelante en el primer punto,
- diferencia hacia atrás en el último punto,
- diferencia centrada en los puntos interiores.

Escribir los valores con 4 cifras decimales.

**e)** Usando la tabla del inciso d), calcular

\[
\int_0^1 f'(x)\,dx
\]

con **trapecio compuesto**. Escribir el resultado con 10 cifras decimales y
calcular el error absoluto exacto.

---

## Inciso a) Derivada exacta

\[
f(x)=(1+x^2)^{1/2}.
\]

Derivando con la regla de la cadena:

\[
f'(x)=\frac{1}{2}(1+x^2)^{-1/2}\cdot 2x
\]

\[
\boxed{f'(x)=\frac{x}{\sqrt{1+x^2}}}
\]

En \(x=1\):

\[
f'(1)=\frac{1}{\sqrt{1+1^2}}=\frac{1}{\sqrt2}.
\]

Racionalizando:

\[
f'(1)=\frac{\sqrt2}{2}.
\]

Numéricamente:

\[
\boxed{f'(1)=0.7071067812}
\]

Con 6 cifras decimales:

\[
\boxed{f'(1)\approx 0.707107}
\]

---

## Inciso b) Diferencia finita centrada con \(h=0.5\)

### Fórmula

La diferencia centrada de orden \(h^2\) es

\[
f'(x_i)\approx \frac{f(x_i+h)-f(x_i-h)}{2h}.
\]

Su error de truncamiento es

\[
E_T=-\frac{f'''(\xi)}{6}h^2,
\]

por lo que el método es de **orden \(h^2\)**.

Para \(x_i=1\) y \(h=0.5\):

\[
f'(1)\approx \frac{f(1.5)-f(0.5)}{2(0.5)}
=\frac{f(1.5)-f(0.5)}{1}.
\]

### Evaluaciones

\[
f(1.5)=\sqrt{1+1.5^2}=\sqrt{3.25}=1.8027756377,
\]

\[
f(0.5)=\sqrt{1+0.5^2}=\sqrt{1.25}=1.1180339887.
\]

### Resultado

\[
f'(1)\approx 1.8027756377-1.1180339887
\]

\[
\boxed{f'(1)\approx 0.6847416490}
\]

Con 6 cifras decimales:

\[
\boxed{f'(1)\approx 0.684742}
\]

### Error absoluto exacto

\[
E=\left|0.6847416490-0.7071067812\right|
\]

\[
\boxed{E=0.0223651322}
\]

Con 6 cifras decimales:

\[
\boxed{E\approx 0.022365}
\]

### Análisis de consistencia

La tercera derivada es

\[
f'''(x)=-\frac{3x}{(1+x^2)^{5/2}}.
\]

Cerca de \(x=1\):

\[
|f'''(1)|=\frac{3}{2^{5/2}}=0.5303300859.
\]

El error teórico aproximado para \(h=0.5\) es

\[
|E_T|\approx \frac{0.5303300859}{6}(0.5)^2
=0.0220970869.
\]

El error obtenido, \(0.0223651322\), es del mismo orden que el error teórico.
Esto confirma que el resultado es consistente con un método de orden \(h^2\).

---

## Inciso c) Diferencia finita centrada con \(h=0.01\)

Para \(x_i=1\) y \(h=0.01\):

\[
f'(1)\approx \frac{f(1.01)-f(0.99)}{2(0.01)}
=\frac{f(1.01)-f(0.99)}{0.02}.
\]

### Evaluaciones

\[
f(1.01)=\sqrt{1+1.01^2}=1.4213022198,
\]

\[
f(0.99)=\sqrt{1+0.99^2}=1.4071602610.
\]

### Resultado

\[
f'(1)\approx \frac{1.4213022198-1.4071602610}{0.02}
\]

\[
f'(1)\approx \frac{0.0141419588}{0.02}
\]

\[
\boxed{f'(1)\approx 0.7070979423}
\]

Con 6 cifras decimales:

\[
\boxed{f'(1)\approx 0.707098}
\]

### Error absoluto exacto

\[
E=\left|0.7070979423-0.7071067812\right|
\]

\[
\boxed{E=0.0000088389}
\]

Con 6 cifras decimales:

\[
\boxed{E\approx 0.000009}
\]

### Comparación con \(h=0.5\)

\[
\frac{E(0.5)}{E(0.01)}
=
\frac{0.0223651322}{0.0000088389}
=2530.31.
\]

Si el método es de orden \(h^2\), al dividir \(h\) por 50 el error debería
dividirse aproximadamente por

\[
\left(\frac{0.5}{0.01}\right)^2=50^2=2500.
\]

El cociente obtenido, \(2530.31\), es muy cercano a \(2500\). La diferencia se
debe a que la tercera derivada se evalúa en puntos distintos y a los errores de
redondeo. Por lo tanto, el comportamiento observado confirma la consistencia
del método de orden \(h^2\).

---

## Inciso d) Tabla de derivadas con \(h=0.1\)

### Fórmulas usadas

Primer punto, diferencia hacia adelante:

\[
f'(x_0)\approx \frac{f(x_1)-f(x_0)}{h}.
\]

Puntos interiores, diferencia centrada:

\[
f'(x_i)\approx \frac{f(x_{i+1})-f(x_{i-1})}{2h}.
\]

Último punto, diferencia hacia atrás:

\[
f'(x_n)\approx \frac{f(x_n)-f(x_{n-1})}{h}.
\]

### Ejemplos de cálculo

**Primer punto, \(x=0.0\):**

\[
f'(0.0)\approx \frac{f(0.1)-f(0.0)}{0.1}
=\frac{1.0049875621-1}{0.1}
=0.0498756211\approx 0.0499.
\]

**Punto interior, \(x=0.5\):**

\[
f'(0.5)\approx \frac{f(0.6)-f(0.4)}{2(0.1)}
=\frac{1.1661903790-1.0770338119}{0.2}
=0.4457828355\approx 0.4458.
\]

**Último punto, \(x=1.0\):**

\[
f'(1.0)\approx \frac{f(1.0)-f(0.9)}{0.1}
=\frac{1.4142135624-1.3453624047}{0.1}
=0.6885115766\approx 0.6885.
\]

### Tabla pedida (4 cifras decimales)

| \(x\) | \(f'(x)\) aproximada | \(f'(x)\) exacta |
|---:|---:|---:|
| 0.0 | 0.0499 | 0.0000000000 |
| 0.1 | 0.0990 | 0.0995037190 |
| 0.2 | 0.1952 | 0.1961161351 |
| 0.3 | 0.2861 | 0.2873478856 |
| 0.4 | 0.3700 | 0.3713906764 |
| 0.5 | 0.4458 | 0.4472135955 |
| 0.6 | 0.5131 | 0.5144957554 |
| 0.7 | 0.5722 | 0.5734623444 |
| 0.8 | 0.6235 | 0.6246950476 |
| 0.9 | 0.6679 | 0.6689647316 |
| 1.0 | 0.6885 | 0.7071067812 |

\[
\boxed{
\begin{array}{c|ccccccccccc}
x & 0.0 & 0.1 & 0.2 & 0.3 & 0.4 & 0.5 & 0.6 & 0.7 & 0.8 & 0.9 & 1.0\\
\hline
f'(x) & 0.0499 & 0.0990 & 0.1952 & 0.2861 & 0.3700 & 0.4458 & 0.5131 & 0.5722 & 0.6235 & 0.6679 & 0.6885
\end{array}}
\]

En \(x=0\) la derivada exacta es \(0\), pero la aproximación hacia adelante da
\(0.0499\). Es el error propio de una fórmula de primer orden en un extremo.

En \(x=1\) la derivada exacta es \(0.707107\), mientras que la aproximación
hacia atrás da \(0.6885\). También es un error de truncamiento mayor porque la
fórmula usada en el extremo es de primer orden.

---

## Inciso e) Trapecio compuesto con la tabla del inciso d

La tabla tiene 11 puntos, por lo tanto hay \(n=10\) subintervalos y

\[
h=0.1.
\]

### Fórmula del trapecio compuesto

\[
I\approx \frac{h}{2}
\left[
y_0+2(y_1+y_2+\cdots+y_9)+y_{10}
\right].
\]

Los valores \(y_i\) son los de la tabla redondeada del inciso d.

### Suma de puntos interiores

\[
\begin{aligned}
S_{\text{int}}
={}&0.0990+0.1952+0.2861+0.3700+0.4458\\
&+0.5131+0.5722+0.6235+0.6679
\end{aligned}
\]

\[
\boxed{S_{\text{int}}=3.7728}
\]

### Sustitución

\[
I_T=\frac{0.1}{2}
\left[
0.0499+2(3.7728)+0.6885
\right]
\]

\[
I_T=0.05\left[0.0499+7.5456+0.6885\right]
\]

\[
I_T=0.05(8.2840)
\]

\[
\boxed{I_T=0.4142000000}
\]

### Valor exacto de la integral

Como el integrando es \(f'(x)\), por el teorema fundamental del cálculo:

\[
\int_0^1 f'(x)\,dx=f(1)-f(0).
\]

\[
f(1)=\sqrt{1+1^2}=\sqrt2=1.4142135624,
\]

\[
f(0)=\sqrt{1+0^2}=1.
\]

\[
I_{\text{exacta}}=\sqrt2-1=0.4142135624.
\]

### Error absoluto exacto

\[
E_T=\left|0.4142000000-0.4142135624\right|
\]

\[
\boxed{E_T=0.0000135624}
\]

Con 10 cifras decimales:

\[
\boxed{I_T\approx 0.4142000000}
\]

\[
\boxed{E_T\approx 0.0000135624}
\]

### Observación importante

Si en lugar de la tabla redondeada se usaran los valores exactos de \(f'\), el
trapecio compuesto telescoparía:

\[
\frac{h}{2}\left[f'(x_0)+2\sum_{i=1}^{9}f'(x_i)+f'(x_{10})\right]
=f(1)-f(0).
\]

Por eso el error \(0.0000135624\) se debe principalmente al redondeo a 4 cifras
decimales de la tabla, no al método del trapecio en sí.

---

## Resumen de resultados

| Cantidad | Resultado |
|---|---:|
| \(f'(1)\) exacta | 0.7071067812 |
| Centrada, \(h=0.5\) | 0.6847416490 |
| Error con \(h=0.5\) | 0.0223651322 |
| Centrada, \(h=0.01\) | 0.7070979423 |
| Error con \(h=0.01\) | 0.0000088389 |
| Trapecio con tabla \(h=0.1\) | 0.4142000000 |
| Integral exacta \(\sqrt2-1\) | 0.4142135624 |
| Error del trapecio | 0.0000135624 |

---

## Salida del programa

El archivo `salida_problema2.txt` contiene la salida completa. Los resultados
principales son:

```text
f'(1) = 1/sqrt(2) = 0.7071067812

f(1.5) = 1.8027756377
f(0.5) = 1.1180339887
f'(1) aproximada = 0.6847416490
f'(1) aproximada (6 decimales) = 0.684742
Error absoluto exacto = 0.0223651322

f(1.01) = 1.4213022198
f(0.99) = 1.4071602610
f'(1) aproximada = 0.7070979423
f'(1) aproximada (6 decimales) = 0.707098
Error absoluto exacto = 0.0000088389
Cociente de errores E(0.5)/E(0.01) = 2530.310052
Cociente esperado (0.5/0.01)^2 = 2500.000000

Tabla con h = 0.1:
    0       0.0000           0.0499     0.0000000000
    1       0.1000           0.0990     0.0995037190
    2       0.2000           0.1952     0.1961161351
    3       0.3000           0.2861     0.2873478856
    4       0.4000           0.3700     0.3713906764
    5       0.5000           0.4458     0.4472135955
    6       0.6000           0.5131     0.5144957554
    7       0.7000           0.5722     0.5734623444
    8       0.8000           0.6235     0.6246950476
    9       0.9000           0.6679     0.6689647316
   10       1.0000           0.6885     0.7071067812

Suma de puntos internos = 3.7728000000
I_Trapecio = 0.4142000000
I_exacta   = f(1) - f(0) = sqrt(2) - 1 = 0.4142135624
Error absoluto exacto = 0.0000135624
```
