# Problema 1 — Segundo Parcial 2025 (Integración numérica y spline)

Resolución completa del problema, con el desarrollo manual paso a paso y la
verificación con el programa `problema1.c++`.

## Archivos de esta carpeta

| Archivo | Qué es |
|---|---|
| `problema1.c++` | Programa: Gauss-Legendre de 2 puntos, Simpson 1/3 simple, spline cúbico natural, tabla de 5 puntos y Simpson compuesto. |
| `problema1` | Ejecutable ya compilado. |
| `salida_problema1.txt` | Salida completa del programa, útil para la captura de consola. |
| `captura_problema1.png` | Captura tipo consola de la salida del programa. |
| `README - Problema 1.md` | Este documento. |

Compilar y ejecutar:

```bash
g++ -Wall -Wextra -O2 -o problema1 problema1.c++ -lm
./problema1
```

---

## Enunciado

Se desea calcular la integral

\[
I=\int_0^1 e^x\sin(x)\,dx.
\]

El valor exacto dado por el examen es

\[
I_{\text{exacta}}=0.9093306736.
\]

**a)** Resolver a mano, usando **cuadratura de Gauss-Legendre de dos puntos**.
Escribir el resultado con 6 cifras decimales y calcular el error absoluto
exacto.

**b)** Repetir el inciso a) con **Simpson 1/3 simple**.

**c)** En la tabla se muestran 7 pares de datos que corresponden a
\(f(x)=e^x\sin(x)\). Construir una nueva tabla de **5 datos equiespaciados**
usando **spline cúbico natural**. Conservar 4 cifras decimales.

Tabla dada:

| \(x\) | 0.00 | 0.20 | 0.43 | 0.55 | 0.85 | 0.90 | 1.00 |
|---|---:|---:|---:|---:|---:|---:|---:|
| \(f(x)\) | 0.0000 | 0.2426 | 0.6408 | 0.9059 | 1.7577 | 1.9267 | 2.2873 |

**d)** A partir de la tabla del inciso c), calcular la integral con **Simpson
compuesto**. Escribir el valor con 10 cifras decimales y calcular el error
absoluto exacto.

---

## Función del problema

La tabla corresponde a

\[
f(x)=e^x\sin(x).
\]

Por ejemplo:

\[
f(0.20)=e^{0.20}\sin(0.20)=0.2426,
\]

\[
f(1.00)=e^{1}\sin(1)=2.2873.
\]

La integral pedida es, entonces,

\[
I=\int_0^1 e^x\sin(x)\,dx.
\]

---

## Inciso a) Gauss-Legendre de 2 puntos (a mano)

### Fórmula

Para dos puntos, en un intervalo \([a,b]\):

\[
\int_a^b f(x)\,dx\approx \frac{b-a}{2}\left[f(x_1)+f(x_2)\right],
\]

donde

\[
x_1=\frac{a+b}{2}-\frac{b-a}{2}\cdot\frac{1}{\sqrt3},
\qquad
x_2=\frac{a+b}{2}+\frac{b-a}{2}\cdot\frac{1}{\sqrt3}.
\]

Como \(a=0\) y \(b=1\):

\[
\frac{a+b}{2}=0.5,\qquad \frac{b-a}{2}=0.5.
\]

### Nodos

\[
x_1=0.5-0.5\cdot\frac{1}{\sqrt3}=0.2113248654,
\]

\[
x_2=0.5+0.5\cdot\frac{1}{\sqrt3}=0.7886751346.
\]

### Evaluación de la función

\[
f(x_1)=e^{0.2113248654}\sin(0.2113248654)=0.2591137921,
\]

\[
f(x_2)=e^{0.7886751346}\sin(0.7886751346)=1.5610642502.
\]

### Resultado

\[
I_{GL}=\frac{1-0}{2}\left(0.2591137921+1.5610642502\right)
\]

\[
I_{GL}=0.5\cdot 1.8201780423
\]

\[
\boxed{I_{GL}=0.9100890212}
\]

Con 6 cifras decimales:

\[
\boxed{I_{GL}\approx 0.910089}
\]

### Error absoluto exacto

\[
E_{GL}=\left|0.9100890212-0.9093306736\right|
\]

\[
\boxed{E_{GL}=0.0007583476}
\]

Con 6 cifras decimales:

\[
\boxed{E_{GL}\approx 0.000758}
\]

---

## Inciso b) Simpson 1/3 simple (a mano)

### Fórmula

Para un solo panel, con \(h=(b-a)/2\):

\[
\int_a^b f(x)\,dx\approx \frac{b-a}{6}\left[f(a)+4f\left(\frac{a+b}{2}\right)+f(b)\right].
\]

Con \(a=0\), \(b=1\):

\[
I_{S1/3}=\frac{1}{6}\left[f(0)+4f(0.5)+f(1)\right].
\]

### Evaluaciones

\[
f(0)=e^0\sin(0)=0,
\]

\[
f(0.5)=e^{0.5}\sin(0.5)=0.7904390832,
\]

\[
f(1)=e^1\sin(1)=2.2873552872.
\]

### Resultado

\[
I_{S1/3}=\frac{1}{6}\left[0+4(0.7904390832)+2.2873552872\right]
\]

\[
I_{S1/3}=\frac{1}{6}\left(3.1617563329+2.2873552872\right)
\]

\[
I_{S1/3}=\frac{5.4491116200}{6}
\]

\[
\boxed{I_{S1/3}=0.9081852700}
\]

Con 6 cifras decimales:

\[
\boxed{I_{S1/3}\approx 0.908185}
\]

### Error absoluto exacto

\[
E_{S1/3}=\left|0.9081852700-0.9093306736\right|
\]

\[
\boxed{E_{S1/3}=0.0011454036}
\]

Con 6 cifras decimales:

\[
\boxed{E_{S1/3}\approx 0.001145}
\]

---

## Inciso c) Spline cúbico natural: nueva tabla de 5 puntos

Se quieren 5 datos equiespaciados en \([0,1]\):

\[
x=0.00,\ 0.25,\ 0.50,\ 0.75,\ 1.00.
\]

Los extremos ya están en la tabla:

\[
f(0.00)=0.0000,\qquad f(1.00)=2.2873.
\]

Los valores interiores se obtienen evaluando el spline cúbico natural que pasa
por los 7 puntos dados.

### Idea del spline cúbico natural

En cada intervalo \([X_i,X_{i+1}]\), el spline es un polinomio cúbico

\[
S_i(x)=a_ix^3+b_ix^2+c_ix+d_i.
\]

Se imponen:

1. que cada cúbica pase por los dos extremos del intervalo;
2. continuidad de la primera derivada en los nodos interiores;
3. continuidad de la segunda derivada en los nodos interiores;
4. condición natural: \(S''(X_0)=0\) y \(S''(X_6)=0\).

Una forma cómoda de construirlo es calcular primero las segundas derivadas
\(M_i=S''(X_i)\). Para los nodos interiores se plantea:

\[
h_{i-1}M_{i-1}+2(h_{i-1}+h_i)M_i+h_iM_{i+1}
=
6\left[
\frac{Y_{i+1}-Y_i}{h_i}
-
\frac{Y_i-Y_{i-1}}{h_{i-1}}
\right],
\]

con \(M_0=M_6=0\).

### Segundas derivadas obtenidas

\[
\begin{aligned}
M_0&=0.0000000000,\\
M_1&=2.9162824378,\\
M_2&=2.6166225673,\\
M_3&=3.0399429605,\\
M_4&=3.0448440172,\\
M_5&=4.0125259971,\\
M_6&=0.0000000000.
\end{aligned}
\]

### Coeficientes de los splines

Resolviendo el sistema se obtienen los siguientes polinomios por tramo:

\[
S_1(x)=2.4302353648x^3+0.0000000000x^2+1.1157905854x+0.0000000000,
\]

\[
S_2(x)=-0.2171448337x^3+1.5884281191x^2+0.7981049616x+0.0211790416,
\]

\[
S_3(x)=0.5879449905x^3+0.5498622459x^2+1.2446882871x-0.0428312351,
\]

\[
S_4(x)=0.0027228093x^3+1.5154788449x^2+0.7135991576x+0.0545351053,
\]

\[
S_5(x)=3.2256066000x^3-6.7028748213x^2+7.6991997739x-1.9247184026,
\]

\[
S_6(x)=-6.6875433286x^3+20.0626299857x^2-16.3897545524x+5.3019678953.
\]

### Evaluación en los puntos pedidos

**Para \(x=0.25\)** se usa \(S_2\):

\[
S_2(0.25)=0.3165891514\approx 0.3166.
\]

**Para \(x=0.50\)** se usa \(S_3\):

\[
S_3(0.50)=0.7904715938\approx 0.7905.
\]

**Para \(x=0.75\)** se usa \(S_4\):

\[
S_4(0.75)=1.4433400090\approx 1.4433.
\]

### Tabla nueva pedida (4 cifras decimales)

| \(x\) | 0.00 | 0.25 | 0.50 | 0.75 | 1.00 |
|---|---:|---:|---:|---:|---:|
| \(y\) | 0.0000 | 0.3166 | 0.7905 | 1.4433 | 2.2873 |

\[
\boxed{
\begin{array}{c|ccccc}
x & 0.00 & 0.25 & 0.50 & 0.75 & 1.00\\
\hline
y & 0.0000 & 0.3166 & 0.7905 & 1.4433 & 2.2873
\end{array}}
\]

---

## Inciso d) Simpson compuesto con la tabla del inciso c

La tabla tiene 5 puntos, por lo tanto hay \(n=4\) subintervalos y

\[
h=\frac{1.00-0.00}{4}=0.25.
\]

### Fórmula de Simpson compuesto

Para \(n\) par:

\[
I\approx \frac{h}{3}\left[
y_0+4y_1+2y_2+4y_3+y_4
\right].
\]

Sustituyendo la tabla redondeada del inciso c:

\[
I_{SC}=\frac{0.25}{3}
\left[
0.0000+4(0.3166)+2(0.7905)+4(1.4433)+2.2873
\right].
\]

Calculando el corchete:

\[
4(0.3166)=1.2664,
\]

\[
2(0.7905)=1.5810,
\]

\[
4(1.4433)=5.7732.
\]

Entonces:

\[
I_{SC}=\frac{0.25}{3}
\left(0.0000+1.2664+1.5810+5.7732+2.2873\right)
\]

\[
I_{SC}=\frac{0.25}{3}(10.9079)
\]

\[
\boxed{I_{SC}=0.9089916667}
\]

Con 10 cifras decimales:

\[
\boxed{I_{SC}\approx 0.9089916667}
\]

### Error absoluto exacto

\[
E_{SC}=\left|0.9089916667-0.9093306736\right|
\]

\[
\boxed{E_{SC}=0.0003390069}
\]

---

## Resumen de resultados

| Método | Resultado | Error absoluto exacto |
|---|---:|---:|
| Gauss-Legendre 2 puntos | 0.910089 | 0.000758 |
| Simpson 1/3 simple | 0.908185 | 0.001145 |
| Simpson compuesto (tabla spline) | 0.9089916667 | 0.0003390069 |
| Valor exacto dado | 0.9093306736 | — |

El método de Simpson compuesto con la tabla del spline es el que más se acerca
al valor exacto, como era de esperar al usar más puntos.

---

## Salida del programa

El archivo `salida_problema1.txt` contiene la salida completa. Los resultados
principales son:

```text
x1 = 0.2113248654
x2 = 0.7886751346
f(x1) = 0.2591137921
f(x2) = 1.5610642502
I_Gauss = 0.9100890212
I_Gauss (6 decimales) = 0.910089
Error absoluto exacto = 0.0007583476

f(0.0) = 0.0000000000
f(0.5) = 0.7904390832
f(1.0) = 2.2873552872
I_Simpson_simple = 0.9081852700
I_Simpson_simple (6 decimales) = 0.908185
Error absoluto exacto = 0.0011454036

Tabla nueva de 5 datos equiespaciados (4 decimales):
    0       0.0000       0.0000
    1       0.2500       0.3166
    2       0.5000       0.7905
    3       0.7500       1.4433
    4       1.0000       2.2873

I_Simpson_compuesto = 0.9089916667
Error absoluto exacto = 0.0003390069
```
