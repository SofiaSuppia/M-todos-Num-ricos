# Problema 1 — Segundo Parcial 2024 (Spline, derivación e integración)

Resolución completa del problema, con el desarrollo manual paso a paso y la
verificación con el programa `problema1.c++`.

## Archivos de esta carpeta

| Archivo | Qué es |
|---|---|
| `problema1.c++` | Programa: spline cúbico natural, tabla equiespaciada, derivadas por diferencias finitas, trapecio, Gauss de 2 y 6 puntos. |
| `problema1` | Ejecutable ya compilado. |
| `salida_problema1.txt` | Salida completa del programa, útil para la captura de consola. |
| `captura_problema1.png` | Captura tipo consola de la salida del programa. |
| `README - Problema 1.md` | Este documento. |

Compilar y ejecutar:

```bash
g++ -std=c++17 -Wall -Wextra -Wpedantic -O2 problema1.c++ -o problema1 -lm
./problema1
```

---

## Enunciado

En la tabla se muestran una serie de datos que provienen de la función
\(f(x)=\ln(1+x^2)-\sin(x)\).

| \(x\) | 1.00 | 1.20 | 1.45 | 1.78 | 2.00 |
|---|---:|---:|---:|---:|---:|
| \(f(x)\) | -0.1483 | -0.0400 | 0.1395 | 0.4493 | 0.7001 |

Mediante interpolación **spline cúbica natural**:

**a)** Escribir en papel el polinomio correspondiente al subintervalo entre
\(1.20\) y \(1.45\).

**b)** Confeccionar una tabla con datos equiespaciados, conservando el número
de puntos, donde figuren \(x\), el valor de la función exacta, el valor
interpolado y el error absoluto exacto.

**c)** A partir de la tabla del inciso b), calcular la derivada en los puntos
equiespaciados usando operadores en diferencias finitas: hacia adelante/atrás
en los extremos, según corresponda, y centradas en los puntos interiores.
Confeccionar una segunda tabla con los valores de las derivadas obtenidas
usando los datos de las aproximaciones spline. Usar los operadores de orden
más bajo en \(h=0.25\).

**d)** Con los datos de la tabla del inciso c), calcular

\[
I=\int_1^2 \frac{df}{dx}\,dx
\]

de forma exacta y mediante el método del trapecio. Indicar el error absoluto
exacto.

**e)** Calcular la integral anterior, a mano y en papel, usando cuadratura de
Gauss de dos puntos y obtener el error porcentual exacto.

**f)** Calcular, mediante el código correspondiente, la integral usando
cuadratura de Gauss de 6 puntos, escribir en papel el resultado y adjuntar la
captura completa de la consola.

---

## Función exacta y derivada exacta

La función es

\[
f(x)=\ln(1+x^2)-\sin(x).
\]

Su derivada exacta, necesaria para comparar y para Gauss, es

\[
\boxed{f'(x)=\frac{2x}{1+x^2}-\cos(x)}.
\]

Una integral importante del problema es

\[
I=\int_1^2 \frac{df}{dx}\,dx
=\int_1^2 f'(x)\,dx
=f(2)-f(1).
\]

Por el teorema fundamental del cálculo:

\[
\boxed{I=f(2)-f(1)}.
\]

Numéricamente:

\[
f(1)=\ln(2)-\sin(1)=-0.1483238042,
\]

\[
f(2)=\ln(5)-\sin(2)=0.7001404856.
\]

Por lo tanto,

\[
I=0.7001404856-(-0.1483238042)
=0.8484642899.
\]

\[
\boxed{I_{\text{exacta}}=0.8484642899}
\]

---

## Inciso a) Spline cúbico natural y polinomio en [1.20, 1.45]

### Datos y espaciados

\[
X=[1.00,\ 1.20,\ 1.45,\ 1.78,\ 2.00],
\]

\[
Y=[-0.1483,\ -0.0400,\ 0.1395,\ 0.4493,\ 0.7001].
\]

Los espaciados son:

\[
h_0=1.20-1.00=0.20,
\]

\[
h_1=1.45-1.20=0.25,
\]

\[
h_2=1.78-1.45=0.33,
\]

\[
h_3=2.00-1.78=0.22.
\]

### Sistema para las segundas derivadas

Para un spline cúbico natural se definen

\[
M_i=S''(X_i),
\]

con las condiciones naturales

\[
M_0=0,\qquad M_4=0.
\]

Para los nodos interiores \(i=1,2,3\) se plantea:

\[
h_{i-1}M_{i-1}+2(h_{i-1}+h_i)M_i+h_iM_{i+1}
=
6\left[
\frac{Y_{i+1}-Y_i}{h_i}
-
\frac{Y_i-Y_{i-1}}{h_{i-1}}
\right].
\]

Entonces el sistema es:

\[
\begin{bmatrix}
2(h_0+h_1) & h_1 & 0\\
h_1 & 2(h_1+h_2) & h_2\\
0 & h_2 & 2(h_2+h_3)
\end{bmatrix}
\begin{bmatrix}
M_1\\M_2\\M_3
\end{bmatrix}
=
\begin{bmatrix}
6\left(\frac{Y_2-Y_1}{h_1}-\frac{Y_1-Y_0}{h_0}\right)\\
6\left(\frac{Y_3-Y_2}{h_2}-\frac{Y_2-Y_1}{h_1}\right)\\
6\left(\frac{Y_4-Y_3}{h_3}-\frac{Y_3-Y_2}{h_2}\right)
\end{bmatrix}.
\]

Reemplazando los valores:

\[
\frac{Y_1-Y_0}{h_0}=\frac{-0.0400+0.1483}{0.20}=0.5415,
\]

\[
\frac{Y_2-Y_1}{h_1}=\frac{0.1395+0.0400}{0.25}=0.7180,
\]

\[
\frac{Y_3-Y_2}{h_2}=\frac{0.4493-0.1395}{0.33}=0.9387878788,
\]

\[
\frac{Y_4-Y_3}{h_3}=\frac{0.7001-0.4493}{0.22}=1.1400000000.
\]

Por lo tanto:

\[
6(0.7180-0.5415)=1.0590,
\]

\[
6(0.9387878788-0.7180)=1.3247272727,
\]

\[
6(1.1400000000-0.9387878788)=1.2072727273.
\]

El sistema queda:

\[
\begin{bmatrix}
0.90 & 0.25 & 0\\
0.25 & 1.16 & 0.33\\
0 & 0.33 & 1.10
\end{bmatrix}
\begin{bmatrix}
M_1\\M_2\\M_3
\end{bmatrix}
=
\begin{bmatrix}
1.0590\\
1.3247272727\\
1.2072727273
\end{bmatrix}.
\]

Resolviendo:

\[
\boxed{M_1=0.98942474}
\]

\[
\boxed{M_2=0.67407094}
\]

\[
\boxed{M_3=0.89529938}
\]

Además,

\[
M_0=0,\qquad M_4=0.
\]

### Polinomio del tramo [1.20, 1.45]

Para el intervalo \([X_i,X_{i+1}]\), el spline se puede escribir como

\[
S_i(x)=
\frac{X_{i+1}-x}{h_i}Y_i
+
\frac{x-X_i}{h_i}Y_{i+1}
+
\left[
\left(\frac{X_{i+1}-x}{h_i}\right)^3
-
\frac{X_{i+1}-x}{h_i}
\right]\frac{h_i^2}{6}M_i
\]

\[
+
\left[
\left(\frac{x-X_i}{h_i}\right)^3
-
\frac{x-X_i}{h_i}
\right]\frac{h_i^2}{6}M_{i+1}.
\]

Para el tramo pedido, \(i=1\), \(X_1=1.20\), \(X_2=1.45\),
\(h_1=0.25\), \(Y_1=-0.0400\), \(Y_2=0.1395\),
\(M_1=0.98942474\), \(M_2=0.67407094\).

Conviene usar la variable local

\[
t=x-1.20.
\]

Entonces el polinomio se puede escribir como

\[
S_1(t)=\frac{M_2-M_1}{6h_1}t^3
+\frac{M_1}{2}t^2
+\left[
\frac{Y_2-Y_1}{h_1}
-\frac{h_1(2M_1+M_2)}{6}
\right]t
+Y_1.
\]

Calculando cada coeficiente:

\[
\frac{M_2-M_1}{6h_1}
=
\frac{0.67407094-0.98942474}{6(0.25)}
=-0.2102358638,
\]

\[
\frac{M_1}{2}=0.4947123691,
\]

\[
\frac{Y_2-Y_1}{h_1}=\frac{0.1395+0.0400}{0.25}=0.7180,
\]

\[
\frac{h_1(2M_1+M_2)}{6}
=
\frac{0.25(2(0.98942474)+0.67407094)}{6}
=0.1105383508,
\]

\[
0.7180-0.1105383508=0.6074616492.
\]

Por lo tanto:

\[
\boxed{
S_1(t)=-0.2102358638\,t^3
+0.4947123691\,t^2
+0.6074616492\,t
-0.0400000000
}
\]

con

\[
\boxed{t=x-1.20}.
\]

Si se prefiere el polinomio en potencias de \(x\), se desarrolla:

\[
S_1(x)=-0.2102358638x^3
+1.2515614787x^2
-1.4880669682x
+0.3067194051.
\]

Verificación:

\[
S_1(1.20)=-0.0400000000,
\]

\[
S_1(1.45)=0.1395000000,
\]

que coinciden con los datos de la tabla.

---

## Inciso b) Tabla de 5 puntos equiespaciados

El intervalo es \([1,2]\) y se conservan 5 puntos. Por lo tanto:

\[
h=\frac{2-1}{4}=0.25.
\]

Los puntos son:

\[
x=1.00,\ 1.25,\ 1.50,\ 1.75,\ 2.00.
\]

Para cada punto se calcula:

- el valor exacto \(f(x)\);
- el valor interpolado \(S(x)\);
- el error absoluto exacto \(|S(x)-f(x)|\).

En la tabla se escribe el spline redondeado a 4 cifras decimales, porque los
datos originales tienen esa precisión. El error se calcula con el valor que
efectivamente se escribe en la tabla.

| \(i\) | \(x\) | \(f(x)\) exacta | \(S(x)\) spline | Error absoluto |
|---:|---:|---:|---:|---:|
| 0 | 1.00 | -0.14832380 | -0.1483 | \(2.38042480\times10^{-5}\) |
| 1 | 1.25 | -0.00800127 | -0.0084 | \(3.98725109\times10^{-4}\) |
| 2 | 1.50 | 0.18116001 | 0.1811 | \(6.00097376\times10^{-5}\) |
| 3 | 1.75 | 0.41781260 | 0.4175 | \(3.12600782\times10^{-4}\) |
| 4 | 2.00 | 0.70014049 | 0.7001 | \(4.04856084\times10^{-5}\) |

### Ejemplo de evaluación en x = 1.25

El punto \(x=1.25\) está en el tramo \([1.20,1.45]\), entonces

\[
t=1.25-1.20=0.05.
\]

\[
S(1.25)=
-0.2102358638(0.05)^3
+0.4947123691(0.05)^2
+0.6074616492(0.05)
-0.0400000000.
\]

\[
S(1.25)=-0.0000262795+0.0012367809+0.0303730825-0.0400000000
\]

\[
S(1.25)=-0.0084164161\approx -0.0084.
\]

El valor exacto es

\[
f(1.25)=\ln(1+1.25^2)-\sin(1.25)
=-0.0080012749.
\]

El error absoluto exacto, usando el valor redondeado de la tabla, es

\[
| -0.0084-(-0.0080012749)|
=3.98725109\times10^{-4}.
\]

**Observación:** si se usara el valor sin redondear \(S(1.25)=-0.0084164161\),
el error sería \(4.151430\times10^{-4}\). En la tabla se informa el error del
valor redondeado, que es el que se escribe y se usa en el inciso c).

---

## Inciso c) Tabla de derivadas por diferencias finitas

Se usan los valores redondeados del spline de la tabla anterior:

\[
y_0=-0.1483,\quad
y_1=-0.0084,\quad
y_2=0.1811,\quad
y_3=0.4175,\quad
y_4=0.7001.
\]

El espaciado es

\[
h=0.25.
\]

### Fórmulas usadas

En el primer punto se usa diferencia hacia adelante:

\[
f'(x_0)\approx \frac{y_1-y_0}{h}.
\]

En los puntos interiores se usa diferencia centrada:

\[
f'(x_i)\approx \frac{y_{i+1}-y_{i-1}}{2h},
\qquad i=1,2,3.
\]

En el último punto se usa diferencia hacia atrás:

\[
f'(x_4)\approx \frac{y_4-y_3}{h}.
\]

### Cálculos

\[
f'(1.00)\approx \frac{-0.0084-(-0.1483)}{0.25}
=\frac{0.1399}{0.25}
=0.5596.
\]

\[
f'(1.25)\approx \frac{0.1811-(-0.1483)}{2(0.25)}
=\frac{0.3294}{0.50}
=0.6588.
\]

\[
f'(1.50)\approx \frac{0.4175-(-0.0084)}{2(0.25)}
=\frac{0.4259}{0.50}
=0.8518.
\]

\[
f'(1.75)\approx \frac{0.7001-0.1811}{2(0.25)}
=\frac{0.5190}{0.50}
=1.0380.
\]

\[
f'(2.00)\approx \frac{0.7001-0.4175}{0.25}
=\frac{0.2826}{0.25}
=1.1304.
\]

### Tabla de derivadas

| \(i\) | \(x\) | \(S(x)\) tabla | \(f'(x)\) aproximada | \(f'(x)\) exacta |
|---:|---:|---:|---:|---:|
| 0 | 1.00 | -0.1483 | 0.5596 | 0.45969769 |
| 1 | 1.25 | -0.0084 | 0.6588 | 0.66028739 |
| 2 | 1.50 | 0.1811 | 0.8518 | 0.85233972 |
| 3 | 1.75 | 0.4175 | 1.0380 | 1.03978452 |
| 4 | 2.00 | 0.7001 | 1.1304 | 1.21614684 |

La columna \(f'(x)\) exacta se incluye solo para comparar; el método pedido es
la diferencia finita aplicada a los datos del spline.

---

## Inciso d) Integral por trapecio compuesto

La integral pedida es

\[
I=\int_1^2 f'(x)\,dx=f(2)-f(1).
\]

Ya se calculó:

\[
\boxed{I_{\text{exacta}}=0.8484642899}.
\]

Con los datos de la tabla del inciso c), el trapecio compuesto con \(h=0.25\)
es:

\[
I_T=\frac{h}{2}
\left[
f'_0+2f'_1+2f'_2+2f'_3+f'_4
\right].
\]

Reemplazando:

\[
I_T=\frac{0.25}{2}
\left[
0.5596+2(0.6588)+2(0.8518)+2(1.0380)+1.1304
\right].
\]

\[
I_T=0.125
\left[
0.5596+1.3176+1.7036+2.0760+1.1304
\right].
\]

\[
I_T=0.125(6.7872)=0.8484.
\]

\[
\boxed{I_T=0.8484000000}
\]

El error absoluto exacto es:

\[
E=|I_T-I_{\text{exacta}}|
=|0.8484-0.8484642899|.
\]

\[
\boxed{E=6.4289856370\times10^{-5}}
\]

---

## Inciso e) Gauss-Legendre de 2 puntos

Se integra

\[
I=\int_1^2 f'(x)\,dx.
\]

La cuadratura de Gauss-Legendre de 2 puntos en \([-1,1]\) usa:

\[
t_1=-\frac{1}{\sqrt3},\qquad t_2=\frac{1}{\sqrt3},
\]

\[
w_1=w_2=1.
\]

Para pasar al intervalo \([1,2]\) se usa la transformación:

\[
x=\frac{a+b}{2}+\frac{b-a}{2}t.
\]

Como \(a=1\), \(b=2\):

\[
\frac{a+b}{2}=1.5,\qquad \frac{b-a}{2}=0.5.
\]

Entonces:

\[
x=1.5+0.5t.
\]

Los nodos son:

\[
x_1=1.5+0.5\left(-\frac{1}{\sqrt3}\right)
=1.2113248654,
\]

\[
x_2=1.5+0.5\left(\frac{1}{\sqrt3}\right)
=1.7886751346.
\]

La fórmula transformada es:

\[
I\approx \frac{b-a}{2}\left[w_1 f'(x_1)+w_2 f'(x_2)\right].
\]

Como \(b-a=1\):

\[
I\approx 0.5\left[f'(1.2113248654)+f'(1.7886751346)\right].
\]

Evaluando:

\[
f'(1.2113248654)
=
\frac{2(1.2113248654)}{1+1.2113248654^2}
-\cos(1.2113248654)
=0.6301205038,
\]

\[
f'(1.7886751346)
=
\frac{2(1.7886751346)}{1+1.7886751346^2}
-\cos(1.7886751346)
=1.0680391982.
\]

Entonces:

\[
I_{G2}=0.5(0.6301205038+1.0680391982)
\]

\[
I_{G2}=0.5(1.6981597020)
\]

\[
\boxed{I_{G2}=0.8490798510}
\]

Error absoluto exacto:

\[
E=|0.8490798510-0.8484642899|
=6.1556113827\times10^{-4}.
\]

Error porcentual exacto:

\[
E_{\%}=
\frac{|I_{G2}-I_{\text{exacta}}|}{|I_{\text{exacta}}|}\times100
\]

\[
E_{\%}=
\frac{0.0006155611}{0.8484642899}\times100
=0.0725500349\%.
\]

\[
\boxed{E_{\%}=0.0725500349\%}
\]

---

## Inciso f) Gauss-Legendre de 6 puntos

Para 6 puntos se usan los nodos y pesos:

| \(i\) | \(t_i\) | \(w_i\) |
|---:|---:|---:|
| 0 | -0.932469514 | 0.1713245 |
| 1 | -0.661209386 | 0.3607616 |
| 2 | -0.238619186 | 0.4679139 |
| 3 | 0.238619186 | 0.4679139 |
| 4 | 0.661209386 | 0.3607616 |
| 5 | 0.932469514 | 0.1713245 |

Con la misma transformación \(x=1.5+0.5t\), los nodos en \([1,2]\) y los
valores de \(f'\) son:

| \(i\) | \(x_i\) | \(f'(x_i)\) | \(w_i\) |
|---:|---:|---:|---:|
| 0 | 1.033765243 | 0.4878616129 | 0.1713245 |
| 1 | 1.169395307 | 0.5971712578 | 0.3607616 |
| 2 | 1.380690407 | 0.7611713226 | 0.4679139 |
| 3 | 1.619309593 | 0.9426060408 | 0.4679139 |
| 4 | 1.830604693 | 1.0983372864 | 0.3607616 |
| 5 | 1.966234757 | 1.1933520681 | 0.1713245 |

La integral aproximada es:

\[
I_{G6}=0.5\sum_{i=0}^{5}w_i f'(x_i).
\]

Realizando la suma:

\[
I_{G6}=0.5(1.6969285793)
\]

\[
\boxed{I_{G6}=0.8484642897}
\]

Error absoluto exacto:

\[
E=|0.8484642897-0.8484642899|
=1.9071944024\times10^{-10}.
\]

Error porcentual exacto:

\[
E_{\%}=
\frac{1.9071944024\times10^{-10}}{0.8484642899}\times100
=2.2478192956\times10^{-8}\%.
\]

\[
\boxed{E=1.9071944024\times10^{-10}}
\]

---

## Resumen comparativo

| Método | Valor de la integral | Error absoluto exacto |
|---|---:|---:|
| Exacta \(f(2)-f(1)\) | 0.8484642899 | 0 |
| Trapecio compuesto | 0.8484000000 | \(6.4289856370\times10^{-5}\) |
| Gauss-Legendre 2 puntos | 0.8490798510 | \(6.1556113827\times10^{-4}\) |
| Gauss-Legendre 6 puntos | 0.8484642897 | \(1.9071944024\times10^{-10}\) |

Conclusión: el trapecio compuesto usando la tabla de derivadas da una
aproximación muy razonable, pero Gauss de 6 puntos es prácticamente exacto para
esta integral, porque la función derivada es suave y la cuadratura de orden
superior captura mucho mejor su forma.
