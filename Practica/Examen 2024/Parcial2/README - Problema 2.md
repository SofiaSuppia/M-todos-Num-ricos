# Problema 2 — Segundo Parcial 2024 (EDO: AB2 y RK4)

Resolución completa del problema, con el desarrollo manual paso a paso y la
verificación con el programa `problema2.c++`.

## Archivos de esta carpeta

| Archivo | Qué es |
|---|---|
| `problema2.c++` | Programa: Euler para iniciar AB2, Adams-Bashforth de 2 pasos, RK4 y escritura del CSV. |
| `problema2` | Ejecutable ya compilado. |
| `salida_problema2.txt` | Salida completa del programa, útil para la captura de consola. |
| `captura_problema2.png` | Captura tipo consola de la salida del programa. |
| `resultados_problema2.csv` | Datos de \(x\), solución exacta, AB2, RK4 y errores. |
| `grafico_problema2.py` | Script que genera el gráfico de las tres soluciones. |
| `grafico_problema2.png` | Gráfico de la solución exacta, AB2 y RK4. |
| `README - Problema 2.md` | Este documento. |

Compilar y ejecutar:

```bash
g++ -std=c++17 -Wall -Wextra -Wpedantic -O2 problema2.c++ -o problema2 -lm
./problema2
python3 grafico_problema2.py
```

---

## Enunciado

Se propone un método de resolución de ecuaciones diferenciales ordinarias de la
forma

\[
y'=f(x,y),
\]

con dato inicial \(y_0\), en el que cada nuevo valor se obtiene a partir de los
dos pasos anteriores, dado la siguiente ecuación:

\[
y_{n+1}=y_n+\frac{h}{2}
\left[
3f(x_n,y_n)-f(x_{n-1},y_{n-1})
\right].
\]

Si \(n\geq 1\). Para \(n=0\) se toma el valor inicial, y para \(n=1\) se podría
usar un método de un solo paso (como Euler) para calcular \(y_1\).

Considere el siguiente problema en el intervalo \(x\in[0,1]\):

\[
\frac{dy}{dx}=\frac{x e^{x^2}}{y},
\qquad y(0)=1.
\]

**a)** Escribir en papel el pseudo-código correspondiente a este método e
implementarlo en un programa.

**b)** Resolver la EDO dada de forma exacta.

**c)** Obtener la solución aproximada usando el método propuesto, con
\(h=0.1\), y escribir en papel el valor de \(y(0.5)\) y de \(y(1)\). Calcular
el error absoluto para ambos.

**d)** Repetir el inciso b), empleando el método de RK4.

**e)** Graficar, en un mismo plot, las tres soluciones (exacta, método
propuesto y RK4). Adjuntar el plot por aula virtual.

---

## Inciso b) Solución exacta de la EDO

La ecuación es separable:

\[
\frac{dy}{dx}=\frac{x e^{x^2}}{y}.
\]

Multiplicando por \(y\,dx\):

\[
y\,dy=x e^{x^2}\,dx.
\]

Integrando ambos miembros:

\[
\int y\,dy=\int x e^{x^2}\,dx.
\]

La integral de la izquierda es

\[
\int y\,dy=\frac{y^2}{2}.
\]

Para la integral de la derecha, se usa la sustitución

\[
u=x^2,\qquad du=2x\,dx.
\]

Entonces:

\[
\int x e^{x^2}\,dx
=
\frac12\int e^u\,du
=
\frac12 e^u
=
\frac12 e^{x^2}.
\]

Por lo tanto:

\[
\frac{y^2}{2}=\frac12 e^{x^2}+C.
\]

Multiplicando por 2:

\[
y^2=e^{x^2}+2C.
\]

Se usa la condición inicial \(y(0)=1\):

\[
1^2=e^{0}+2C=1+2C,
\]

\[
2C=0.
\]

Entonces:

\[
y^2=e^{x^2}.
\]

Como \(y(0)=1>0\), se toma la rama positiva:

\[
\boxed{y(x)=e^{x^2/2}}.
\]

La derivada de esta solución es

\[
y'(x)=e^{x^2/2}\cdot x
=\frac{x e^{x^2}}{y},
\]

que coincide con la EDO dada.

---

## Inciso a) Pseudo-código del método propuesto

El método es Adams-Bashforth de 2 pasos. El pseudo-código es:

```text
Entrada: x0, y0, xf, h, f(x,y)

n = (xf - x0) / h

x[0] = x0
y[0] = y0

// Para poder usar dos pasos anteriores, se calcula y[1] con Euler:
x[1] = x[0] + h
y[1] = y[0] + h * f(x[0], y[0])

// Desde i = 1 hasta n-1:
para i = 1, 2, ..., n-1:
    x[i+1] = x[i] + h
    y[i+1] = y[i] + (h/2) * (3*f(x[i], y[i]) - f(x[i-1], y[i-1]))
fin para

Salida: x[0..n], y[0..n]
```

En este problema:

\[
x_0=0,\qquad y_0=1,\qquad x_f=1,\qquad h=0.1,
\]

por lo que

\[
n=\frac{1-0}{0.1}=10.
\]

Se calculan 11 puntos:

\[
x=0.0,0.1,0.2,\ldots,1.0.
\]

La función es

\[
f(x,y)=\frac{x e^{x^2}}{y}.
\]

---

## Inciso c) Solución aproximada con AB2 y h = 0.1

### Inicialización con Euler

\[
y_1=y_0+h f(x_0,y_0).
\]

Como \(x_0=0\) y \(y_0=1\):

\[
f(0,1)=\frac{0\cdot e^{0}}{1}=0.
\]

Entonces:

\[
y_1=1+0.1(0)=1.
\]

### Segundo paso

\[
y_2=y_1+\frac{h}{2}
\left[
3f(x_1,y_1)-f(x_0,y_0)
\right].
\]

Con \(x_1=0.1\), \(y_1=1\):

\[
f(0.1,1)=\frac{0.1e^{0.01}}{1}=0.1010050167.
\]

Como \(f(0,1)=0\):

\[
y_2=1+\frac{0.1}{2}
\left[
3(0.1010050167)-0
\right].
\]

\[
y_2=1+0.05(0.3030150501)
\]

\[
y_2=1.0151507525.
\]

El programa continúa de la misma manera hasta \(x=1\).

### Tabla completa

| \(i\) | \(x\) | \(y\) exacta | \(y\) AB2 | Error AB2 | \(y\) RK4 | Error RK4 |
|---:|---:|---:|---:|---:|---:|---:|
| 0 | 0.0 | 1.0000000000 | 1.0000000000 | \(0\) | 1.0000000000 | \(0\) |
| 1 | 0.1 | 1.0050125209 | 1.0000000000 | \(5.012521\times10^{-3}\) | 1.0050125418 | \(2.091177\times10^{-8}\) |
| 2 | 0.2 | 1.0202013400 | 1.0151507525 | \(5.050588\times10^{-3}\) | 1.0202014258 | \(8.576757\times10^{-8}\) |
| 3 | 0.3 | 1.0460278599 | 1.0408588133 | \(5.169047\times10^{-3}\) | 1.0460280611 | \(2.012303\times10^{-7}\) |
| 4 | 0.4 | 1.0832870677 | 1.0779110587 | \(5.376009\times10^{-3}\) | 1.0832874471 | \(3.794512\times10^{-7}\) |
| 5 | 0.5 | 1.1331484531 | 1.1274641134 | \(5.684340\times10^{-3}\) | 1.1331490929 | \(6.397849\times10^{-7}\) |
| 6 | 0.6 | 1.1972173631 | 1.1911049245 | \(6.112439\times10^{-3}\) | 1.1972183746 | \(1.011499\times10^{-6}\) |
| 7 | 0.7 | 1.2776213132 | 1.2709358934 | \(6.685420\times10^{-3}\) | 1.2776228510 | \(1.537844\times10^{-6}\) |
| 8 | 0.8 | 1.3771277643 | 1.3696909578 | \(7.436806\times10^{-3}\) | 1.3771300464 | \(2.282041\times10^{-6}\) |
| 9 | 0.9 | 1.4993025001 | 1.4908915853 | \(8.410915\times10^{-3}\) | 1.4993058361 | \(3.336036\times10^{-6}\) |
| 10 | 1.0 | 1.6487212707 | 1.6390551041 | \(9.666167\times10^{-3}\) | 1.6487261040 | \(4.833266\times10^{-6}\) |

### Valores pedidos en x = 0.5

La solución exacta es

\[
y(0.5)=e^{0.5^2/2}=e^{0.125}=1.1331484531.
\]

La aproximación AB2 es

\[
y_{AB2}(0.5)=1.1274641134.
\]

El error absoluto es

\[
E_{AB2}(0.5)
=
|1.1274641134-1.1331484531|
=
5.6843396437\times10^{-3}.
\]

\[
\boxed{y_{AB2}(0.5)=1.1274641134}
\]

\[
\boxed{E_{AB2}(0.5)=5.6843396437\times10^{-3}}
\]

### Valores pedidos en x = 1

\[
y(1)=e^{1^2/2}=e^{0.5}=1.6487212707.
\]

\[
y_{AB2}(1)=1.6390551041.
\]

\[
E_{AB2}(1)
=
|1.6390551041-1.6487212707|
=
9.6661665842\times10^{-3}.
\]

\[
\boxed{y_{AB2}(1)=1.6390551041}
\]

\[
\boxed{E_{AB2}(1)=9.6661665842\times10^{-3}}
\]

---

## Inciso d) Método de RK4

El método de Runge-Kutta de cuarto orden es:

\[
k_1=f(x_n,y_n),
\]

\[
k_2=f\left(x_n+\frac{h}{2},y_n+\frac{h}{2}k_1\right),
\]

\[
k_3=f\left(x_n+\frac{h}{2},y_n+\frac{h}{2}k_2\right),
\]

\[
k_4=f(x_n+h,y_n+h k_3),
\]

\[
y_{n+1}=y_n+\frac{h}{6}(k_1+2k_2+2k_3+k_4).
\]

Se usa el mismo paso \(h=0.1\) y la misma condición inicial \(y(0)=1\).
Los resultados se muestran en la última parte de la tabla anterior.

### Valores en x = 0.5

\[
y_{RK4}(0.5)=1.1331490929.
\]

\[
E_{RK4}(0.5)
=
|1.1331490929-1.1331484531|
=
6.3978486531\times10^{-7}.
\]

\[
\boxed{y_{RK4}(0.5)=1.1331490929}
\]

\[
\boxed{E_{RK4}(0.5)=6.3978486531\times10^{-7}}
\]

### Valores en x = 1

\[
y_{RK4}(1)=1.6487261040.
\]

\[
E_{RK4}(1)
=
|1.6487261040-1.6487212707|
=
4.8332664457\times10^{-6}.
\]

\[
\boxed{y_{RK4}(1)=1.6487261040}
\]

\[
\boxed{E_{RK4}(1)=4.8332664457\times10^{-6}}
\]

---

## Inciso e) Gráfico de las tres soluciones

El programa `problema2.c++` guarda los datos en
`resultados_problema2.csv`. Luego, el script `grafico_problema2.py` lee ese
archivo y genera `grafico_problema2.png` con:

- la solución exacta \(y=e^{x^2/2}\);
- la solución aproximada por Adams-Bashforth de 2 pasos;
- la solución aproximada por RK4.

El gráfico se genera con:

```bash
python3 grafico_problema2.py
```

---

## Resumen de resultados

| \(x\) | Exacta | AB2 | Error AB2 | RK4 | Error RK4 |
|---:|---:|---:|---:|---:|---:|
| 0.5 | 1.1331484531 | 1.1274641134 | \(5.684340\times10^{-3}\) | 1.1331490929 | \(6.397849\times10^{-7}\) |
| 1.0 | 1.6487212707 | 1.6390551041 | \(9.666167\times10^{-3}\) | 1.6487261040 | \(4.833266\times10^{-6}\) |

Conclusión: el método propuesto (Adams-Bashforth de 2 pasos) aproxima
razonablemente la solución, pero acumula un error del orden de \(10^{-2}\) al
llegar a \(x=1\). RK4, al ser de cuarto orden, reduce el error a valores del
orden de \(10^{-6}\) con el mismo paso \(h=0.1\).
