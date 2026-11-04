# `sistemas_edo.c++`

## 1. Qué resuelve

Este programa resuelve un **sistema de \(N\) ecuaciones diferenciales ordinarias de primer orden**:

\[
\begin{aligned}
y_1' &= f_1(x,y_1,\dots,y_N)\\
y_2' &= f_2(x,y_1,\dots,y_N)\\
&\vdots\\
y_N' &= f_N(x,y_1,\dots,y_N)
\end{aligned}
\]

con condiciones iniciales:

\[
y_1(x_0)=y_{10},\quad y_2(x_0)=y_{20},\quad\dots,\quad y_N(x_0)=y_{N0}
\]

Métodos disponibles:

1. Euler
2. Heun
3. Punto Medio
4. Runge-Kutta de orden 4 (RK4)

El programa calcula la solución numérica, la solución exacta, los errores absolutos y, opcionalmente, el factor de convergencia \(Q\) por componente.

---

## 2. Diferencia entre \(N\) y \(n\)

Es muy importante no confundirlos:

- **\(N\)** = número de ecuaciones del sistema. También es el orden de la EDO original cuando el sistema proviene de reducir una EDO de orden \(N\).
- **\(n\)** = número de subintervalos que usa el método numérico en \([x_0,x_f]\).

Por ejemplo, para resolver \(y'''=6\) se usa un sistema de \(N=3\) ecuaciones. Si se elige \(n=10\), el método da 10 pasos con \(h=(x_f-x_0)/10\).

---

## 3. Ejemplos incluidos

El archivo tiene bloques para \(N=2\), \(N=3\) y \(N=4\). Se selecciona cambiando la línea:

```cpp
#define N 2
```

por `#define N 3` o `#define N 4`, y luego se vuelve a compilar.

### Ejemplo \(N=2\)

Sistema:

\[
\begin{aligned}
y_1' &= 3x+y_2\\
y_2' &= x^2-y_1-1
\end{aligned}
\]

Condiciones iniciales:

\[
y_1(0)=0,\qquad y_2(0)=1
\]

Solución exacta:

\[
y_1=x^2+\sin(x),\qquad y_2=\cos(x)-x
\]

Datos de entrada para \(x\in[0,1]\), \(n=10\), método RK4, sin convergencia:

```text
0 1
1
10
0
1
4
2
n
```

Interpretación:

| Línea | Valor | Significado |
|---|---|---|
| 1 | `0 1` | intervalo \([0,1]\) |
| 2 | `1` | se ingresará `n` |
| 3 | `10` | \(n=10\) |
| 4 | `0` | \(y_1(0)=0\) |
| 5 | `1` | \(y_2(0)=1\) |
| 6 | `4` | RK4 |
| 7 | `2` | no calcular convergencia |
| 8 | `n` | no hacer otro cálculo |

### Ejemplo \(N=3\)

El sistema representa la EDO de tercer orden:

\[
y'''=6
\]

Definiendo:

\[
y_1=y,\qquad y_2=y',\qquad y_3=y''
\]

se obtiene:

\[
\begin{aligned}
y_1' &= y_2\\
y_2' &= y_3\\
y_3' &= 6
\end{aligned}
\]

Condiciones iniciales:

\[
y_1(0)=0,\qquad y_2(0)=0,\qquad y_3(0)=0
\]

Solución exacta:

\[
y_1=x^3,\qquad y_2=3x^2,\qquad y_3=6x
\]

En \(x=1\):

\[
y_1(1)=1,\qquad y_2(1)=3,\qquad y_3(1)=6
\]

Datos de entrada para \(x\in[0,1]\), \(n=10\), método RK4, sin convergencia:

```text
0 1
1
10
0
0
0
4
2
n
```

### Ejemplo \(N=4\)

El sistema representa la EDO de cuarto orden:

\[
y''''=24
\]

Definiendo:

\[
y_1=y,\qquad y_2=y',\qquad y_3=y'',\qquad y_4=y'''
\]

se obtiene:

\[
\begin{aligned}
y_1' &= y_2\\
y_2' &= y_3\\
y_3' &= y_4\\
y_4' &= 24
\end{aligned}
\]

Condiciones iniciales:

\[
y_1(0)=0,\qquad y_2(0)=0,\qquad y_3(0)=0,\qquad y_4(0)=0
\]

Solución exacta:

\[
y_1=x^4,\qquad y_2=4x^3,\qquad y_3=12x^2,\qquad y_4=24x
\]

En \(x=1\):

\[
y_1(1)=1,\qquad y_2(1)=4,\qquad y_3(1)=12,\qquad y_4(1)=24
\]

Datos de entrada para \(x\in[0,1]\), \(n=10\), método RK4, sin convergencia:

```text
0 1
1
10
0
0
0
0
4
2
n
```

> **Nota:** las soluciones exactas corresponden a las condiciones iniciales del ejemplo. Si se ingresan otras condiciones iniciales, la función `exacta` ya no representará la solución del problema y los errores mostrados no serán válidos.

---

## 4. Compilación

Para \(N=2\), el archivo ya tiene `#define N 2`. Compilar desde `8. EDO`:

```bash
g++ -Wall -Wextra -o sistemas_edo sistemas_edo.c++
```

Para \(N=3\) o \(N=4\):

1. Editar `#define N 2` y poner `#define N 3` o `#define N 4`.
2. Compilar de nuevo:

```bash
g++ -Wall -Wextra -o sistemas_edo sistemas_edo.c++
```

En Windows:

```bash
g++ -Wall -Wextra -o sistemas_edo.exe sistemas_edo.c++
```

Si se pone un \(N\) distinto de 2, 3 o 4, el compilador muestra un error indicando que hay que agregar un bloque `#elif N == ...` con las funciones `f` y `exacta`.

---

## 5. Ejecución

```bash
./sistemas_edo
```

El programa guarda los resultados en `resultados_sistema.txt` en la carpeta desde la que se ejecuta.

---

## 6. Datos de entrada

El orden de entrada es:

1. **`x0` y `xf`**: intervalo de integración.
2. **Tipo de paso**:
   - `1` para ingresar `n`.
   - `2` para ingresar `h`.
3. **`n` o `h`**, según la opción anterior.
4. **Condiciones iniciales**: \(y_1(x_0), y_2(x_0), \dots, y_N(x_0)\), una por línea o separadas por espacios.
5. **Método**:
   - `1` Euler
   - `2` Heun
   - `3` Punto Medio
   - `4` RK4
6. **Factor de convergencia**:
   - `1` calcularlo.
   - `2` no calcularlo.
7. **Continuar**: `s` para otro cálculo, `n` para salir.

### Validaciones

- Debe cumplirse `xf > x0`.
- Si se ingresa `n`, debe cumplirse `1 <= n <= 200`.
- Si se ingresa `h`, debe ser `h > 0` y el `n` calculado no puede superar 200.
- El método debe estar entre 1 y 4.

---

## 7. Qué hace el programa paso a paso

1. Define `f` y `exacta` según el valor de `N`.
2. Lee el intervalo, el tipo de paso, `n` o `h`, las condiciones iniciales y el método.
3. Opcionalmente calcula el factor de convergencia.
4. Ejecuta `resolver`, que aplica `paso_metodo` desde `i=0` hasta `n-1`:

   ```cpp
   X[0] = x0;
   Y[0][j] = y0[j];
   for (i = 0; i < n; i++) {
       paso_metodo(metodo, X[i], Y[i], h, Y[i+1]);
       X[i+1] = X[i] + h;
   }
   ```

5. Muestra una tabla con:
   - `i`, `x`,
   - valores numéricos \(y_1,\dots,y_N\),
   - valores exactos \(y_1^{ex},\dots,y_N^{ex}\).
6. Muestra una tabla de errores absolutos por componente.
7. Muestra un resumen final en \(x_f\) con:
   - valor numérico,
   - valor exacto,
   - error absoluto,
   - error relativo porcentual.
8. Guarda `resultados_sistema.txt`.

---

## 8. Métodos implementados

En forma vectorial, con \(\mathbf{y}=(y_1,\dots,y_N)^T\) y \(\mathbf{f}(x,\mathbf{y})=(f_1,\dots,f_N)^T\):

### Euler

\[
\mathbf{y}_{i+1}=\mathbf{y}_i+h\mathbf{f}(x_i,\mathbf{y}_i)
\]

### Heun

\[
\mathbf{y}_p=\mathbf{y}_i+h\mathbf{f}(x_i,\mathbf{y}_i)
\]

\[
\mathbf{y}_{i+1}=\mathbf{y}_i+\frac{h}{2}\left[\mathbf{f}(x_i,\mathbf{y}_i)+\mathbf{f}(x_{i+1},\mathbf{y}_p)\right]
\]

### Punto Medio

\[
\mathbf{k}_1=\mathbf{f}(x_i,\mathbf{y}_i)
\]

\[
\mathbf{k}_2=\mathbf{f}\left(x_i+\frac{h}{2},\,\mathbf{y}_i+\frac{h}{2}\mathbf{k}_1\right)
\]

\[
\mathbf{y}_{i+1}=\mathbf{y}_i+h\mathbf{k}_2
\]

### RK4

\[
\mathbf{k}_1=\mathbf{f}(x_i,\mathbf{y}_i)
\]

\[
\mathbf{k}_2=\mathbf{f}\left(x_i+\frac{h}{2},\,\mathbf{y}_i+\frac{h}{2}\mathbf{k}_1\right)
\]

\[
\mathbf{k}_3=\mathbf{f}\left(x_i+\frac{h}{2},\,\mathbf{y}_i+\frac{h}{2}\mathbf{k}_2\right)
\]

\[
\mathbf{k}_4=\mathbf{f}(x_i+h,\,\mathbf{y}_i+h\mathbf{k}_3)
\]

\[
\mathbf{y}_{i+1}=\mathbf{y}_i+\frac{h}{6}(\mathbf{k}_1+2\mathbf{k}_2+2\mathbf{k}_3+\mathbf{k}_4)
\]

Cada componente se actualiza con la misma fórmula, usando el vector de pendientes correspondiente.

---

## 9. Errores y factor de convergencia

### Error absoluto por componente

\[
e_{i,j}=|y_j^{ex}(x_i)-y_{i,j}|
\]

### Error relativo final

\[
\text{error relativo}_j=\frac{|y_j^{ex}(x_f)-y_{n,j}|}{|y_j^{ex}(x_f)|}\cdot 100\%
\]

Si el valor exacto es prácticamente cero, el error relativo se deja en 0 para evitar división por cero.

### Factor de convergencia \(Q\)

Para cada componente \(j\), se resuelve el sistema con pasos `h`, `h/2` y `h/4` y se calcula:

\[
Q_j=\frac{\log\left(\frac{|y_{h,j}-y_{h/2,j}|}{|y_{h/2,j}-y_{h/4,j}|}\right)}{\log(2)}
\]

Interpretación:

- \(Q_j\approx 1\): orden 1 (Euler).
- \(Q_j\approx 2\): orden 2 (Heun, Punto Medio).
- \(Q_j\approx 4\): orden 4 (RK4).

Si el numerador o el denominador son casi cero, \(Q_j\) no es confiable y se muestra como `N/A`.

---

## 10. Archivos de salida

### `resultados_sistema.txt`

Tiene encabezado y columnas:

```text
x    y1 ... yN    y1_exacta ... yN_exacta    error1 ... errorN
```

Ejemplo para \(N=2\):

```text
x    y1    y2    y1_exacta    y2_exacta    error1    error2
```

### `convergence_*.txt`

Si se pide convergencia, se genera un archivo según el método:

- `convergence_euler.txt`
- `convergence_heun.txt`
- `convergence_punto_medio.txt`
- `convergence_rk4.txt`

El formato es:

```text
i    x    Q1 ... QN
```

A diferencia del código escalar, este programa no llama a Python; los archivos quedan disponibles para analizarlos o graficarlos por separado.

---

## 11. Cómo cambiar el sistema

Dentro del bloque activo (`#if N == ...`) editar:

```cpp
void f(double x, const double y[], double dy[]);
void exacta(double x, double ye[]);
```

- `f` debe llenar `dy[0]` a `dy[N-1]` con las derivadas.
- `exacta` debe llenar `ye[0]` a `ye[N-1]` con la solución exacta.

Si no se conoce la solución exacta:

- se puede dejar una función `exacta` aproximada o nula,
- pero los errores y el resumen final no tendrán significado,
- o se puede modificar el código para ocultar/eliminar el cálculo de errores.

Si se quiere otro valor de \(N\):

1. Agregar un bloque `#elif N == nuevo_valor`.
2. Definir `f` y `exacta` para ese tamaño.
3. Cambiar `#define N` y recompilar.

---

## 12. Ejemplos de resultados

### \(N=2\), RK4, \(n=10\)

Entrada:

```text
0 1
1
10
0
1
4
2
n
```

Resultados finales:

| Componente | Numérica | Exacta | Error absoluto |
|---|---:|---:|---:|
| \(y_1\) | 1.8414696085 | 1.8414709848 | \(1.376259\times 10^{-6}\) |
| \(y_2\) | -0.4596988324 | -0.4596976941 | \(1.138240\times 10^{-6}\) |

### \(N=3\), RK4, \(n=10\)

Entrada:

```text
0 1
1
10
0
0
0
4
2
n
```

Resultados finales:

| Componente | Numérica | Exacta |
|---|---:|---:|
| \(y_1\) | 1.0000000000 | 1.0000000000 |
| \(y_2\) | 3.0000000000 | 3.0000000000 |
| \(y_3\) | 6.0000000000 | 6.0000000000 |

El error es del orden de \(10^{-15}\), es decir, prácticamente error de redondeo.

### \(N=4\), RK4, \(n=10\)

Entrada:

```text
0 1
1
10
0
0
0
0
4
2
n
```

Resultados finales:

| Componente | Numérica | Exacta |
|---|---:|---:|
| \(y_1\) | 1.0000000000 | 1.0000000000 |
| \(y_2\) | 4.0000000000 | 4.0000000000 |
| \(y_3\) | 12.0000000000 | 12.0000000000 |
| \(y_4\) | 24.0000000000 | 24.0000000000 |

Para este ejemplo polinómico, RK4 resulta exacto salvo por redondeo. Por eso, al pedir convergencia con RK4, \(Q\) puede aparecer como `N/A` o 0. Para observar \(Q\approx 1\), usar Euler; para \(Q\approx 2\), usar Heun o Punto Medio.

---

## 13. Por qué el ejemplo de sistema no funciona en el código escalar

El archivo `euler_neum_midpoint_RK4.c++` está hecho para una sola EDO escalar:

- tiene un solo `Y0`;
- usa un arreglo `Y[i]` de una dimensión;
- su función es `f(double x, double y)`;
- su solución exacta es `y_exacta(double x)`, que devuelve un solo número.

El sistema \(N=2\) necesita dos variables dependientes, dos condiciones iniciales y una función que devuelva dos derivadas. Eso no cabe en la estructura escalar. Por eso el ejemplo \(y_1'=3x+y_2,\;y_2'=x^2-y_1-1\) debe resolverse con `sistemas_edo.c++`, que usa `f(x, y[], dy[])`, `Y[i][j]` y `y0[N]`.

---

## 14. Limitaciones

- \(N\) se fija en tiempo de compilación mediante `#define N`.
- `n` no puede superar 200 por el tamaño de los arreglos (`MAX_PASOS`).
- El cálculo de convergencia usa hasta `4*n` subintervalos; por eso también se valida `n <= 200`.
- Se requiere una solución exacta coherente con las condiciones iniciales para que los errores sean válidos.
- Los ejemplos \(N=3\) y \(N=4\) usan condiciones iniciales nulas; si se cambian, también hay que cambiar `exacta`.
- Este programa no genera gráficas automáticamente.
