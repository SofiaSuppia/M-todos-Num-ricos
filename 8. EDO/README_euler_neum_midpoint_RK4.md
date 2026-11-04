# `euler_neum_midpoint_RK4.c++`

## 1. Qué resuelve

Este programa resuelve **una sola ecuación diferencial ordinaria (EDO) escalar** de primer orden:

\[
y' = f(x,y), \qquad y(x_0)=Y_0
\]

en un intervalo \([x_0, x_f]\). Puede usar cuatro métodos numéricos:

1. Euler (Runge-Kutta de orden 1)
2. Heun (Runge-Kutta de orden 2)
3. Punto Medio (Runge-Kutta de orden 2)
4. Runge-Kutta de orden 4 (RK4)

El programa también calcula la solución exacta, el error absoluto, una estimación del error de truncamiento local, el factor de convergencia \(Q\) y guarda los resultados en archivos de texto.

> **Importante:** este código **no resuelve sistemas de EDOs**. Para un sistema como
> \(y_1'=3x+y_2,\; y_2'=x^2-y_1-1\), usar `sistemas_edo.c++`.

---

## 2. Ecuación de ejemplo implementada

Actualmente el código resuelve:

\[
y' = (2x+1)\sqrt{y}, \qquad y(0)=1
\]

La solución exacta es:

\[
y(x)=\frac{(x^2+x+2)^2}{4}
\]

Por eso, en el código:

- `f(x, y)` devuelve \((2x+1)\sqrt{y}\).
- `y_exacta(x)` devuelve \(\frac{(x^2+x+2)^2}{4}\).
- `f_prima(x, y)` calcula \(y''\), usada para el error local de Euler.
- `y3_prima(x, y)` calcula \(y'''=6x+3\), usada para los errores locales de Heun y Punto Medio.

Si se cambia el problema, hay que actualizar **las cuatro funciones**. En particular, `y3_prima` ya no será \(6x+3\) para otra solución exacta.

---

## 3. Archivos relacionados

| Archivo | Función |
|---|---|
| `euler_neum_midpoint_RK4.c++` | Código fuente principal. |
| `graph_points.py` | Gráfica la solución numérica contra la exacta leyendo `results.txt`. |
| `graph_convergence.py` | Gráfica el factor de convergencia leyendo un archivo `convergence_*.txt`. |
| `results.txt` | Se genera al ejecutar el programa. |
| `convergence_euler.txt` | Se genera si se pide convergencia con Euler. |
| `convergence_heun.txt` | Se genera si se pide convergencia con Heun. |
| `convergence_midpoint.txt` | Se genera si se pide convergencia con Punto Medio. |
| `convergence_rk4.txt` | Se genera si se pide convergencia con RK4. |

Los scripts de Python necesitan `numpy` y `matplotlib`.

---

## 4. Compilación

Desde la carpeta `8. EDO`:

```bash
g++ -Wall -Wextra -o euler_neum_midpoint_RK4 euler_neum_midpoint_RK4.c++
```

En Windows con MinGW:

```bash
g++ -Wall -Wextra -o euler_neum_midpoint_RK4.exe euler_neum_midpoint_RK4.c++
```

También se puede compilar dentro de `8. EDO/output`:

```bash
g++ -Wall -Wextra -o output/euler_neum_midpoint_RK4 euler_neum_midpoint_RK4.c++
```

El programa busca los scripts de Python primero en la carpeta actual y, si no los encuentra, en la carpeta padre. Por eso funciona tanto ejecutándolo desde `8. EDO` como desde `8. EDO/output`.

---

## 5. Ejecución

Desde `8. EDO`:

```bash
./euler_neum_midpoint_RK4
```

Desde `8. EDO/output`:

```bash
./euler_neum_midpoint_RK4
```

Los archivos `results.txt` y `convergence_*.txt` se guardan en la carpeta desde la que se ejecuta el programa.

---

## 6. Datos de entrada

El programa pide los datos en este orden:

1. **`X0` y `Xf`**: inicio y fin del intervalo.
2. **`Y0`**: valor inicial \(y(x_0)\).
3. **Tipo de paso**:
   - `1` para ingresar el número de subintervalos `n`.
   - `2` para ingresar el tamaño de paso `h`.
4. **`n` o `h`**, según la opción anterior.
5. **Método**:
   - `1` Euler
   - `2` Heun
   - `3` Punto Medio
   - `4` RK4
6. **Factor de convergencia**:
   - `1` calcularlo.
   - `2` no calcularlo.
7. **Continuar**: `s` para otro cálculo, `n` para salir.

### Ejemplo de entrada: RK4 con \(n=10\), sin convergencia

```text
0 1
1
1
10
4
2
n
```

Interpretación:

| Línea | Valor | Significado |
|---|---|---|
| 1 | `0 1` | intervalo \([0,1]\) |
| 2 | `1` | \(y(0)=1\) |
| 3 | `1` | se ingresará `n` |
| 4 | `10` | \(n=10\), entonces \(h=0.1\) |
| 5 | `4` | método RK4 |
| 6 | `2` | no calcular convergencia |
| 7 | `n` | no hacer otro cálculo |

### Ejemplo de entrada: Euler con \(n=10\), con convergencia

```text
0 1
1
1
10
1
1
n
```

La línea `1` después del método pide calcular el factor de convergencia.

### Validaciones

- Debe cumplirse `Xf > X0`.
- `Y0` debe ser un número real.
- Si se ingresa `n`, debe cumplirse `1 <= n <= 200`.
- Si se ingresa `h`, debe ser `h > 0` y el `n` calculado no puede superar 200.
- En la rama de `h`, se calcula primero `n_real = (Xf-X0)/h` como `double` y luego se redondea, evitando desbordamientos.

---

## 7. Qué hace el programa paso a paso

1. Lee el intervalo, la condición inicial y el tipo de paso.
2. Calcula `h` y `n`.
3. Lee el método y si se desea calcular convergencia.
4. Inicializa:

   ```cpp
   X[0] = X0;
   Y[0] = Y0;
   ```

5. Avanza desde `i = 0` hasta `n-1` aplicando el método elegido.
6. Para cada punto calcula:
   - solución exacta `y_exacta(X[i])`,
   - error absoluto `|y_exacta(X[i]) - Y[i]|`,
   - estimación del error de truncamiento local.
7. Muestra una tabla con columnas:
   - `i`
   - `X[i]`
   - `Y Exacta`
   - `Y` del método
   - `Error Exacto`
   - `Error Trunc. Local`
8. Muestra un resumen final en `X[n]` con:
   - valor numérico,
   - valor exacto,
   - error absoluto,
   - error relativo porcentual.
9. Guarda `results.txt`.
10. Si se pidió convergencia, ejecuta el cálculo con `h`, `h/2` y `h/4`, muestra la tabla de \(Q\), guarda `convergence_*.txt` y llama al script de Python.

---

## 8. Métodos implementados

### Euler

\[
y_{i+1}=y_i+h f(x_i,y_i)
\]

Usa la pendiente al inicio del intervalo. Es de orden 1.

### Heun

\[
y_p=y_i+h f(x_i,y_i)
\]

\[
y_{i+1}=y_i+\frac{h}{2}\left[f(x_i,y_i)+f(x_{i+1},y_p)\right]
\]

Primero predice con Euler y luego corrige promediando las pendientes. Es de orden 2.

### Punto Medio

\[
k_1=f(x_i,y_i)
\]

\[
k_2=f\left(x_i+\frac{h}{2},\,y_i+\frac{h}{2}k_1\right)
\]

\[
y_{i+1}=y_i+h k_2
\]

Evalúa la pendiente en el punto medio del intervalo. Es de orden 2.

### RK4

\[
k_1=f(x_i,y_i)
\]

\[
k_2=f\left(x_i+\frac{h}{2},\,y_i+\frac{h}{2}k_1\right)
\]

\[
k_3=f\left(x_i+\frac{h}{2},\,y_i+\frac{h}{2}k_2\right)
\]

\[
k_4=f(x_i+h,\,y_i+h k_3)
\]

\[
y_{i+1}=y_i+\frac{h}{6}(k_1+2k_2+2k_3+k_4)
\]

Es de orden 4.

---

## 9. Errores y factor de convergencia

### Error exacto

\[
e_i = |y_{\text{exacta}}(x_i)-y_i|
\]

### Error de truncamiento local mostrado

- Euler: \(\frac{h^2}{2}y''(x_i)\). Se muestra con signo.
- Heun: \(\frac{h^3}{12}y'''(x_i)\). Se muestra con signo.
- Punto Medio: \(\frac{h^3}{24}y'''(x_i)\). Se muestra con signo.
- RK4: se estima comparando un paso `h` con dos pasos `h/2`:

\[
\text{LTE} \approx \frac{|y_{h/2}-y_h|}{2^4-1}=\frac{|y_{h/2}-y_h|}{15}
\]

### Factor de convergencia \(Q\)

Para cada punto se resuelve el problema con pasos `h`, `h/2` y `h/4`. Luego:

\[
Q=\frac{\log\left(\frac{|y_h-y_{h/2}|}{|y_{h/2}-y_{h/4}|}\right)}{\log(2)}
\]

Interpretación:

- \(Q\approx 1\): orden 1 (Euler).
- \(Q\approx 2\): orden 2 (Heun, Punto Medio).
- \(Q\approx 4\): orden 4 (RK4).

Si el numerador o el denominador son prácticamente cero, \(Q\) no es confiable y se muestra como `N/A`.

---

## 10. Archivos de salida

### `results.txt`

Tiene encabezado y cuatro columnas separadas por tabulaciones:

```text
x    y_exacta    y_numerica    error_exacto
```

`graph_points.py` lee este archivo con `skiprows=1` y usa:

- columna 0: `x`,
- columna 1: solución exacta,
- columna 2: solución numérica.

### `convergence_*.txt`

Tiene dos columnas:

```text
x    Q
```

`graph_convergence.py` recibe el nombre del archivo como argumento. Por ejemplo:

```bash
python3 graph_convergence.py convergence_euler.txt
```

El programa C++ ya pasa automáticamente el archivo correcto según el método.

---

## 11. Cómo cambiar el problema

Editar estas funciones:

```cpp
double f(double x, double y);
double y_exacta(double x);
double f_prima(double x, double y);
double y3_prima(double x, double y);
```

- `f`: lado derecho de la EDO.
- `y_exacta`: solución analítica. Si no se conoce, los errores calculados no tendrán sentido; habría que modificar la impresión/guardado de errores.
- `f_prima`: \(y''\), necesaria para el error local de Euler.
- `y3_prima`: \(y'''\), necesaria para los errores locales de Heun y Punto Medio.

También se debe verificar que la solución no salga del dominio de `f`. Por ejemplo, el ejemplo actual usa `sqrt(y)`, así que se necesita \(y\ge 0\).

---

## 12. Por qué el ejemplo de sistema no funciona aquí

El código escalar está diseñado para **una sola variable dependiente**:

- tiene un único `Y0`;
- tiene un arreglo `Y[i]` de una dimensión;
- `f` recibe y devuelve un solo `double`;
- `y_exacta` devuelve un solo `double`.

El sistema

\[
y_1'=3x+y_2,\qquad y_2'=x^2-y_1-1
\]

necesita dos variables dependientes y dos condiciones iniciales. No se puede representar con un solo `Y0` ni con una función `f(double x, double y)`. Para ese caso se debe usar `sistemas_edo.c++`, que trabaja con arreglos `y[N]`, `Y[i][j]` y una función `f(x, y[], dy[])`.

---

## 13. Ejemplo de resultados

Para el problema actual, \([0,1]\), \(n=10\), \(y(0)=1\):

| Método | \(y(1)\) aproximada | Error absoluto |
|---|---:|---:|
| Euler | 3.678572 | \(3.214285\times 10^{-1}\) |
| Heun | 3.995105 | \(4.895394\times 10^{-3}\) |
| Punto Medio | 3.991063 | \(8.936679\times 10^{-3}\) |
| RK4 | 3.999995 | \(4.552175\times 10^{-6}\) |
| Exacta | 4.000000 | 0 |

El orden esperado se refleja en la magnitud del error: RK4 es mucho más preciso que Euler para el mismo `n`.

---

## 14. Limitaciones

- `n` no puede superar 200 porque los arreglos tienen tamaño `MAX_SIZE + 1`.
- El cálculo de convergencia usa hasta `4*n` subintervalos; por eso también se valida `n <= 200`.
- Se requiere una solución exacta implementada para calcular errores.
- Las gráficas requieren `numpy` y `matplotlib`.
- En Linux, si no hay entorno gráfico, `plt.show()` puede no abrir una ventana; el programa igualmente termina.
