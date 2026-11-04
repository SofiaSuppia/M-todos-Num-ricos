# Integrar \( \int_0^3 (x^2+1)\,dx \) con Trapecio y Simpson

## Resultado exacto

\[
\int_0^3 (x^2+1)\,dx
=
\left[\frac{x^3}{3}+x\right]_0^3
=
\frac{27}{3}+3
=
12
\]

Por lo tanto, el valor exacto que se debe usar en las opciones simples es:

```text
12
```

## 1. ¿Hay que cambiar la función?

No. Los dos archivos ya tienen definida la función pedida:

```c
double f(double x) {
    return (x * x) + 1.0;
}
```

- En `NewtonCoartsTrapecio.c++` está en la función `f`.
- En `NewtonCoartsSimpson.c++` está en la función `f`.

Solamente habría que cambiar esa función si se quisiera integrar otra expresión.

## 2. Compilar

Desde la carpeta `6. Metodos de Integracion Numerica`:

```bash
g++ -Wall -Wextra -o trapecio NewtonCoartsTrapecio.c++
g++ -Wall -Wextra -o simpson NewtonCoartsSimpson.c++
```

Luego ejecutar:

```bash
./trapecio
./simpson
```

---

# Trapecio

## 3. Trapecio compuesto con función

Elegir:

```text
1
1
0 3
6
n
```

Interpretación de cada dato:

| Dato | Significado |
|---|---|
| `1` | Trapecio compuesto |
| `1` | Tengo una función |
| `0 3` | Límites de integración: \(a=0\), \(b=3\) |
| `6` | Número de subintervalos \(n\) |
| `n` | No realizar otro cálculo |

Resultado esperado con \(n=6\):

```text
La integral es: 12.125000
```

Con otros valores de \(n\):

| \(n\) | Resultado aproximado |
|---|---|
| 3 | 12.500000 |
| 6 | 12.125000 |
| 10 | 12.030000 |

El trapecio compuesto no exige que \(n\) sea par.

## 4. Trapecio simple con errores

Elegir:

```text
2
0 3
1.5
12
n
```

Interpretación:

| Dato | Significado |
|---|---|
| `2` | Trapecio simple |
| `0 3` | Límites de integración: \(a=0\), \(b=3\) |
| `1.5` | Valor \(c\) dentro de \([0,3]\) para el error aproximado |
| `12` | Valor exacto de la integral |
| `n` | No realizar otro cálculo |

Resultado esperado:

```text
La integral aproximada es: 16.500000
El error aproximado es: 4.500000
El error exacto es: 4.500000
El error porcentual es: 37.500000%
```

---

# Simpson

## 5. Simpson compuesto con función

Elegir:

```text
1
1
0 3
2
n
```

Interpretación:

| Dato | Significado |
|---|---|
| `1` | Simpson compuesto |
| `1` | Tengo una función |
| `0 3` | Límites de integración: \(a=0\), \(b=3\) |
| `2` | Número de subintervalos \(n\), debe ser par |
| `n` | No realizar otro cálculo |

Resultado esperado con \(n=2\):

```text
La integral es: 12.000000
```

Como \(x^2+1\) es un polinomio de grado 2, Simpson 1/3 lo integra exactamente. Por eso, con cualquier \(n\) par se obtiene 12:

| \(n\) | Resultado |
|---|---|
| 2 | 12.000000 |
| 4 | 12.000000 |
| 6 | 12.000000 |

Si se ingresa un \(n\) impar, el programa muestra:

```text
ERROR: El numero de subintervalos debe ser PAR.
```

## 6. Simpson 1/3 simple con errores

Elegir:

```text
2
0 3
1.5
12
n
```

Interpretación:

| Dato | Significado |
|---|---|
| `2` | Simpson 1/3 simple |
| `0 3` | Límites de integración: \(a=0\), \(b=3\) |
| `1.5` | Valor \(c\) dentro de \([0,3]\) para el error aproximado |
| `12` | Valor exacto de la integral |
| `n` | No realizar otro cálculo |

Resultado esperado:

```text
La integral aproximada es: 12.000000
El error aproximado es: 0.000000
El error exacto es: 0.000000
El error porcentual es: 0.000000%
```

---

# Uso con tabla de datos

Los dos programas leen la tabla desde:

```c
"../data.txt"
```

Eso significa que:

- Si se ejecuta desde `6. Metodos de Integracion Numerica`, busca `../data.txt` en la carpeta anterior.
- Si se ejecuta desde `6. Metodos de Integracion Numerica/output`, busca `../data.txt`, que es el archivo `data.txt` de la carpeta 6.

El `data.txt` actual de la carpeta 6 contiene:

```text
2
1.0 2.0
3.0 4.0
```

Esos puntos **no** corresponden a \(f(x)=x^2+1\), así que no sirven para comprobar esta integral.

Para probar la tabla con \(x^2+1\), se puede usar este archivo:

```text
7
0.0 1.0
0.5 1.25
1.0 2.0
1.5 3.25
2.0 5.0
2.5 7.25
3.0 10.0
```

Ese archivo tiene 7 puntos y 6 subintervalos de paso \(h=0.5\).

## 7. Trapecio compuesto con tabla

Elegir:

```text
1
2
n
```

Los datos se leen automáticamente desde `../data.txt`.

Resultado esperado con la tabla anterior:

```text
La integral es: 12.125000
```

## 8. Simpson compuesto con tabla

Elegir:

```text
1
2
n
```

Resultado esperado con la tabla anterior:

```text
La integral es: 12.000000
```

Simpson necesita una cantidad par de subintervalos. Si la tabla tiene una cantidad impar de intervalos y los puntos están equiespaciados, el programa agrega un punto usando el spline cúbico para poder aplicar Simpson.

---

# Resumen rápido de entradas

## Trapecio compuesto

```text
1
1
0 3
6
n
```

## Trapecio simple

```text
2
0 3
1.5
12
n
```

## Simpson compuesto

```text
1
1
0 3
2
n
```

## Simpson simple

```text
2
0 3
1.5
12
n
```

## Tabla para \(x^2+1\)

```text
7
0.0 1.0
0.5 1.25
1.0 2.0
1.5 3.25
2.0 5.0
2.5 7.25
3.0 10.0
```
