# Métodos de Localización de Raíces

Programas en C++ para encontrar raíces de `f(x) = 0` (valores de `x` donde
la función corta al eje `x`). Cada archivo resuelve, por defecto, la función
del Parcial 1 – 2024:

```
f(x) = ln(x^2 + 1) - sin(x)      raiz en [1, 2]  ->  x ~ 1.2620326547
```

Para usar otra función, solo hay que cambiar la definición de `f()` (y en
punto fijo, la de `g()`).

---

## Compilar y ejecutar

Cada método es un programa independiente. En la terminal, desde esta carpeta:

```bash
g++ -o biseccion      biseccion.cpp      -lm   && ./biseccion
g++ -o falsaPosicion  falsaPosicion.cpp  -lm   && ./falsaPosicion
g++ -o newtonRapshon  newtonRapshon.cpp  -lm   && ./newtonRapshon
g++ -o puntoFijo      puntoFijo.cpp      -lm   && ./puntoFijo
g++ -o secante        secante.cpp        -lm   && ./secante
```

También existe `TodosMetodos.cpp`, un menú con los cinco métodos juntos
(resuelve `f(x) = x^3 - 9x + 3`).

---

## Resumen de métodos

| Método            | Archivo           | Tipo     | Necesita                    | Convergencia |
|-------------------|-------------------|----------|-----------------------------|--------------|
| Bisección         | `biseccion.cpp`   | Cerrado  | intervalo `[a, b]` con cambio de signo | Lenta, pero segura |
| Regula Falsi      | `falsaPosicion.cpp` | Cerrado| intervalo `[a, b]` con cambio de signo | Media |
| Punto Fijo        | `puntoFijo.cpp`   | Abierto  | punto inicial `x0` y una `g(x)` bien elegida | Si \|g'(x)\|<1 |
| Newton–Raphson    | `newtonRapshon.cpp` | Abierto | punto inicial `x0`          | Rápida (cuadrática) |
| Secante           | `secante.cpp`     | Abierto  | dos puntos `x0` y `x1`      | Rápida |

- **Métodos cerrados**: siempre convergen si `f(a)·f(b) < 0`.
- **Métodos abiertos**: convergen más rápido, pero el punto inicial debe
  estar cerca de la raíz, si no pueden divergir.

---

## Bisección — `biseccion.cpp`

**Qué hace:** parte del intervalo `[a, b]` donde `f(a)·f(b) < 0` y lo parte
siempre por la mitad, quedándose con la mitad que contiene la raíz.

**Cómo calcula:**
```
c = (a + b) / 2
si f(a)·f(c) < 0  ->  la raiz esta en [a, c]  (b = c)
si no             ->  la raiz esta en [c, b]  (a = c)
```

**Qué ingresar:**
- `a` y `b`: extremos del intervalo (debe haber cambio de signo).
- `tolerancia` en **%**: condición de corte (error relativo entre
  iteraciones).

**Notas:** no usa derivada. Cada iteración reduce el intervalo a la mitad,
entonces después de `n` iteraciones el error es `(b−a)/2ⁿ`.

---

## Regula Falsi — `falsaPosicion.cpp`

**Qué hace:** como bisección, pero en vez del punto medio usa la intersección
de la **cuerda** (recta entre `(a, f(a))` y `(b, f(b))`) con el eje `x`.

**Cómo calcula:**
```
c = (a·f(b) - b·f(a)) / (f(b) - f(a))
```
(si `f(b) == f(a)` la cuerda es horizontal y el método falla).

Luego acota el intervalo igual que bisección, según el signo de `f(a)·f(c)`.

**Qué ingresar:** `a`, `b` y `tolerancia` (en %).

**Notas:** más rápido que bisección, pero puede quedar un extremo "clavado"
(no se mueve), como se ve en la tabla cuando `b` queda fijo en `2`.

---

## Punto Fijo — `puntoFijo.cpp`

**Qué hace:** reescribe `f(x) = 0` como `x = g(x)` y genera la sucesión
`x_{k+1} = g(x_k)`. Si converge, el límite es la raíz.

Para el parcial se despejó:
```
ln(x^2+1) - sin(x) = 0  ->  x^2 + 1 = e^{sin(x)}
                        ->  g(x) = sqrt(e^{sin(x)} - 1)
```

**Cómo calcula:**
```
x_{k+1} = g(x_k)
error = |x_{k+1} - x_k|
```

**Qué ingresar:**
- `xviejo`: punto inicial.
- `tolerancia`: condición de corte (error absoluto).

**Notas:**
- La derivada es **numérica**, sirve para cualquier `g`.
- El criterio de convergencia es `|g'(x)| < 1` cerca de la raíz. Es
  condición **suficiente, no necesaria**: el código avisa si `|g'(x0)| ≥ 1`
  pero igual itera (a veces converge).
- **Importante:** la `g` debe elegirse bien. Si la `g` está mal elegida (o
  el despeje es otro), el método puede no converger a la raíz buscada.

---

## Newton–Raphson — `newtonRapshon.cpp`

**Qué hace:** usa la tangente a `f` en `x0` para aproximar el cruce por cero
y repite con la nueva aproximación.

**Cómo calcula:**
```
x_{n+1} = x_n - f(x_n) / f'(x_n)
error = |x_{n+1} - x_n|
```

La derivada `f'` se aproxima numéricamente con `h = 0.001`:
```
f'(x) ~ (f(x + h) - f(x)) / h
```
por eso funciona con **cualquier** función, sin calcular la derivada a mano.

**Qué ingresar:** `x0` (punto inicial) y `tolerancia`.

**Notas:**
- Convergencia cuadrática (muy rápida) si `x0` está cerca de la raíz.
- Falla si `f'(x) ≈ 0` (tangente horizontal): la fórmula divide por
  (casi) cero. El código avisa y sigue.
- Si querés la derivada analítica, usá la de `problema1.c++`:
  `f'(x) = 2x/(x^2+1) - cos(x)`.

---

## Secante — `secante.cpp`

**Qué hace:** igual que Newton, pero reemplaza `f'(x)` por el cociente
incremental entre las dos últimas aproximaciones. No usa derivada.

**Cómo calcula:**
```
x_{n+1} = x_n - f(x_n)·(x_n - x_{n-1}) / (f(x_n) - f(x_{n-1}))
error = |x_{n+1} - x_n|
```

**Qué ingresar:** `x0` y `x1` (dos puntos iniciales) y `tolerancia`.

**Notas:**
- Falla si `f(x1) == f(x0)` (secante horizontal); el código avisa.
- Velocidad casi tan buena como Newton, sin necesidad de derivada.

---

## Consejos generales

1. **Antes de correr un método cerrado**, comprobá `f(a)·f(b) < 0`
   (bisección y regula falsi ya lo verifican solos).
2. **Tolerancia**: en bisección/regula falsi se ingresa en **%** (ej. `0.01`
   = 0.01 %); en los demás es un error absoluto (ej. `0.0000001`).
3. Si el método **no converge**, probá otro punto/intervalo inicial más
   cerca de la raíz, o chequeá que la `g` del punto fijo esté bien despejada.