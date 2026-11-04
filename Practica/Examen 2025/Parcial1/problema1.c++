// ============================================================================
//  ENUNCIADO  (Problema n.° 1 — Examen 2025)
//
//  Dada la funcion
//
//       f(x) = (x + 1)/(x + 4) − 0.25·x
//
//  a)  Determinar la raiz positiva de f analiticamente.
//
//  b)  Obtener la raiz del inciso a) empleando el metodo de biseccion,
//      tomando como intervalo inicial [1.75, 3.2], con un error porcentual
//      de 0.01%. Utilice como criterio de error el siguiente:
//
//           ε = ( |x_nuevo − x_viejo| ) / ( 1/2·|x_nuevo + x_viejo| ) · 100   (#)
//
//      donde x_viejo en la primera iteracion es a = 3.2. Escriba el valor
//      obtenido de la raiz (con 10 cifras decimales) asi como el error
//      estimado, el error porcentual exacto y el numero de iteraciones
//      empleadas.
//
//  c1) Realice a mano, paso a paso y en papel, la primera iteracion del
//      metodo de Newton-Raphson para localizar la raiz de f(x), usando como
//      valor inicial x0 = 1.0. Escriba el valor obtenido de la raiz asi como
//      el error exacto y el error estimado (segun el mismo criterio (#)).
//
//  c2) Usando el codigo apropiado, escriba el valor de la raiz obtenida
//      (con 10 cifras decimales) y el error estimado (#), para la cuarta
//      iteracion.
// ============================================================================
#include <cstdio>
#include <cmath>
#include <iostream>

using namespace std;

// ============================================================
// f(x) = (x + 1)/(x + 4) - 0.25*x
// ============================================================
double f(double x) {
    return (x + 1) / (x + 4) - 0.25 * x;
}

// ============================================================
// Derivada NUMERICA (igual que newtonRapshon.cpp del tema
// "1. LocalizacionRaices"): sirve para CUALQUIER funcion f,
// no hay que calcular la derivada a mano.
//   f'(x) ≈ (f(x + h) - f(x)) / h     con h = 0.001
// ============================================================
double df(double x) {
    double h = 0.001;
    return (f(x + h) - f(x)) / h;
}

// ============================================================
// Criterio de error (#) del enunciado:
//   ε = |x_nuevo - x_viejo| / ( 1/2·|x_nuevo + x_viejo| ) · 100
// ============================================================
double errorPorcentual(double xNuevo, double xViejo) {
    return fabs(xNuevo - xViejo) / (0.5 * fabs(xNuevo + xViejo)) * 100.0;
}

int main() {

    // Raiz exacta del inciso a): x = 2 (se demuestra abajo)
    const double raizExacta = 2.0;

    printf("===========================================================\n");
    printf(" PROBLEMA 1 - RAICES DE ECUACIONES NO LINEALES\n");
    printf(" f(x) = (x + 1)/(x + 4) - 0.25·x\n");
    printf("===========================================================\n\n");

    // ============================================================
    // PARTE a) RAIZ POSITIVA ANALITICA
    // ============================================================
    printf("=========== PARTE a) RAIZ ANALITICA ===========\n\n");
    printf("f(x) = 0  ->  (x + 1)/(x + 4) - 0.25·x = 0\n");
    printf("Multiplico por (x + 4) (x != -4):\n");
    printf("  x + 1 - 0.25·x·(x + 4) = 0\n");
    printf("  x + 1 - 0.25·x^2 - x    = 0\n");
    printf("  1 - 0.25·x^2            = 0\n");
    printf("  x^2 = 4  ->  x = +-2\n");
    printf("Raiz POSITIVA:  x = 2.0000000000\n\n");

    // ============================================================
    // PARTE b) BISECCION
    // Intervalo inicial [1.75, 3.2], error porcentual 0.01%.
    // Criterio (#): x_viejo de la primera iteracion es 3.2.
    // ============================================================
    printf("=========== PARTE b) BISECCION ===========\n\n");

    double a = 1.75, b = 3.2;
    printf("f(a) = f(1.75) = %.6lf\n", f(a));
    printf("f(b) = f(3.20) = %.6lf\n", f(b));
    printf("f(a)·f(b) < 0  ->  hay raiz en [1.75, 3.2]\n\n");

    double xViejo = 3.2;   // x_viejo en la primera iteracion (segun el enunciado)
    double xNuevo;
    double error_abs;      // error absoluto estimado (b - a)/2 (como biseccion.cpp)
    double epsEst;         // error estimado, criterio (#), en %
    double epsExacto;      // error porcentual exacto: |x - 2| / 2 * 100
    int iter = 0;

    // Tabla igual al estilo de biseccion.cpp, agregando la columna del
    // error (#) que exige el enunciado y el error porcentual exacto.
    printf("\n iter |      a      |      b      |    f(a)     |      c      |    f(c)     |  error abs  |   error (#)%%  |  error exacto%%\n");
    printf("------+-------------+-------------+-------------+-------------+-------------+-------------+---------------+---------------\n");

    do {
        iter++;

        // Punto medio del intervalo (aproximacion de la raiz)
        xNuevo = (a + b) / 2.0;

        // Error absoluto estimado: la raiz esta a menos de la mitad de la
        // longitud del intervalo respecto de c (igual que biseccion.cpp)
        error_abs = (b - a) / 2.0;

        // Criterio de error (#) del enunciado, contra el valor anterior
        epsEst = errorPorcentual(xNuevo, xViejo);
        // Error porcentual exacto (se conoce la raiz exacta: 2)
        epsExacto = fabs(xNuevo - raizExacta) / fabs(raizExacta) * 100.0;

        printf("%5d | %11.8f | %11.8f | %11.8f | %11.8f | %11.8f | %11.8f | %13.8f | %13.8f\n",
               iter, a, b, f(a), xNuevo, f(xNuevo), error_abs, epsEst, epsExacto);

        // Si f(c) = 0, la raiz es exacta: terminar
        if (f(xNuevo) == 0) {
            epsEst = 0;
            break;
        }

        // Determinar en que subintervalo esta la raiz
        if (f(a) * f(xNuevo) > 0) {
            a = xNuevo;      // la raiz esta en [xNuevo, b]
        } else {
            b = xNuevo;      // la raiz esta en [a, xNuevo]
        }

        xViejo = xNuevo;     // guardar aprox. anterior para el error (#)
    } while (epsEst >= 0.01 && iter < 200);

    printf("\nRESULTADOS DE BISECCION:\n");
    printf("  Raiz obtenida (10 cifras decimales):  %.10lf\n", xNuevo);
    printf("  Error estimado (criterio #):          %.6lf %%\n", epsEst);
    printf("  Error porcentual exacto:              %.6lf %%\n", epsExacto);
    printf("  Numero de iteraciones empleadas:      %d\n", iter);
    printf("  (Raiz exacta: x = 2)\n\n");

    // ============================================================
    // PARTE c1) NEWTON-RAPHSON, PRIMERA ITERACION (x0 = 1.0)
    // La cuenta se hace paso a paso (como se haria en papel).
    // ============================================================
    printf("=========== PARTE c1) NEWTON-RAPHSON - 1ra ITERACION ===========\n");
    printf("(en papel se usa la derivada ANALITICA: f'(x) = 3/(x+4)^2 - 0.25)\n\n");

    double x0 = 1.0;
    // En papel, con x0 = 1:
    //   f(1)   = 2/5 - 1/4 = 3/20        = 0.15
    //   f'(1)  = 3/25 - 1/4 = -13/100    = -0.13
    //   x1 = x0 - f(x0)/f'(x0) = 1 - (3/20)/(-13/100) = 28/13
    double fx0Papel  = 3.0 / 20.0;      // f(1) exacto (analitico)
    double dfx0Papel = -13.0 / 100.0;   // f'(1) exacto (analitico)
    double x1        = x0 - fx0Papel / dfx0Papel;

    printf("  x0 = %.10lf\n", x0);
    printf("  f(x0)  = f(1)   = %.10lf\n", fx0Papel);
    printf("  f'(x0) = f'(1)  = %.10lf\n", dfx0Papel);
    printf("  x1 = x0 - f(x0)/f'(x0)\n");
    printf("     = 1 - (%.10lf)/(%.10lf)\n", fx0Papel, dfx0Papel);
    printf("     = %.10lf\n", x1);
    printf("\n");

    double errExactoC1   = fabs(x1 - raizExacta);                 // error absoluto
    double errPorcExacto = fabs(x1 - raizExacta) / raizExacta * 100.0;
    double errEstC1      = errorPorcentual(x1, x0);               // criterio (#)

    printf("  Valor obtenido (1ra iteracion):  %.10lf\n", x1);
    printf("  Error exacto (absoluto):         %.10lf\n", errExactoC1);
    printf("  Error porcentual exacto:         %.6lf %%\n", errPorcExacto);
    printf("  Error estimado (criterio #):     %.6lf %%\n\n", errEstC1);

    // ============================================================
    // PARTE c2) NEWTON-RAPHSON: CUARTA ITERACION
    // ============================================================
    printf("=========== PARTE c2) NEWTON-RAPHSON - 4ta ITERACION ===========\n");
    printf("(con el codigo se usa la derivada NUMERICA, como newtonRapshon.cpp)\n\n");

    double xAnt = x0;     // para el error estimado (#) entre iteraciones
    double xAct = x0;

    printf(" iter |     x_nuevo     |   f(x_nuevo)     f'(x_nuevo)   "
           "|  error est.(#)\n");
    printf("--------------------------------------------------------------"
           "------------\n");
    printf("    0  |   %.10lf   |   %10.6lf   %10.6lf   |       ----\n",
           x0, f(x0), df(x0));

    double eps4 = 0.0;
    for (int i = 1; i <= 4; i++) {
        xAct = xAnt - f(xAnt) / df(xAnt);
        eps4 = errorPorcentual(xAct, xAnt);   // guarda el de la 4ta iteracion

        printf("  %3d  |   %.10lf   |   %10.6lf   %10.6lf   |   %12.8lf\n",
               i, xAct, f(xAct), df(xAct), eps4);

        xAnt = xAct;
    }

    printf("\nRESULTADOS DE NEWTON-RAPHSON (4ta iteracion):\n");
    printf("  Raiz obtenida (10 cifras decimales):  %.10lf\n", xAct);
    printf("  Error estimado (criterio #):          %.6lf %%\n", eps4);
    printf("  Error porcentual exacto:              %.6lf %%\n",
           fabs(xAct - raizExacta) / raizExacta * 100.0);
    printf("  (Raiz exacta: x = 2)\n");

    return 0;
}