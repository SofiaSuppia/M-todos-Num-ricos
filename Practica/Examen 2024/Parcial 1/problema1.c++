#include <cstdio>
#include <math.h>

// ============================================================
// PARCIAL 1 - 2024 - PROBLEMA 1
// f(x) = ln(x^2 + 1) - sin(x)   (raiz en [1.0 , 2.0])
// ============================================================

// f(x) = ln(x^2 + 1) - sin(x)
double f(double x) {
    return log(x * x + 1) - sin(x);
}

// f'(x) = 2x/(x^2 + 1) - cos(x) (derivada analitica, para Newton-Raphson)
double df(double x) {
    return 2.0 * x / (x * x + 1) - cos(x);
}

// g(x) elegida para el metodo de punto fijo:
// f(x)=0 -> ln(x^2+1) = sin(x) -> x = sqrt(e^{sin(x)} - 1)
double g(double x) {
    return sqrt(exp(sin(x)) - 1);
}

// g'(x) = e^{sin(x)}*cos(x) / (2*sqrt(e^{sin(x)} - 1))
double dg(double x) {
    return exp(sin(x)) * cos(x) / (2.0 * sqrt(exp(sin(x)) - 1));
}

int main() {
    printf("=============================================================\n");
    printf("  PARCIAL 1 - 2024 - PROBLEMA 1\n");
    printf("  f(x) = ln(x^2 + 1) - sin(x)\n");
    printf("=============================================================\n\n");

    // =============================================================
    // PARTE a) Existencia de raiz en [1.0, 2.0] (Bolzano)
    // =============================================================
    double aE = 1.0, bE = 2.0;
    printf("--- PARTE a) ---\n");
    printf("f(1.0) = %.10f\n", f(aE));
    printf("f(2.0) = %.10f\n", f(bE));
    printf("f(1.0)*f(2.0) = %.10f  <  0  ->  existe raiz en [1,2] (Bolzano)\n\n",
           f(aE) * f(bE));

    // =============================================================
    // PARTE b) BISECCION
    // =============================================================
    printf("--- PARTE b) BISECCION ---\n\n");

    // b1) Primera iteracion a mano, [1, 2]
    double a1 = 1.0, b1 = 2.0;
    double c1 = (a1 + b1) / 2.0;
    printf("b1) Primera iteracion:\n");
    printf("    [a,b] = [1.0000000000, 2.0000000000]\n");
    printf("    c1 = (a+b)/2 = %.10f\n", c1);
    printf("    f(a) = %.10f , f(b) = %.10f\n", f(a1), f(b1));
    printf("    f(c1) = %.10f\n", f(c1));
    if (f(a1) * f(c1) < 0)
        printf("    f(a)*f(c1) < 0  ->  nueva raiz en [a, c1] = [1.0000000000, %.10f]\n", c1);
    else
        printf("    f(a)*f(c1) > 0  ->  nueva raiz en [c1, b]\n");
    printf("    Error absoluto estimado: |b - a|/2 = %.10f\n\n", (b1 - a1) / 2.0);

    // b2) Iteracion numero 10
    printf("b2) Diez iteraciones de biseccion:\n");
    double aB = 1.0, bB = 2.0, cB;
    for (int i = 1; i <= 10; i++) {
        cB = (aB + bB) / 2.0;
        if (f(aB) * f(cB) < 0)
            bB = cB;
        else
            aB = cB;
        printf("    it %2d: raiz = %.10f   error = %.10f\n", i, cB, (bB - aB));
    }
    printf("    -> Iteracion 10: raiz = %.10f, error absoluto = %.10f\n\n", cB, (bB - aB));

    // =============================================================
    // PARTE c) REGULA FALSI
    // =============================================================
    printf("--- PARTE c) REGULA FALSI ---\n\n");

    // c1) Primera iteracion a mano, [1, 2]
    double aR = 1.0, bR = 2.0;
    double cR = (aR * f(bR) - bR * f(aR)) / (f(bR) - f(aR));
    printf("c1) Primera iteracion:\n");
    printf("    [a,b] = [1.0000000000, 2.0000000000]\n");
    printf("    c1 = (a*f(b) - b*f(a)) / (f(b) - f(a))\n");
    printf("       = (1*%.5f - 2*(%.5f)) / (%.5f - (%.5f))\n", f(bR), f(aR), f(bR), f(aR));
    printf("       = %.10f\n", cR);
    printf("    f(c1) = %.10f\n", f(cR));
    if (f(aR) * f(cR) < 0)
        printf("    f(a)*f(c1) < 0  ->  nueva raiz en [a, c1] = [1.0000000000, %.10f]\n", cR);
    else
        printf("    f(a)*f(c1) > 0  ->  nueva raiz en [c1, b] = [%.10f, 2.0000000000]\n", cR);
    printf("    Error absoluto estimado: |c1 - a| = %.10f\n\n", fabs(cR - aR));

    // c2) Iteracion numero 10
    printf("c2) Diez iteraciones de Regula Falsi:\n");
    double cr = 0.0, cr_viejo = aR, errR = 0.0;
    aR = 1.0; bR = 2.0;
    for (int i = 1; i <= 10; i++) {
        cr = (aR * f(bR) - bR * f(aR)) / (f(bR) - f(aR));
        if (f(aR) * f(cr) < 0)
            bR = cr;
        else
            aR = cr;
        errR = fabs(cr - cr_viejo);
        printf("    it %2d: raiz = %.10f   error = %.10f\n", i, cr, errR);
        cr_viejo = cr;
    }
    printf("    -> Iteracion 10: raiz = %.10f, error absoluto = %.10f\n\n", cr, errR);

    // =============================================================
    // PARTE d) METODO DE PUNTO FIJO
    // =============================================================
    printf("--- PARTE d) PUNTO FIJO ---\n");
    double alpha = 1.26203265474; // raiz aproximada
    printf("    f(x)=0  =>  x^2 + 1 = e^{sin(x)}  =>  x = sqrt(e^{sin(x)} - 1) = g(x)\n");
    printf("    g'(x) = e^{sin(x)}*cos(x) / (2*sqrt(e^{sin(x)} - 1))\n");
    printf("    |g'(alpha)| = %.10f  <  1  ->  SI converge (|g'| < 1)\n\n", fabs(dg(alpha)));

    // =============================================================
    // PARTE e) NEWTON - RAPHSON
    // =============================================================
    printf("--- PARTE e) NEWTON-RAPHSON ---\n\n");

    // e1) Primera iteracion con x0 = 1.0
    double x0 = 1.0;
    double x1 = x0 - f(x0) / df(x0);
    printf("e1) Primera iteracion (x0 = 1.0):\n");
    printf("    f(x0) = %.10f\n", f(x0));
    printf("    f'(x0) = 2*x0/(x0^2+1) - cos(x0) = %.10f\n", df(x0));
    printf("    x1 = x0 - f(x0)/f'(x0) = %.10f\n", x1);
    printf("    Error absoluto estimado: |x1 - x0| = %.10f\n\n", fabs(x1 - x0));

    // e2) Cuarta iteracion
    printf("e2) Cuatro iteraciones de Newton-Raphson:\n");
    double x = 1.0, xn, err;
    for (int i = 1; i <= 4; i++) {
        xn = x - f(x) / df(x);
        err = fabs(xn - x);
        printf("    it %d: raiz = %.10f   error = %.10f\n", i, xn, err);
        x = xn;
    }
    printf("    -> Iteracion 4: raiz = %.10f, error absoluto = %.10f\n", xn, err);

    return 0;
}