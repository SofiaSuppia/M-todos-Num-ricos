#include <cstdio>
#include <math.h>

// ============================================================
// PARCIAL 1 - 2024 - PROBLEMA 2
// Datos de la misma funcion del problema 1: f(x) = ln(x^2+1) - sin(x)
//   x      |  1.0   |  1.2  |  1.5  |  1.75 |  2.0
//   y      | -0.148 | -0.040|  0.181|  0.419|  0.700
// ============================================================

#define MAX_PUNTOS 50

// Datos de la tabla
double X[] = {1.0, 1.2, 1.5, 1.75, 2.0};
double Y[] = {-0.148, -0.040, 0.181, 0.419, 0.700};
int n = 5;

// Funcion real (problema 1)
double f(double x) {
    return log(x * x + 1) - sin(x);
}

// ============================================================
// VALOR DEL POLINOMIO DE LAGRANGE
// Pn(x) = SUM y[k] * PROD(i!=k) (x-X[i])/(X[k]-X[i])
// ============================================================
double lagrangeValor(double x) {
    double suma = 0.0;
    for (int k = 0; k < n; k++) {
        double producto = 1.0;
        for (int i = 0; i < n; i++) {
            if (i != k)
                producto *= (x - X[i]) / (X[k] - X[i]);
        }
        suma += Y[k] * producto;
    }
    return suma;
}

// ============================================================
// COEFICIENTES DEL POLINOMIO DE LAGRANGE EN FORMA EXPANDIDA
// Pn(x) = coef[0] + coef[1]*x + coef[2]*x^2 + ...
// ============================================================
void lagrangeCoeficientes(double coef[]) {
    for (int i = 0; i < n; i++) coef[i] = 0.0;

    for (int k = 0; k < n; k++) {
        double base[MAX_PUNTOS] = {1.0};
        int grado = 0;

        // Denominador de L_k: PROD (X[k]-X[i]) con i != k
        double denominador = 1.0;
        for (int i = 0; i < n; i++) {
            if (i == k) continue;

            // base(x) = base(x) * (x - X[i])
            double nuevaBase[MAX_PUNTOS] = {0};
            for (int j = 0; j <= grado; j++) {
                nuevaBase[j] += base[j] * (-X[i]);
                nuevaBase[j + 1] += base[j];
            }
            for (int j = 0; j <= grado + 1; j++)
                base[j] = nuevaBase[j];

            grado++;
            denominador *= (X[k] - X[i]);
        }

        double factor = Y[k] / denominador;
        for (int j = 0; j <= grado; j++)
            coef[j] += base[j] * factor;
    }
}

// ============================================================
// MOSTRAR POLINOMIO EN FORMA LEGIBLE
// ============================================================
void mostrarPolinomio(double coef[], int grado) {
    printf("P_%d(x) = ", grado);
    int primero = 1;
    for (int i = 0; i <= grado; i++) {
        if (fabs(coef[i]) < 1e-12) continue;
        if (!primero) {
            if (coef[i] >= 0) printf(" + ");
            else printf(" - ");
        } else if (coef[i] < 0) {
            printf("-");
        }
        printf("%.8f", fabs(coef[i]));
        if (i > 0) printf("*x^%d", i);
        primero = 0;
    }
    printf("\n");
}

// ============================================================
// ELIMINACION GAUSSIANA CON PIVOTEO (para las ecuaciones
// normales del ajuste cuadratico: sistema 3x3)
// ============================================================
void gauss(int m, double A[][3], double b[], double sol[]) {
    // Triangulacion
    for (int i = 0; i < m - 1; i++) {
        // Pivoteo parcial
        int p = i;
        if (fabs(A[i][i]) < 1e-5) {
            for (int l = i + 1; l < m; l++)
                if (fabs(A[l][i]) > fabs(A[p][i])) p = l;
            for (int mm = i; mm < m; mm++) {
                double aux = A[p][mm];
                A[p][mm] = A[i][mm];
                A[i][mm] = aux;
            }
            double aux = b[p];
            b[p] = b[i];
            b[i] = aux;
        }
        // Eliminacion
        for (int j = i + 1; j < m; j++) {
            double factor = A[j][i] / A[i][i];
            for (int k = i; k < m; k++)
                A[j][k] -= factor * A[i][k];
            b[j] -= factor * b[i];
        }
    }
    // Retrosustitucion
    sol[m - 1] = b[m - 1] / A[m - 1][m - 1];
    for (int i = m - 2; i >= 0; i--) {
        double suma = b[i];
        for (int j = i + 1; j < m; j++)
            suma -= A[i][j] * sol[j];
        sol[i] = suma / A[i][i];
    }
}

int main() {
    printf("=============================================================\n");
    printf("  PARCIAL 1 - 2024 - PROBLEMA 2\n");
    printf("=============================================================\n\n");

    // --- Tabla de datos ---
    printf("Datos de la tabla:\n");
    printf("   i  |      Xi      |      Yi\n");
    printf("------|-------------|------------\n");
    for (int i = 0; i < n; i++)
        printf("%4d  | %12.4f | %12.3f\n", i + 1, X[i], Y[i]);
    printf("\n");

    // =============================================================
    // PARTE a) POLINOMIO INTERPOLADOR
    // =============================================================
    printf("--- PARTE a) POLINOMIO INTERPOLADOR ---\n");
    double coef[MAX_PUNTOS];
    lagrangeCoeficientes(coef);
    mostrarPolinomio(coef, n - 1); // P4(x)
    printf("    Verificacion en los nodos: P(xi) = yi\n");
    for (int i = 0; i < n; i++)
        printf("      P(%.2f) = %.6f  (tabla: %.3f)\n", X[i], lagrangeValor(X[i]), Y[i]);
    printf("\n");

    // =============================================================
    // PARTE b) APROXIMACION DE f(1.6) Y ERROR EXACTO
    // =============================================================
    printf("--- PARTE b) VALOR APROXIMADO DE f(1.6) ---\n");
    double xs = 1.6;
    double p16 = lagrangeValor(xs);
    double real16 = f(xs);
    printf("    P(1.6) = %.10f\n", p16);
    printf("    f(1.6) = %.10f  (verdadero)\n", real16);
    printf("    Error absoluto exacto = |P(1.6) - f(1.6)| = %.10f\n\n",
           fabs(p16 - real16));

    // =============================================================
    // PARTE c) GRAFICO DE E(x)=|P(x)-f(x)|
    // =============================================================
    printf("--- PARTE c) FUNCION ERROR E(x) = |P(x) - f(x)| ---\n");
    printf("    Se grafica con problema2_graficas.py (guarda el PNG).\n");
    printf("    E(x_i) = |P(x_i) - f(x_i)| = |y_i - f(x_i)|:\n");
    for (int i = 0; i < n; i++)
        printf("      E(%.2f) = %.6f\n", X[i], fabs(Y[i] - f(X[i])));
    printf("    E(x) vale 0 en los nodos: el polinomio interpola exacto\n");
    printf("    y crece entre nodo y nodo (maximo entre dos nodos).\n\n");

    // =============================================================
    // PARTE d) REGRESION LINEAL
    // =============================================================
    printf("--- PARTE d) REGRESION LINEAL  y = a0 + a1*x ---\n");
    double sx = 0, sy = 0, sxx = 0, sxy = 0, syy = 0;
    for (int i = 0; i < n; i++) {
        sx += X[i];        sy += Y[i];
        sxx += X[i] * X[i]; sxy += X[i] * Y[i];
        syy += Y[i] * Y[i];
    }
    // Minimos cuadrados (formula cerrada)
    double a1 = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    double a0 = (sy - a1 * sx) / n;

    // Coeficiente de correlacion r = sqrt((St - Sr)/St)
    double ym = sy / n;
    double St = 0, Sr = 0, ajuste;
    for (int i = 0; i < n; i++) {
        St += pow(Y[i] - ym, 2);
        ajuste = a0 + a1 * X[i];
        Sr += pow(Y[i] - ajuste, 2);
    }
    double r = sqrt((St - Sr) / St);

    // Raiz de la recta: y = 0 -> x = -a0/a1
    double raiz_lin = -a0 / a1;
    double alpha = 1.2620326547; // valor exacto de la raiz (problema 1)

    printf("    a0 = %.10f\n", a0);
    printf("    a1 = %.10f\n", a1);
    printf("    y = %.6f + %.6f*x\n", a0, a1);
    printf("    St = %.10f , Sr = %.10f\n", St, Sr);
    printf("    r  = sqrt((St-Sr)/St) = %.10f\n", r);
    printf("    R2 = %.10f\n", r * r);
    if (fabs(r - 1) < 0.1)
        printf("    -> r muy cercano a 1: ajuste MUY bueno\n");
    else if (fabs(r - 1) < 0.25)
        printf("    -> ajuste bueno\n");
    else
        printf("    -> ajuste pobre\n");
    printf("    Raiz de la recta (y=0): x = %.10f\n", raiz_lin);
    printf("    Error |raiz_lin - alpha| = %.10f\n\n", fabs(raiz_lin - alpha));

    // =============================================================
    // PARTE e) REGRESION CUADRATICA
    // =============================================================
    printf("--- PARTE e) REGRESION CUADRATICA  y = a0 + a1*x + a2*x^2 ---\n");
    double sx3 = 0, sx4 = 0, sx2y = 0;
    for (int i = 0; i < n; i++) {
        sx3 += pow(X[i], 3);
        sx4 += pow(X[i], 4);
        sx2y += X[i] * X[i] * Y[i];
    }

    // Ecuaciones normales:  A * a = b  (sistema 3x3)
    // [ n    Sx    Sx2 ] [a0]   [ Sy ]
    // [ Sx   Sx2   Sx3 ] [a1] = [Sxy]
    // [ Sx2  Sx3   Sx4 ] [a2]   [Sx2y]
    double A[3][3] = {
        { (double)n, sx, sxx },
        { sx,         sxx, sx3 },
        { sxx,        sx3, sx4 }
    };
    double b[3] = { sy, sxy, sx2y };
    double a[3];

    gauss(3, A, b, a);

    // Coeficiente de correlacion del ajuste cuadratico
    double Sr2 = 0;
    for (int i = 0; i < n; i++) {
        ajuste = a[0] + a[1] * X[i] + a[2] * X[i] * X[i];
        Sr2 += pow(Y[i] - ajuste, 2);
    }
    double r2 = sqrt((St - Sr2) / St);

    // Raices de la parabola a0 + a1*x + a2*x^2 = 0
    double disc = a[1] * a[1] - 4 * a[0] * a[2];

    printf("    a0 = %.10f\n", a[0]);
    printf("    a1 = %.10f\n", a[1]);
    printf("    a2 = %.10f\n", a[2]);
    printf("    y = %.6f + %.6f*x + %.6f*x^2\n", a[0], a[1], a[2]);
    printf("    Sr = %.10f\n", Sr2);
    printf("    r  = sqrt((St-Sr)/St) = %.10f\n", r2);
    printf("    R2 = %.10f\n", r2 * r2);

    if (disc >= 0) {
        double r1 = (-a[1] - sqrt(disc)) / (2 * a[2]);
        double r2r = (-a[1] + sqrt(disc)) / (2 * a[2]);
        printf("    Raices: x1 = %.10f , x2 = %.10f\n", r1, r2r);
        double e1 = fabs(r1 - alpha), e2 = fabs(r2r - alpha);
        printf("    |x1 - alpha| = %.10f , |x2 - alpha| = %.10f\n", e1, e2);
        if (e1 < e2)
            printf("    -> Raiz que se busca = %.10f , error = %.10f\n", r1, e1);
        else
            printf("    -> Raiz que se busca = %.10f , error = %.10f\n", r2r, e2);
    } else {
        printf("    La parabola no tiene raices reales (discriminante < 0).\n");
    }

    printf("\n--- COMPARACION DE AJUSTES ---\n");
    printf("    Lineal:    r  = %.10f , Sr = %.10f\n", r, Sr);
    printf("    Cuadratico: r = %.10f , Sr = %.10f\n", r2, Sr2);
    if (r2 > r)
        printf("    -> El ajuste CUADRATICO es mejor (r mas cercano a 1 y Sr menor).\n");
    else
        printf("    -> El ajuste LINEAL es mejor.\n");

    return 0;
}