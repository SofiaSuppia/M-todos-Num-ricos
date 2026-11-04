// ============================================================================
//  ENUNCIADO  (Problema n.° 4 — Examen 2025)
//
//  Mediante un radar Doppler se mide la velocidad de caida de un paracaidista
//  como funcion del tiempo. Los datos obtenidos en una prueba particular se
//  muestran en la tabla.
//
//     Tiempo (s)   Velocidad medida (cm/s)
//        1                800
//        3               2310
//        5               3090
//        7               3940
//       13               4755
//
//  a) Empleando interpolacion de Lagrange, estime las velocidades para
//     t = 4 s y t = 6 s. (Adjunte captura de pantalla).
//
//  b) Obtenga el polinomio de Lagrange que interpola los datos de la tabla.
//
//  c) Use interpolacion cubica (normal) para estimar la velocidad en
//     t = 10 s, y escriba el polinomio usado para calcular dicha estimacion.
//     Escriba ambas cosas en papel, y adjunte captura de pantalla de la consola.
// ============================================================================
#include <cstdio>
#include <cmath>
#include <cstdlib>

#define MAX_PUNTOS 50   // cantidad maxima de puntos de datos

// ============================================================
// SECCION: LECTURA DE ARCHIVO
// Lee los pares (x, y) desde 'data.txt' (mismo criterio que
// 3. Interpolacion/LagrangePolinomio.c++): una linea por punto "xi yi",
// sin encabezados. La cantidad de puntos se cuenta por lineas.
// ============================================================
bool leerDatosPuntos(const char* nombreArchivo, double X[], double Y[], int* n) {
    FILE* archivo = fopen(nombreArchivo, "r");
    if (archivo == NULL) {
        printf("ERROR: no se puede abrir el archivo '%s'\n", nombreArchivo);
        return false;
    }
    printf("Archivo '%s' abierto correctamente\n\n", nombreArchivo);

    int filas = 0;
    char buffer[1024];
    while (fgets(buffer, sizeof(buffer), archivo) != NULL) filas++;
    *n = filas;
    printf("Cantidad de puntos: %d\n", *n);

    if (*n <= 0 || *n > MAX_PUNTOS) {
        printf("ERROR: cantidad de puntos invalida (%d)\n", *n);
        fclose(archivo);
        return false;
    }

    fclose(archivo);
    archivo = fopen(nombreArchivo, "r");

    for (int i = 0; i < *n; i++) {
        if (fscanf(archivo, "%lf %lf", &X[i], &Y[i]) != 2) {
            printf("ERROR: no se pudo leer el punto numero %d\n", i + 1);
            fclose(archivo);
            return false;
        }
    }
    printf("Se leyeron %d puntos de datos correctamente\n\n", *n);
    fclose(archivo);
    return true;
}

// ============================================================
// SECCION: MOSTRAR PUNTOS DE DATOS
// ============================================================
void mostrarPuntosDatos(double X[], double Y[], int n) {
    printf("Puntos de datos:\n");
    printf("=================\n");
    printf("   i  |      Xi      |      Yi\n");
    printf("------|--------------|-----------\n");
    for (int i = 0; i < n; i++) {
        printf("%4d  | %12.6f | %12.6f\n", i + 1, X[i], Y[i]);
    }
    printf("\n");
}

// ============================================================
// SECCION: POLINOMIO DE LAGRANGE (valor numerico)
//   Pn(x) = Σ(k) y[k]·Lk(x)
//   Lk(x) = Π(i!=k) (x - X[i]) / (X[k] - X[i])
// Imprime cada coeficiente Lk(x) (la suma de todos da 1) y devuelve Pn(x).
// ============================================================
double lagrangeValor(double x, double X[], double Y[], int n) {
    double suma = 0.0;
    for (int k = 0; k < n; k++) {
        double producto = 1.0;
        for (int i = 0; i < n; i++) {
            if (i != k) producto = producto * ((x - X[i]) / (X[k] - X[i]));
        }
        printf("L%d%d(%.3f) = %.6f\n", n - 1, k, x, producto);
        suma = suma + Y[k] * producto;
    }
    return suma;
}

// ============================================================
// SECCION: COEFICIENTES DEL POLINOMIO DE LAGRANGE (expandido)
// Convierte la forma con productos Lk(x) en la forma clasica
//   P(x) = coef[0] + coef[1]*x + coef[2]*x^2 + ...
// (misma tecnica que lagrangeCoeficientes de 3. Interpolacion).
// ============================================================
void lagrangeCoeficientes(double X[], double Y[], int n, double coef[]) {
    for (int i = 0; i < n; i++) coef[i] = 0.0;

    for (int k = 0; k < n; k++) {
        double base[MAX_PUNTOS] = {1.0};
        int grado = 0;
        double denominador = 1.0;

        for (int i = 0; i < n; i++) {
            if (i == k) continue;
            double nuevaBase[MAX_PUNTOS] = {0};
            for (int j = 0; j <= grado; j++) {
                nuevaBase[j] += base[j] * (-X[i]);
                nuevaBase[j + 1] += base[j];
            }
            for (int j = 0; j <= grado + 1; j++) base[j] = nuevaBase[j];
            grado++;
            denominador *= (X[k] - X[i]);
        }

        double factor = Y[k] / denominador;
        for (int j = 0; j <= grado; j++) coef[j] += base[j] * factor;
    }
}

// ============================================================
// SECCION: MOSTRAR POLINOMIO EXPANDIDO
// ============================================================
void mostrarPolinomio(double coef[], int n, const char* nombre) {
    printf("%s = ", nombre);
    int primero = 1;
    for (int i = 0; i < n; i++) {
        if (fabs(coef[i]) < 1e-9) continue;
        if (!primero) {
            if (coef[i] >= 0) printf(" + ");
            else printf(" - ");
        } else if (coef[i] < 0) {
            printf("-");
        }
        printf("%.6f", fabs(coef[i]));
        if (i > 0) printf("*t^%d", i);
        primero = 0;
    }
    printf("\n\n");
}

// ============================================================
// SECCION: ELECCION DE LOS NODOS MAS CERCANOS (para la cubica)
// Dado x, elige los m nodos de la tabla mas cercanos en distancia.
// Para t = 10 (entre 7 y 13) los 4 mas cercanos resultan 3, 5, 7 y 13.
// ============================================================
void elegirNodosCercanos(double x, double X[], int n, int m, int idx[]) {
    for (int i = 0; i < m; i++) idx[i] = i;
    for (int i = m; i < n; i++) {
        int lejos = 0;
        for (int j = 1; j < m; j++) {
            if (fabs(X[idx[j]] - x) > fabs(X[idx[lejos]] - x)) lejos = j;
        }
        if (fabs(X[i] - x) < fabs(X[idx[lejos]] - x)) idx[lejos] = i;
    }
}

// ============================================================
// PROGRAMA PRINCIPAL
// ============================================================
int main() {
    int n;
    double X[MAX_PUNTOS], Y[MAX_PUNTOS];

    printf("=========================================================\n");
    printf(" PROBLEMA 4 - RADAR DOPPLER (INTERPOLACION DE LAGRANGE)\n");
    printf("=========================================================\n\n");

    if (!leerDatosPuntos("data.txt", X, Y, &n)) {
        printf("No se pudo leer el archivo de datos. Fin del programa.\n");
        return 1;
    }
    mostrarPuntosDatos(X, Y, n);

    // ============================================================
    // PARTE a) LAGRANGE: estimar v(4) y v(6) con todos los nodos
    // ============================================================
    double puntos_a[] = {4.0, 6.0};

    printf("============== PARTE a) INTERPOLACION DE LAGRANGE ==============\n\n");
    for (int p = 0; p < 2; p++) {
        double t = puntos_a[p];
        printf("Estimar v(t) para t = %.1f s (polinomio de grado 4, 5 nodos):\n", t);
        double v = lagrangeValor(t, X, Y, n);
        printf("\nv(%.1f) = 800*L40 + 2310*L41 + 3090*L42 + 3940*L43 + 4755*L44 = %.6f\n\n",
               t, v);
    }

    // ============================================================
    // PARTE b) POLINOMIO DE LAGRANGE QUE INTERPOLA LA TABLA (grado 4)
    // ============================================================
    double coef4[MAX_PUNTOS];

    printf("============== PARTE b) POLINOMIO DE LAGRANGE (GRADO 4) ==============\n\n");
    lagrangeCoeficientes(X, Y, n, coef4);
    mostrarPolinomio(coef4, n, "P4(t)");

    printf("Verificacion: el polinomio debe reproducir los datos (P4(ti) = vi):\n");
    for (int i = 0; i < n; i++) {
        double valor = 0.0;
        for (int j = 0; j < n; j++) valor += coef4[j] * pow(X[i], j);
        printf("  P4(%5.1f) = %10.6f   (dato: %.1f)\n", X[i], valor, Y[i]);
    }

    // ============================================================
    // PARTE c) INTERPOLACION CUBICA para estimar v(10)
    // Se usan los 4 nodos mas cercanos a t = 10: 3, 5, 7 y 13 s.
    // ============================================================
    int idx[4];
    double Xc[4], Yc[4];
    double t_c = 10.0;

    printf("\n============== PARTE c) INTERPOLACION CUBICA EN t = 10 s ==============\n\n");
    elegirNodosCercanos(t_c, X, n, 4, idx);

    printf("t = %.0f s esta entre 7 y 13 s. Se toman los 4 nodos mas cercanos:\n", t_c);
    printf("  nodo de t = %.0f s  -> v = %.0f cm/s\n", X[idx[0]], Y[idx[0]]);
    printf("  nodo de t = %.0f s  -> v = %.0f cm/s\n", X[idx[1]], Y[idx[1]]);
    printf("  nodo de t = %.0f s  -> v = %.0f cm/s\n", X[idx[2]], Y[idx[2]]);
    printf("  nodo de t = %.0f s  -> v = %.0f cm/s\n\n", X[idx[3]], Y[idx[3]]);

    for (int i = 0; i < 4; i++) {
        Xc[i] = X[idx[i]];
        Yc[i] = Y[idx[i]];
    }

    printf("Polinomio cubico de Lagrange con esos nodos:\n");
    double v_cubic = lagrangeValor(t_c, Xc, Yc, 4);
    printf("\nv(10) = 2310*L30 + 3090*L31 + 3940*L32 + 4755*L33 = %.6f cm/s\n\n",
           v_cubic);

    double coef3[MAX_PUNTOS];
    lagrangeCoeficientes(Xc, Yc, 4, coef3);
    printf("Polinomio cubico usado (grado 3, forma clasica):\n");
    mostrarPolinomio(coef3, 4, "P3(t)");

    return 0;
}