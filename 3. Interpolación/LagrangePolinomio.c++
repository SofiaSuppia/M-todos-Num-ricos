#include <cstdio>
#include <cmath>
#include <iostream>

#define MAX_PUNTOS 50   // cantidad maxima de puntos de datos
#define MAX_TAM 100     // tamanio maximo para la matriz del sistema

using namespace std;

// ============================================================
// SECCION: LECTURA DE ARCHIVO
// Funcion para leer los pares de datos (x, y) desde un archivo.
// Formato esperado del archivo 'data.txt' (igual criterio que la
// lectura de matrices en EliminGaussiana):
//   Una linea por punto: xi yi
//
// Ejemplo para 3 puntos:
//   1.0 3.0
//   2.0 3.0
//   2.5 3.3
//
// La cantidad de puntos se obtiene contando las lineas del archivo.
//
// @param nombreArchivo  nombre del archivo a leer
// @param X              arreglo donde se guardan las x
// @param Y              arreglo donde se guardan las y
// @param n              puntero donde se guarda la cantidad de puntos
// @return               true si se leyó correctamente, false si hubo error
// ============================================================
bool leerDatosPuntos(const char* nombreArchivo, double X[], double Y[], int* n) {
    // Se abre el archivo en modo lectura
    FILE* archivo;
    archivo = fopen(nombreArchivo, "r");

    // Si no se pudo abrir, se avisa y se termina
    if (archivo == NULL) {
        printf("ERROR: no se puede abrir el archivo '%s'\n", nombreArchivo);
        return false;
    }

    printf("Archivo '%s' abierto correctamente\n\n", nombreArchivo);

    // Se cuentan las lineas del archivo = cantidad de puntos
    // (mismo criterio que la lectura de matrices en EliminGaussiana)
    int filas = 0;
    char buffer[1024];
    while (fgets(buffer, sizeof(buffer), archivo) != NULL) {
        filas++;
    }

    // La cantidad de puntos es la cantidad de lineas leidas
    *n = filas;
    printf("Cantidad de puntos: %d\n", *n);

    // Se verifica que la cantidad de puntos sea valida
    if (*n <= 0 || *n > MAX_PUNTOS) {
        printf("ERROR: cantidad de puntos invalida (%d)\n", *n);
        fclose(archivo);
        return false;
    }

    // Se cierra y reabre el archivo para reiniciar el puntero de lectura
    // (mismo criterio que la lectura de matrices en EliminGaussiana)
    fclose(archivo);
    archivo = fopen(nombreArchivo, "r");

    // Se leen los pares (x, y) uno por uno
    for (int i = 0; i < *n; i++) {
        if (fscanf(archivo, "%lf %lf", &X[i], &Y[i]) != 2) {
            printf("ERROR: no se pudo leer el punto numero %d\n", i + 1);
            fclose(archivo);
            return false;
        }
    }

    // Se cierra el archivo y se confirma la lectura
    fclose(archivo);
    printf("Se leyeron %d puntos de datos correctamente\n\n", *n);
    return true;
}

// ============================================================
// SECCION: MOSTRAR PUNTOS DE DATOS
// Funcion que imprime en pantalla los pares (x, y) leídos
// para que el usuario verifique que se cargaron bien.
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
// SECCION: POLINOMIO DE LAGRANGE (valor numérico)
// Calcula el valor interpolado Pn(xsombrero) usando la formula
// de Lagrange:
//
//   Pn(x) = SUMATORIA(k=0..n)  y[k] * PROD(i=0..n, i!=k)
//            (x - x[i]) / (x[k] - x[i])
//
// @param x   el punto donde se quiere interpolar (xsombrero)
// @param X   arreglo de abscisas
// @param Y   arreglo de ordenadas
// @param n   cantidad de puntos
// @return    el valor interpolado Pn(x)
// ============================================================
double lagrangeValor(double x, double X[], double Y[], int n) {
    double suma = 0.0;   // acumula el resultado de Pn(x)

    // Lazo externo: recorre cada termino k del polinomio
    for (int k = 0; k < n; k++) {
        double producto = 1.0;  // acumula el producto de Lagrange

        // Lazo interno: calcula el producto Lk(x) = PROD (x-x[i])/(x[k]-x[i])
        for (int i = 0; i < n; i++) {
            // Se saltea el termino i == k (no se puede dividir por cero)
            if (i != k) {
                // Se multiplica cada factor (x - x[i]) / (x[k] - x[i])
                producto = producto * ((x - X[i]) / (X[k] - X[i]));
            }
        }

        // Se muestra el coeficiente Cnk(x) para verificar que suman 1
        printf("L%d%d(%.3f) = %.6f\n", n - 1, k, x, producto);

        // Se acumula el termino y[k] * Lk(x)
        suma = suma + Y[k] * producto;
    }

    // Se devuelve el valor interpolado
    return suma;
}

// ============================================================
// SECCION: PIVOTEO PARA ELIMINACION GAUSSIANA
// Si el elemento de la diagonal a[i][i] es muy chico (casi cero),
// se busca en la columna i la fila con el mayor valor absoluto
// y se intercambian las filas completas (matriz y vector b).
//
// @param n   tamanio del sistema
// @param A   matriz de coeficientes
// @param b   vector de terminos independientes
// @param i   columna actual (donde esta el pivote)
// ============================================================
void pivoteo(int n, double A[][MAX_TAM], double b[], int i) {
    int filaMayor = i;    // fila que tiene el mejor pivote
    double auxiliar;      // variable para intercambiar   valores

    // Se busca la fila con el elemento de mayor valor absoluto
    // en la columna i, desde la fila i+1 hasta la ultima
    for (int l = i + 1; l < n; l++) {
        if (fabs(A[l][i]) > fabs(A[filaMayor][i])) {
            filaMayor = l;   // se guarda la mejor fila encontrada
        }
    }

    // Si se encontro una mejor fila que la actual, se intercambian
    if (filaMayor != i) {
        // Se intercambian las filas i y filaMayor en la matriz A
        for (int m = i; m < n; m++) {
            auxiliar = A[filaMayor][m];
            A[filaMayor][m] = A[i][m];
            A[i][m] = auxiliar;
        }

        // Se intercambian las filas i y filaMayor en el vector b
        auxiliar = b[filaMayor];
        b[filaMayor] = b[i];
        b[i] = auxiliar;
    }
}

// ============================================================
// SECCION: ELIMINACION GAUSSIANA CON PIVOTEO
// Resuelve el sistema de ecuaciones lineales A * x = b
// mediante el metodo de eliminacion gaussiana con pivoteo.
//
// Pasos:
//   1. Pivotear (intercambiar filas si el pivote es muy chico)
//   2. Triangulacion (hacer ceros debajo de la diagonal)
//   3. Retrosustitucion (resolver de abajo hacia arriba)
//   4. Calcular el determinante (producto de la diagonal)
//
// @param n         tamanio del sistema
// @param A         matriz de coeficientes (se modifica)
// @param b         vector de terminos independientes (se modifica)
// @param solucion  vector donde se guarda el resultado x
// ============================================================
void gaussPivoteo(int n, double A[][MAX_TAM], double b[], double solucion[]) {
    double factor;    // multiplicador para anular elementos bajo el pivote
    double suma;      // acumulador para la retrosustitucion
    double det = 1.0; // determinante de la matriz

    // --- PASO 1 Y 2: PIVOTEO + TRIANGULACION ---
    // Se recorre cada columna (pivote)
    for (int i = 0; i < n - 1; i++) {

        // PRIMERO SE PIVOTEA: se evita dividir por cero o numeros muy chicos
        if (fabs(A[i][i]) < 1e-5) {
            pivoteo(n, A, b, i);
        }

        // LUEGO SE ELIMINA: se hacen ceros los elementos debajo del pivote
        for (int j = i + 1; j < n; j++) {
            // factor = a[j][i] / a[i][i]  (relacion entre fila actual y pivote)
            factor = A[j][i] / A[i][i];

            // Se resta 'factor' veces la fila pivote de la fila j
            for (int k = i; k < n; k++) {
                A[j][k] = A[j][k] - factor * A[i][k];
            }
            // Se hace lo mismo con el termino independiente
            b[j] = b[j] - factor * b[i];
        }
    }

    // --- PASO 3: CALCULO DEL DETERMINANTE ---
    // El determinante de una matriz triangular es el producto
    // de los elementos de la diagonal principal.
    for (int i = 0; i < n; i++) {
        det = det * A[i][i];
    }

    printf("------------------DETERMINANTE------------------\n");
    printf("El determinante de la matriz es: %lf\n\n", det);

    // Si el determinante es cero, la matriz no tiene solucion unica
    if (fabs(det) < 1e-10) {
        printf("ERROR: el determinante es cero, el sistema no tiene solucion unica.\n");
        exit(1);
    }

    // --- PASO 4: RETROSUSTITUCION ---
    // Se resuelve de la ultima ecuacion hacia la primera.
    //   x[n-1] = b[n-1] / a[n-1][n-1]
    //   x[i]   = (b[i] - SUMA(a[i][j]*x[j], j=i+1..n)) / a[i][i]

    // Se calcula la ultima incognita directamente
    solucion[n - 1] = b[n - 1] / A[n - 1][n - 1];

    // Se recorren las restantes incognitas de abajo hacia arriba
    for (int i = n - 2; i >= 0; i--) {
        suma = b[i];  // se arranca con el termino independiente

        // Se restan los terminos que ya se conocen
        for (int j = i + 1; j < n; j++) {
            suma = suma - A[i][j] * solucion[j];
        }

        // Se divide por el elemento de la diagonal
        solucion[i] = suma / A[i][i];
    }
}

// ============================================================
// SECCION: MOSTRAR POLINOMIO INTERPOLANTE
// Imprime el polinomio con formato legible:
//   Pn(x) = a0 + a1*x + a2*x^2 + ... + a(n-1)*x^(n-1)
// ============================================================
void mostrarPolinomio(double coef[], int n) {
    printf("Pn(x) = ");

    // Se imprime el primer termino (el constante a0)
    if (fabs(coef[0]) > 1e-10) {
        printf("%.6f", coef[0]);
    } else {
        printf("0");
    }

    // Se imprimen los demas terminos
    for (int i = 1; i < n; i++) {
        if (fabs(coef[i]) > 1e-10) {
            // Se imprime el signo
            if (coef[i] > 0) {
                printf(" + ");
            } else {
                printf(" - ");
            }

            // Se imprime el coeficiente (sin el signo, que ya se imprimio)
            double coeficiente = fabs(coef[i]);
            if (coeficiente != 1.0) {
                printf("%.6f", coeficiente);
            }

            // Se imprime la variable con su exponente
            if (i == 1) {
                printf("x");
            } else {
                printf("x^%d", i);
            }
        }
    }
    printf("\n\n");
}

// ============================================================
// SECCION: COEFICIENTES DEL POLINOMIO DE LAGRANGE (expandido)
// Calcula los coeficientes del polinomio de Lagrange en su
// forma expandida:  a0 + a1*x + a2*x^2 + ...
// Es decir, convierte la forma con productos (Lagrange) en
// la forma de polinomio clasico.
// ============================================================
void lagrangeCoeficientes(double X[], double Y[], int n, double coef[]) {
    // Se inicializan todos los coeficientes en cero
    for (int i = 0; i < n; i++) {
        coef[i] = 0.0;
    }

    // Se recorre cada termino k del polinomio de Lagrange
    for (int k = 0; k < n; k++) {
        // 'base' almacena el polinomio base Lk(x) que va creciendo
        // multiplicando por (x - X[i]) para cada i != k
        double base[MAX_PUNTOS] = {1.0};
        int grado = 0;

        // Se calcula el denominador: prod (X[k] - X[i]) para i != k
        double denominador = 1.0;
        for (int i = 0; i < n; i++) {
            if (i == k) continue;   // se saltea el termino donde i == k

            // Se multiplica el polinomio base por (x - X[i]):
            //   base_j * (-X[i])  -> termino constante
            //   base_j * x        -> termino con x
            double nuevaBase[MAX_PUNTOS] = {0};
            for (int j = 0; j <= grado; j++) {
                nuevaBase[j] += base[j] * (-X[i]);
                nuevaBase[j + 1] += base[j];
            }

            // Se copia la nueva base y se aumenta el grado
            for (int j = 0; j <= grado + 1; j++) {
                base[j] = nuevaBase[j];
            }
            grado++;

            // Se multiplica el denominador
            denominador *= (X[k] - X[i]);
        }

        // factor = y[k] / denominador
        double factor = Y[k] / denominador;

        // Se acumulan los coeficientes de este termino
        for (int j = 0; j <= grado; j++) {
            coef[j] += base[j] * factor;
        }
    }
}

// ============================================================
// SECCION: MOSTRAR LAGRANGE EXPANDIDO
// Imprime el polinomio de Lagrange ya expandido (forma clasica)
// ============================================================
void mostrarLagrangeExpandido(double coef[], int n) {
    printf("P_%d(x) = ", n - 1);
    int primero = 1;   // controla si es el primer termino impreso

    for (int i = 0; i < n; i++) {
        // Se saltean los coeficientes despreciables
        if (fabs(coef[i]) < 1e-12) continue;

        // Se imprime el signo
        if (!primero) {
            if (coef[i] >= 0) printf(" + ");
            else printf(" - ");
        } else if (coef[i] < 0) {
            printf("-");
        }

        // Se imprime el coeficiente y la parte de x
        printf("%.6f", fabs(coef[i]));
        if (i > 0) printf("*x^%d", i);
        primero = 0;
    }
    printf("\n");
}

// ============================================================
// PROGRAMA PRINCIPAL
//
// Menu con las distintas formas de interpolar:
//   1. Lagrange numerico  -> da el valor Pn(xsombrero)
//   2. Polinomio interpolante -> arma el sistema de Vandermonde
//      y lo resuelve con gaussPivoteo (pseudocodigo 2 de la profe)
//   3. Lagrange expandido -> muestra el polinomio en forma clasica
// ============================================================
int main(int argc, char const *argv[]) {

    // --- LECTURA DE LOS DATOS DESDE EL ARCHIVO ---
    int n;                      // cantidad de puntos
    double X[MAX_PUNTOS];       // arreglo de abscisas xi
    double Y[MAX_PUNTOS];       // arreglo de ordenadas yi

    // Se lee la data desde el archivo 'data.txt'
    if (!leerDatosPuntos("data.txt", X, Y, &n)) {
        printf("No se pudo leer el archivo de datos. Fin del programa.\n");
        return 1;
    }

    // Se muestran los puntos leidos para verificar
    mostrarPuntosDatos(X, Y, n);

    // --- ELEGIR METODO DE INTERPOLACION ---
    int opcion;

    do {
        printf("\n--- MENU DE INTERPOLACION ---\n");
        printf("1. Lagrange (valor interpolado)\n");
        printf("2. Polinomio interpolante (con Gauss y pivoteo)\n");
        printf("3. Lagrange expandido (coeficientes del polinomio)\n");
        printf("0. Salir\n");
        printf("OPCION: ");
        scanf("%d", &opcion);

        switch (opcion) {

            case 1: {
                // ============================================================
                // OPCION 1: LAGRANGE NUMERICO
                // Se pide xsombrero y se calcula Pn(xsombrero)
                // Si xsombrero esta dentro del rango de datos, es interpolacion
                // Si esta fuera, es extrapolacion
                // ============================================================
                double xsombrero;

                printf("\nIngrese el valor x (xsombrero) a interpolar: ");
                scanf("%lf", &xsombrero);

                // Se calcula Pn(xsombrero) con la formula de Lagrange
                double resultado = lagrangeValor(xsombrero, X, Y, n);

                printf("\nEl valor interpolado en x = %lf es: %lf\n", xsombrero, resultado);
                break;
            }

            case 2: {
                // ============================================================
                // OPCION 2: POLINOMIO INTERPOLANTE (VANDERMONDE + GAUSS)
                // Se arma la matriz de Vandermonde:
                //   A[i][j] = x[i]^j   (potencia j de la abscisa i)
                //   b[i]    = y[i]
                // Luego se resuelve A * a = b con eliminacion gaussiana
                // y pivoteo. El vector 'a' son los coeficientes del polinomio.
                // ============================================================
                double A[MAX_TAM][MAX_TAM] = {0};  // matriz de Vandermonde
                double b[MAX_TAM] = {0};           // vector de terminos independientes
                double solucion[MAX_TAM] = {0};    // vector de soluciones (coeficientes)
                double xsombrero, resultado = 0.0;

                // SE ARMA LA MATRIZ DE VANDERMONDE Y EL VECTOR b
                // A[i][j] = x[i]^j   con j = 0, 1, ..., n-1
                for (int i = 0; i < n; i++) {
                    for (int j = 0; j < n; j++) {
                        A[i][j] = pow(X[i], j);
                    }
                    b[i] = Y[i];
                }

                // SE RESUELVE EL SISTEMA CON ELIMINACION GAUSSIANA Y PIVOTEO
                // (adentro se hace el pivoteo y la triangulacion)
                gaussPivoteo(n, A, b, solucion);

                // SE MUESTRAN LOS COEFICIENTES OBTENIDOS
                printf("------------------SOLUCION------------------\n");
                printf("Los coeficientes del polinomio son:\n");
                for (int i = 0; i < n; i++) {
                    printf("a[%d] = %lf\n", i, solucion[i]);
                }

                // SE MUESTRA EL POLINOMIO INTERPOLANTE
                printf("\n------------------POLINOMIO INTERPOLANTE------------------\n");
                mostrarPolinomio(solucion, n);

                // SE EVALUA EL POLINOMIO EN XSOMBRERO (si se quiere)
                printf("Ingrese el valor x (xsombrero) a evaluar: ");
                scanf("%lf", &xsombrero);

                // SE EVALUA: resultado = sum a[i] * xsombrero^i
                for (int i = 0; i < n; i++) {
                    resultado = resultado + solucion[i] * pow(xsombrero, i);
                }

                printf("El valor interpolado en x = %lf es: %lf\n", xsombrero, resultado);
                break;
            }

            case 3: {
                // ============================================================
                // OPCION 3: LAGRANGE EXPANDIDO
                // Calcula los coeficientes del polinomio de Lagrange
                // en su forma clasica y lo muestra
                // ============================================================
                double xsombrero;
                double coef[MAX_PUNTOS];  // coeficientes del polinomio expandido

                printf("\nIngrese el valor x (xsombrero) a interpolar: ");
                scanf("%lf", &xsombrero);

                // Se calcula el valor interpolado (para mostrarlo al final)
                double pn = lagrangeValor(xsombrero, X, Y, n);

                // SE CALCULAN LOS COEFICIENTES DEL POLINOMIO EXPANDIDO
                lagrangeCoeficientes(X, Y, n, coef);

                // SE MUESTRA EL POLINOMIO DE LAGRANGE EN FORMA EXPANDIDA
                printf("------------------POLINOMIO DE LAGRANGE------------------\n");
                mostrarLagrangeExpandido(coef, n);

                printf("\nEl valor interpolado en x = %.2f es: %.6f\n", xsombrero, pn);
                break;
            }

            case 0:
                printf("Saliendo...\n");
                break;

            default:
                printf("Opcion invalida.\n");
        }

    } while (opcion != 0);

    return 0;
}