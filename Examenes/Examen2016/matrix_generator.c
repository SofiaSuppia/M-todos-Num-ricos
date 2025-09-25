#include <stdio.h>
#include <stdlib.h>

#define NOMBRE_ARCHIVO_MATRIZ "matriz.txt"
#define N_TAMANO 100 // Tamaño de la matriz N x N

int main(int argc, char const *argv[])
{
    // Abrir el archivo para escritura
    FILE *archivo = fopen(NOMBRE_ARCHIVO_MATRIZ, "w");
    if (archivo == NULL) {
        printf("No se pudo abrir el archivo para escritura.\n");
        return 1;
    }

    // Asignación de memoria correcta para una matriz 2D (punteros a punteros)
    double **A = (double **)malloc(N_TAMANO * sizeof(double *));
    // Asignación de memoria para el vector de términos independientes
    double *b = (double *)malloc(N_TAMANO * sizeof(double));
    if (!A || !b)
    {
        printf("Error de asignación de memoria\n");
        fclose(archivo);
        return 1;
    }

    // Inicializar la matriz A (asignar memoria para cada fila)
    for (size_t i = 0; i < N_TAMANO; i++)
    {
        A[i] = (double *)malloc(N_TAMANO * sizeof(double));
        if (!A[i]) {
            printf("Error de memoria en la fila %zu\n", i);
            // Aquí se necesitaría liberar la memoria ya asignada antes de salir.
            fclose(archivo);
            return 1;
        }
        for (size_t j = 0; j < N_TAMANO; j++)
        {
            A[i][j] = 0;
        }
    }

    // Construir la matriz tridiagonal A y el vector b
    for (size_t i = 0; i < N_TAMANO; i++)
    {
        // Diagonal Principal: A[i][i] = 2
        A[i][i] = 2;

        // Vector Independiente: b[i] = 6 (valor por defecto)
        b[i] = 6;

        // Subdiagonal: A[i][i-1] = 1
        if (i > 0)
            A[i][i-1] = 1;

        // Superdiagonal: A[i][i+1] = 1
        if (i < N_TAMANO-1)
            A[i][i+1] = 1;
    }

    // Modificaciones en los bordes del vector independiente
    b[0] = 4.5;
    b[N_TAMANO-1] = 4.5;

    // Guardar la matriz y el vector en el archivo
    printf("Guardando matriz de %dx%d en %s...\n", N_TAMANO, N_TAMANO, NOMBRE_ARCHIVO_MATRIZ);
    for (size_t i = 0; i < N_TAMANO; i++)
    {
        for (size_t j = 0; j < N_TAMANO; j++)
        {
            // Escribe el elemento de la matriz con un decimal y un espacio
            fprintf(archivo, "%.1lf ", A[i][j]);
        }
        // Escribe el término independiente y un salto de línea
        // Cuando quiera leer el programa, recuerde quitar el | y el tabulador (si es que los usa).
        fprintf(archivo, "%.1lf\n", b[i]);
    }
    printf("Matriz guardada exitosamente.\n");

    // Liberar memoria
    fclose(archivo);
    // Liberar las filas
    for (size_t i = 0; i < N_TAMANO; i++) {
        free(A[i]);
    }
    // Liberar el puntero principal de la matriz
    free(A);
    // Liberar el vector independiente
    free(b);

    return 0;
}