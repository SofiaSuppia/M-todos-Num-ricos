#include <cstdio>
#include <stdlib.h>
#include <cmath>

/* This file have bug fixes done with Github Copilot in Gauss-Siedel and Relaxation Method. */

// Now you can handle up to 50x50 matrixs
#define MAX_SIZE 100  

/**
 * Reads an augmented matrix from a file text called data.dat
 * The expected format is:
 *   a11 a12 ... a1n b1
 *   a21 a22 ... a2n b2
 *   ...
 *   an1 an2 ... ann bn
 * @param filename Name of the file to read
 * @param a Matrix of coefficients
 * @param b Vector of independent terms
 * @param n Pointer to the size of the system (number of equations)
 * @return true if the file was read successfully, false otherwise
 *  */
bool read_array_file(const char* filename, double a[][MAX_SIZE+1], double b[], int* n);

/**
 * Checks if the matrix is ​​diagonally dominant and that there are no zeros on the diagonal.
 * Returns 0 if it's everything ok, otherwise returns 1
 * @param a Matrix of coefficients
 * @param n Size of the matrix
 * @return 0 if everything is ok, 1 if there is a zero on the diagonal
 *         or a warning if the matrix is not diagonally dominant
 */
int checkDiagonalDominance(double a[][MAX_SIZE+1], int n);

/** Function to initialize the guess vector with zeros
 * @param Xv Vector to initialize == Old X
 * @param n Size of the vector
 * 
 **/ 
void initializeGuess(double Xv[], int n);

/** Function to compute the error between Xn and Xv
 * @param Xn New X = Solution of the matrix
 * @param Xv Old X = Previous iteration of the solution
 * @param n Size of the vectors
 * @return The computed error (Euclidean norm)
 *  */
double computeError(double Xn[], double Xv[], int n);


/** Function to print the solution
 * @param methodName Name of the method used
 * @param Xn Solution vector
 * @param n Size of the vector
 * @param iterations Number of iterations taken to converge
 * @param error Final error
 *  */ 
void printSolution(const char* methodName, double Xn[], int n, int iterations, double error);


/**
 * Implementation of the Jacobi method for solving linear systems
 * @param a Matrix of coefficients
 * @param b Vector of independent terms
 * @param Xv Old X = Previous iteration of the solution
 * @param Xn New X = Solution of the matrix
 * @param n Size of the vectors
 *  */
void jacobiMethod(double a[][MAX_SIZE+1], double b[], double Xv[], double Xn[], int n);

/**
 * Implementation of the Jacobi method for solving linear systems
 * @param a Matrix of coefficients
 * @param b Vector of independent terms
 * @param Xv Old X = Previous iteration of the solution
 * @param Xn New X = Solution of the matrix
 * @param n Size of the vectors
 *  */
void gaussSeidel(double a[][MAX_SIZE+1], double b[], double Xv[], double Xn[], int n);


/**
 * Implementation of the Relaxation method for solving linear systems
 * It's similar to Gauss-Seidel but with a relaxation factor omega and one additional line of code
 * @param a Matrix of coefficients
 * @param b Vector of independent terms
 * @param Xv Old X = Previous iteration of the solution
 * @param Xn New X = Solution of the matrix
 * @param n Size of the vectors
 *  */
void relaxationMethod(double a[][MAX_SIZE+1], double b[], double Xv[], double Xn[], int n);

/* Conviene utilizar el metodo de relajacion ya que es preferible utilizarlo para ecuaciones muy grandes y que requieren de una
convergencia rapida, esto se logra utilizando un factor de relajacion de valor entre 1 y 2.
Gauss no conviene ya que es para sistemas pequeños a medianos, la ventaja es que te devuelve los valores de la solucion exacta
Jacobi conviene si el sistema es grande y la convergencia esta asegurada.
Gauss-Seidel conviene si Jacobi converge muy lentamente 
Metodo Relajacion usado:
    tolerance = 1e-4 = error deseado
    omega = 1.6
    iterations number = 117
    error = 0.0001 = error obtenido
    The solution of the system is:
    Xn[1] =   1.500073
    Xn[2] =   1.499855
    Xn[3] =   1.500215
    Xn[4] =   1.499716
    Xn[5] =   1.500352
    Xn[6] =   1.499582
    Xn[7] =   1.500482
    Xn[8] =   1.499456
    Xn[9] =   1.500604
    Xn[10] =   1.499338
    Xn[11] =   1.500717
    Xn[12] =   1.499229
    Xn[13] =   1.500822
    Xn[14] =   1.499129
    Xn[15] =   1.500918
    Xn[16] =   1.499037
    Xn[17] =   1.501006
    Xn[18] =   1.498953
    Xn[19] =   1.501086
    Xn[20] =   1.498876
    Xn[21] =   1.501160
    Xn[22] =   1.498806
    Xn[23] =   1.501227
    Xn[24] =   1.498741
    Xn[25] =   1.501290
    Xn[26] =   1.498680
    Xn[27] =   1.501349
    Xn[28] =   1.498622
    Xn[29] =   1.501407
    Xn[30] =   1.498565
    Xn[31] =   1.501464
    Xn[32] =   1.498507
    Xn[33] =   1.501522
    Xn[34] =   1.498449
    Xn[35] =   1.501581
    Xn[36] =   1.498389
    Xn[37] =   1.501643
    Xn[38] =   1.498325
    Xn[39] =   1.501708
    Xn[40] =   1.498258
    Xn[41] =   1.501776
    Xn[42] =   1.498188
    Xn[43] =   1.501848
    Xn[44] =   1.498115
    Xn[45] =   1.501922
    Xn[46] =   1.498040
    Xn[47] =   1.501998
    Xn[48] =   1.497963
    Xn[49] =   1.502075
    Xn[50] =   1.497886
    Xn[51] =   1.502151
    Xn[52] =   1.497811
    Xn[53] =   1.502225
    Xn[54] =   1.497740
    Xn[55] =   1.502294
    Xn[56] =   1.497674
    Xn[57] =   1.502357
    Xn[58] =   1.497615
    Xn[59] =   1.502411
    Xn[60] =   1.497567
    Xn[61] =   1.502453
    Xn[62] =   1.497530
    Xn[63] =   1.502483
    Xn[64] =   1.497507
    Xn[65] =   1.502499
    Xn[66] =   1.497500
    Xn[67] =   1.502497
    Xn[68] =   1.497510
    Xn[69] =   1.502478
    Xn[70] =   1.497539
    Xn[71] =   1.502440
    Xn[72] =   1.497586
    Xn[73] =   1.502382
    Xn[74] =   1.497654
    Xn[75] =   1.502305
    Xn[76] =   1.497741
    Xn[77] =   1.502208
    Xn[78] =   1.497848
    Xn[79] =   1.502092
    Xn[80] =   1.497974
    Xn[81] =   1.501957
    Xn[82] =   1.498117
    Xn[83] =   1.501806
    Xn[84] =   1.498276
    Xn[85] =   1.501639
    Xn[86] =   1.498450
    Xn[87] =   1.501458
    Xn[88] =   1.498637
    Xn[89] =   1.501266
    Xn[90] =   1.498834
    Xn[91] =   1.501065
    Xn[92] =   1.499039
    Xn[93] =   1.500856
    Xn[94] =   1.499250
    Xn[95] =   1.500643
    Xn[96] =   1.499464
    Xn[97] =   1.500428
    Xn[98] =   1.499680
    Xn[99] =   1.500213
    Xn[100] =   1.499894
*/

int main(int argc, char const *argv[]) {
    int n, p;
    double factor, product, sum, aux;
    double Xv[MAX_SIZE+1], Xn[MAX_SIZE+1]; // Old X, New X
    double tolerance, old_error, new_error;
    int iterations;
    double omega;

    // Defining arrays using the global MAX_SIZE
    double a[MAX_SIZE+1][MAX_SIZE+1], b[MAX_SIZE+1], X[MAX_SIZE+1];

    // Read array from file using function
    if(!read_array_file("C:\\Users\\sofim\\OneDrive\\Documentos\\M-todos-Num-ricos\\Examenes\\Examen2016\\data.dat", a, b, &n)) {
        return 1;
    }
    
    printf("Original system of equations:\n");
    printf("===============================\n");
    for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= n; j++) {
            printf("%10.6lf ", a[i][j]);  // 10 total spaces, 6 decimal places
        }
        printf("| %10.6lf\n", b[i]);       // 10 total spaces, 6 decimal places
    }
    printf("\n");

    // Verificación de diagonal dominante
    if (checkDiagonalDominance(a, n) != 0) {
        printf("The method cannot continue. The matrix has zeros on the diagonal. The program exits.");
        return 1;
    }
    printf("Verificación completada.\n");

    printf("Choose a method to solve the system:\n");
    printf("1. Jacobi Method\n");
    printf("2. Gauss-Seidel\n");
    printf("3. Relaxation Method\n");
    printf("Enter option: ");

    int option;
    scanf("%d", &option);

    switch(option) {
        case 1:
            jacobiMethod(a, b, Xv, Xn, n);
            break;
        case 2:
            gaussSeidel(a, b, Xv, Xn, n);
            break;
        case 3:
            relaxationMethod(a, b, Xv, Xn, n);
            break;
        default:
            printf("❌ Invalid option.\n");
    }


    return 0;
}

int checkDiagonalDominance(double a[][MAX_SIZE+1], int n) {
    for (int i = 1; i <= n; i++) {
        double sum = 0.0;

        // We first check if there is zero on the diagonal
        if (fabs(a[i][i]) == 0.0) {
            printf("❌ Error: Zero element on the diagonal at position a[%d][%d].\n", i, i);
            return 1;
        }

        // Sum of off-diagonal elements (elementos fuera de la diagonal)
        for (int j = 1; j <= n; j++) {
            if (j != i) {
                sum += fabs(a[i][j]);
            }
        }

        // We check the dominance condition
        if (fabs(a[i][i]) < sum) {
            printf("⚠️  Warning: The matrix is not diagonally dominant at row %d.\n", i);
        }
    }

    return 0; // Everything is OK
}

void initializeGuess(double Xv[], int n) {
    for(int i = 1; i <= n; i++) {
        Xv[i] = 0.0;
    }
}

double computeError(double Xn[], double Xv[], int n) {
    double error = 0.0;
    for(int i = 1; i <= n; i++) {
        error += pow(Xn[i] - Xv[i], 2);
    }
    return sqrt(error);
}

void printSolution(const char* methodName, double Xn[], int n, int iterations, double error) {
    printf("------------------SOLUTION OF %s------------------\n", methodName);
    printf("The solution of the system is:\n");
    for(int i = 1; i <= n; i++) {
        printf("Xn[%d] = %10.6lf\n", i, Xn[i]);
    }
    printf("The method converged in %d iterations with an error of %10.6lf\n", iterations, error);
}


void jacobiMethod(double a[][MAX_SIZE+1], double b[], double Xv[], double Xn[], int n) {
    double sum, tolerance, old_error, new_error;
    int iterations;

    initializeGuess(Xv, n);

    printf("Please enter tolerance:");
    scanf("%lf", &tolerance);

    old_error = 1000;
    iterations = 0;

    do {
        iterations++;
        for(int i = 1; i <= n; i++) {
            sum = 0.0;
            for(int j = 1; j <= n; j++) {
                if(j != i) {
                    sum += a[i][j] * Xv[j];
                }
            }
            Xn[i] = (b[i] - sum) / a[i][i];
        }

        new_error = computeError(Xn, Xv, n);

        if(new_error > old_error) {
            printf("The method does not converge, we stop the process.\n");
            return;
        }

        old_error = new_error;

        for(int i = 1; i <= n; i++) {
            Xv[i] = Xn[i];
        }
    } while(new_error > tolerance);

    printSolution("JACOBI METHOD", Xn, n, iterations, new_error);
}

void gaussSeidel(double a[][MAX_SIZE+1], double b[], double Xv[], double Xn[], int n) {
    double sum, tolerance, old_error, new_error;
    int iterations;

    initializeGuess(Xv, n);

    printf("Please enter tolerance:");
    scanf("%lf", &tolerance);

    old_error = 1000;
    iterations = 0;

    do {
        iterations++;
        for(int i = 1; i <= n; i++) {
            sum = 0.0;  // Reset sum for each row
            
            // Sum elements before diagonal (using NEW values Xn)
            for(int j = 1; j <= i-1; j++) {
                sum += a[i][j] * Xn[j];
            }
            
            // Sum elements after diagonal (using OLD values Xv)
            for(int j = i+1; j <= n; j++) {
                sum += a[i][j] * Xv[j];
            }
            
            Xn[i] = (b[i] - sum) / a[i][i];
        }

        new_error = computeError(Xn, Xv, n);

        if(new_error > old_error) {
            printf("The method does not converge, we stop the process.\n");
            return;
        }

        old_error = new_error;

        for(int i = 1; i <= n; i++) {
            Xv[i] = Xn[i];
        }
    } while(new_error > tolerance);

    printSolution("GAUSS-SEIDEL", Xn, n, iterations, new_error);
}

void relaxationMethod(double a[][MAX_SIZE+1], double b[], double Xv[], double Xn[], int n) {
    double sum, tolerance, old_error, new_error, omega;
    int iterations;

    initializeGuess(Xv, n);

    printf("Please enter tolerance:");
    scanf("%lf", &tolerance);

    printf("Please enter relaxation factor (0 < omega < 2):");
    scanf("%lf", &omega);

    old_error = 1000;
    iterations = 0;

    do {
        iterations++;
        for(int i = 1; i <= n; i++) {
            sum = 0.0;  // Reset sum for each row
            
            // Sum elements before diagonal (using NEW values Xn)
            for(int j = 1; j <= i-1; j++) {
                sum += a[i][j] * Xn[j];
            }
            
            // Sum elements after diagonal (using OLD values Xv)
            for(int j = i+1; j <= n; j++) {
                sum += a[i][j] * Xv[j];
            }
            
            // Calculate Gauss-Seidel step
            double gauss_seidel = (b[i] - sum) / a[i][i];
            
            // Apply relaxation factor (SOR)
            Xn[i] = omega * gauss_seidel + (1.0 - omega) * Xv[i];
        }

        new_error = computeError(Xn, Xv, n);

        if(new_error > old_error) {
            printf("The method does not converge, we stop the process.\n");
            return;
        }

        old_error = new_error;

        for(int i = 1; i <= n; i++) {
            Xv[i] = Xn[i];
        }
    } while(new_error > tolerance);

    printSolution("RELAXATION METHOD", Xn, n, iterations, new_error);
}


bool read_array_file(const char* filename, double a[][MAX_SIZE+1], double b[], int* n) {
    FILE *fp;
    char c;
    
    // Open data file
    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("❌ Error: Cannot open file '%s'\n", filename);
        return false;
    }
    
    printf("✅ File '%s' opened\n\n", filename);

    // Count rows in file
    /* int rows = 0;
    while((c = fgetc(fp)) != EOF) {
        if(c == '\n') {
            rows++;
        }
    } */

    // This alternative works better than the previous one
    int rows = 0;
    while(!feof(fp)) {
        char buffer[1024];
        if(fgets(buffer, sizeof(buffer), fp) != NULL) {
            rows++;
        }
    }
    
    // System must be square, so n = rows
    *n = rows;
    printf("📊 Size of the system: %d x %d\n", *n, *n);

    // Close and reopen the file to reset the pointer
    fclose(fp);
    fp = fopen(filename, "r");
    
    // Check maximum size using global MAX_SIZE
    if(*n > MAX_SIZE) {
        printf("❌ Error: System too big (%d). Maximum allowed: %d\n", *n, MAX_SIZE);
        fclose(fp);
        return false;
    }

    // Read augmented matrix from file
    // Expected format: each row contains n coefficients + 1 independent term
    // Example for 3x3: a11 a12 a13 b1
    //                   a21 a22 a23 b2  
    //                   a31 a32 a33 b3
    
    int i, j;
    for(i = 1; i <= *n; i++) {
        // Reading the matrix coefficients
        for(j = 1; j <= *n; j++) {
            if(fscanf(fp, "%lf", &a[i][j]) != 1) {
                printf("❌ Error reading element a[%d][%d]\n", i, j);
                fclose(fp);
                return false;
            }
        }
        // Read the term independent
        if(fscanf(fp, "%lf", &b[i]) != 1) {
            printf("❌ Error reading independent term b[%d]\n", i);
            fclose(fp);
            return false;
        }
    }
    
    fclose(fp);
    printf("✅ Array read successfully from file\n\n");
    
    return true;
}
