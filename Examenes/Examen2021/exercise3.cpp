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
 *  */ 
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

/*  w = 1.8, tolerance = 1e-7
    Iterations number = 2373
    error = 0.0000001
    X1 = 1
    X2 = -194.9999999
    X50 = -4898.9999985
    X99 = -195.0000000
    X100 = 1.0000000
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
    if(!read_array_file("data.dat", a, b, &n)) {
        return 1;
    }
    
    printf("Original system of equations:\n");
    printf("===============================\n");
    for(int i = 1; i <= n; i++) {
        for(int j = 1; j <= n; j++) {
            printf("%10.7lf ", a[i][j]);  // 10 total spaces, 6 decimal places
        }
        printf("| %10.7lf\n", b[i]);       // 10 total spaces, 6 decimal places
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
        printf("Xn[%d] = %10.7lf\n", i, Xn[i]);
    }
    printf("The method converged in %d iterations with an error of %10.7lf\n", iterations, error);
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
