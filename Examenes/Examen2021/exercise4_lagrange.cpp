#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_POINTS 20
#define MAX_SIZE 100 

/**
 * Function to read Xi, Yi data pairs from a file
 * Expected format: first line contains number of points, then xi yi pairs
 * @param filename Name of the file to read
 * @param X Array to store X values
 * @param Y Array to store Y values
 * @param n Pointer to store the number of data points read
 * @return 1 if successful, 0 otherwise
 */
int read_data_points(const char* filename, double X[], double Y[], int* n);

/**
 * Function to display the data points
 * @param X Array of X values
 * @param Y Array of Y values
 * @param n Number of data points
 */
void print_data_points(double X[], double Y[], int n);

/**
 * Function to print the interpolating polynomial Pn(x)
 * @param a Array of coefficients
 * @param n Degree of polynomial (n-1 is the highest power)
 */
void print_polynomial(double a[], int n);

/**
 * Function to calculate the absolute error between actual and estimated values
 * @param fx The actual function value valuated in X̂
 * @param Pn The estimated polynomial value
 * @return The absolute error
 */
double calculate_error(double fx, double Pn);

/**
 * Function to evaluate the actual function f(x)
 * @param x The point at which to evaluate the function (in our case, X̂)
 * @return The value of the function at x
 */
double valuate_function(double x);

/**
 * Function to compute the Lagrange polynomial at a given point
 * @param x The point at which to evaluate the polynomial
 * @param X Array of X values
 * @param Y Array of Y values
 * @param n Number of data points
 * @return The value of the Lagrange polynomial at x
 */
double lagrange_polinomial(double x, double X[], double Y[], int n);

/**
 * Function to compute the coefficients of the Lagrange polynomial
 * @param X Array of X values
 * @param Y Array of Y values
 * @param n Number of data points
 * @param coeffs Array to store the computed coefficients
 */
void lagrange_coefficients(double X[], double Y[], int n, double coeffs[]);

/**
 * Function to print the symbolic form of the Lagrange polynomial
 * @param coeffs Array of coefficients
 * @param n Number of data points
 */
void print_lagrange_expanded(double coeffs[], int n);

/**
 * Function to perform Gaussian elimination with partial pivoting
 * @param n Size of the system
 * @param a Coefficient matrix
 * @param b Independent terms vector
 * @param X Solution vector
 *  */ 
void gauss_elimination(int n, double a[MAX_SIZE+1][MAX_SIZE+1], double b[MAX_SIZE+1], double X[MAX_SIZE+1]);

/* Lagrange Interpolation:
    X̂ = 2.5
    P3(X̂) = 2.100816
    error = 1.199184
    Pn(x) = 3.96 + 4.510476*x - 2.455147*x² - 2.401723*x³
 */

int main(int argc, char const *argv[]) {
    double X_hat, sum, product, error, fx, Pn;
    double X[MAX_POINTS], Y[MAX_POINTS];
    // Arrays for polynomial coefficients calculation
    double A[MAX_SIZE+1][MAX_SIZE+1], b[MAX_SIZE+1], solution[MAX_SIZE+1];
    // Solution of interpolating polynomial
    double a[MAX_SIZE+1];
    int n, option;

    // Read data points from file
    if (!read_data_points("data4.txt", X, Y, &n)) {
        printf("Failed to read data from file. Exiting.\n");
        return 1;
    }
    
    // Print the data points
    print_data_points(X, Y, n);

    printf("Choose an option:\n");
    printf("1. Lagrange Interpolation\n");
    printf("2. Interpolating Polynomial (using Gaussian elimination)\n");
    printf("3. Polynomial Lagrange Printed\n");
    printf("0. Exit\n");
    scanf("%d", &option);
    
    switch(option) {
        case 1:
            /* If X_hat is between the range of the data points, it's an interpolation
            If X_hat is outside the range, it's an extrapolation */
            // Lagrange Interpolation
            printf("Introduce your X̂ that is to be interpolated: ");
            scanf("%lf", &X_hat);

            // Calculate Pn(X̂)
            sum = 0.0;
            for(int k = 0; k < n; k++) {        
                product = 1.0;
                for(int i = 0; i < n; i++) {    
                    if(i != k) {
                        // Cnk(X̂) = (X̂ - Xi) / (Xk - Xi)
                        product = product * ((X_hat - X[i]) / (X[k] - X[i])); 
                    }
                }
                // Print coefficient Cnk(X̂)
                // We can see that sumCnk = 1
                printf("C%d%d(%.3f) = %.6f\n", n-1, k, X_hat, product);
                
                // Pn(X̂) = Σ Yk * Cnk(X̂)
                sum = sum + (Y[k] * product);
            }

            // error = |f(X̂) - Pn(X̂)|
            fx = 0.0; // Here you can define the actual function f(X̂) if known
            error = 0.0; // If fx is known, calculate error

            fx = valuate_function(X_hat);
            error = calculate_error(fx, sum);

            printf("\nThe interpolated value at X̂ = %lf is: %lf\n", X_hat, sum);
            printf("\nerror = |f(X̂) - Pn(X̂)| = %lf\n", error);
            break;
        case 2:
            // Interpolating Polynomial
            for(int i = 0; i < n; i++) {
                for(int j = 0; j < n; j++) {
                    A[i][j] = pow(X[i], j);
                }
                b[i] = Y[i];
            }
        
            // We use the function from gauss.h to solve the system with Gaussian elimination
            gauss_elimination(n, A, b, solution);
        
            // We copy the solution to a[i] to give relevance to our context
            for(int i = 0; i < n; i++) {
                a[i] = solution[i];
            }
        
            printf("------------------SOLUTION------------------\n");
            printf("The solution of the system is:\n");
            for(int i = 0; i < n; i++) {
                printf("a[%d] = %lf\n", i, a[i]);
            }
            
            printf("\n------------------INTERPOLATING POLYNOMIAL------------------\n");
            // Pn(x) = a0 + a1*x + a2*x^2 + ... + a(n-1)*x^(n-1)
            print_polynomial(a, n);

            // f(x) is defined in the function valuate_function
            fx = 0.0;

            // error = |f(X̂) - Pn(X̂)|
            error = 0.0; // If fx is known, calculate error

            printf("Introduce your X̂ that is to be interpolated: ");
            scanf("%lf", &X_hat);
            sum = 0.0;
            for(int i = 0; i < n; i++) {
                sum = sum + a[i] * pow(X_hat, i);
            }
            fx = valuate_function(X_hat);
            error = calculate_error(fx, sum);

            printf("The interpolated value at X̂ = %lf is: %lf\n", X_hat, sum);

            printf("\nerror = |f(X̂) - Pn(X̂)| = %lf\n", error);
            break;
        case 3:
            /* If X_hat is between the range of the data points, it's an interpolation
            If X_hat is outside the range, it's an extrapolation */
            // Lagrange Interpolation
            printf("Introduce your X̂ that is to be interpolated: ");
            scanf("%lf", &X_hat);

            Pn = lagrange_polinomial(X_hat, X, Y, n);

            // Calculate coefficients of the expanded polynomial
            double coeffs[MAX_POINTS];
            lagrange_coefficients(X, Y, n, coeffs);
            
            printf("------------------LAGRANGE POLYNOMIAL------------------\n");
            print_lagrange_expanded(coeffs, n);

            // error = |f(T) - Pn(T)| 
            // Nota: Para datos experimentales no hay función teórica exacta
            fx = valuate_function(X_hat);
            error = calculate_error(fx, Pn);

            printf("\nThe interpolated value at x_hat = %.2f is: %.6f \n", X_hat, Pn);
            if (fx != 0.0) {
                printf("\nerror = |f(x) - Pn(x)| = %lf\n", error);
            } else {
                printf("\nNote: Error is not calculated - experimental data is used without an exact theoretical function\n");
            }
            break;
        default: 
            printf("Exiting program.\n");
            return 0;
    }

    return 0;
}

int read_data_points(const char* filename, double X[], double Y[], int* n) {
    FILE *fp;
    
    fp = fopen(filename, "r");
    if (fp == NULL) {
        printf("Error: Cannot open file '%s'\n", filename);
        return 0;
    }
    
    printf("File '%s' opened successfully\n", filename);
    
    // Read number of data points
    if (fscanf(fp, "%d", n) != 1) {
        printf("Error: Cannot read number of data points\n");
        fclose(fp);
        return 0;
    }
    
    if (*n <= 0 || *n > MAX_POINTS) {
        printf("Error: Invalid number of points (%d)\n", *n);
        fclose(fp);
        return 0;
    }
    
    // Read data points
    for (int i = 0; i < *n; i++) {
        if (fscanf(fp, "%lf %lf", &X[i], &Y[i]) != 2) {
            printf("Error: Cannot read data point %d\n", i + 1);
            fclose(fp);
            return 0;
        }
    }
    
    fclose(fp);
    printf("Successfully read %d data points\n\n", *n);
    return 1;
}


void print_data_points(double X[], double Y[], int n) {
    printf("Data Points:\n");
    printf("=============\n");
    printf("   i  |      Xi      |      Yi      \n");
    printf("------|--------------|-------------\n");
    for (int i = 0; i < n; i++) {
        printf("%4d  | %12.6f | %12.6f\n", i + 1, X[i], Y[i]);
    }
    printf("\n");
}

void print_polynomial(double a[], int n) {
    printf("Pn(x) = ");
    
    // Handle the first term (constant term a0)
    if (fabs(a[0]) > 1e-10) {  // Avoid printing very small values as zero
        printf("%.6f", a[0]);
    } else {
        printf("0");
    }
    
    // Handle the rest of the terms
    for (int i = 1; i < n; i++) {
        if (fabs(a[i]) > 1e-10) {  // Only print if coefficient is significant
            // Print sign
            if (a[i] > 0) {
                printf(" + ");
            } else {
                printf(" - ");
            }
            
            // Print coefficient (absolute value since sign is already printed)
            double coeff = fabs(a[i]);
            if (coeff != 1.0) {
                printf("%.6f", coeff);
            }
            
            // Print variable part
            if (i == 1) {
                printf("x");
            } else {
                printf("x^%d", i);
            }
        }
    }
    printf("\n\n");
}

double lagrange_polinomial(double x, double X[], double Y[], int n) {
    double sum, product;

    // Calculate Pn(X̂)
    sum = 0.0;
    for(int k = 0; k < n; k++) {        
        product = 1.0;
        for(int i = 0; i < n; i++) {    
            if(i != k) {
                // Cnk(X̂) = (X̂ - Xi) / (Xk - Xi)
                product = product * ((x - X[i]) / (X[k] - X[i])); 
            }
        }
        // Print coefficient Cnk(X̂)
        // We can see that sumCnk = 1
        printf("C%d%d(%.3f) = %.6f\n", n-1, k, x, product);
        
        // Pn(X̂) = Σ Yk * Cnk(X̂)
        sum = sum + (Y[k] * product);
    }

    return sum;
}



void lagrange_coefficients(double X[], double Y[], int n, double coeffs[]) {
    for (int i = 0; i < n; i++)
        coeffs[i] = 0.0;

    for (int k = 0; k < n; k++) {
        double base[MAX_POINTS] = {1.0}; 
        int degree = 0;

        double denom = 1.0;
        for (int i = 0; i < n; i++) {
            if (i == k) continue;

            // Multiplies (x - X[i]) to the base polynomial
            for (int j = degree; j >= 0; j--)
                base[j+1] += base[j];
            degree++;

            for (int j = degree; j >= 0; j--)
                base[j] = (j>0 ? base[j-1] : 0) - X[i]*(j>=0?base[j]:0);

            denom *= (X[k] - X[i]);
        }

        double factor = Y[k] / denom;
        for (int j = 0; j <= degree; j++)
            coeffs[j] += base[j] * factor;
    }
}

void print_lagrange_expanded(double coeffs[], int n) {
    printf("P_%d(x) = ", n-1);
    int first = 1;
    for (int i = 0; i < n; i++) {
        if (fabs(coeffs[i]) < 1e-12) continue; 
        if (!first) {
            if (coeffs[i] >= 0) printf(" + ");
            else printf(" - ");
        } else if (coeffs[i] < 0) {
            printf("-");
        }
        printf("%.6f", fabs(coeffs[i]));
        if (i > 0) printf("*x^%d", i);
        first = 0;
    }
    printf("\n");
}

double calculate_error(double fx, double Pn) {
    return fabs(fx - Pn);
}

double valuate_function(double x) {
    return x + (2/x);
}

void gauss_elimination(int n, double a[MAX_SIZE+1][MAX_SIZE+1], double b[MAX_SIZE+1], double X[MAX_SIZE+1]) {
    int p;
    double factor, product, sum, aux;

    // Walk the rows of the matrix (Gaussian elimination)
    for(int i = 0; i <= n-2; i++) {
        // Partial pivoting
        p = i;
        if(fabs(a[i][i]) < 1e-5) {
            for(int l = i+1; l <= n-1; l++) {
                if(fabs(a[l][i]) > fabs(a[p][i])) {
                    p = l;
                }
            }
            for(int m = i; m <= n-1; m++) {
                aux = a[p][m];
                a[p][m] = a[i][m];
                a[i][m] = aux;
            }
            aux = b[p];
            b[p] = b[i];
            b[i] = aux;
        }

        // Make zero the elements below the diagonal
        for(int j = i+1; j <= n-1; j++) {
            factor = a[j][i] / a[i][i];
            for(int k = i; k <= n-1; k++) {
                a[j][k] -= factor * a[i][k];
            }
            b[j] -= factor * b[i];
        }
    }

    // We print the matrix after Gaussian elimination, It's more convenient
    printf("The matrix after Gaussian elimination is:\n");
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            printf("%lf ", a[i][j]);
        }
        printf("\n");
    }
    printf("The vector b after Gaussian elimination is:\n");
    for(int i = 0; i < n; i++) {
        printf("%lf ", b[i]);
    }
    printf("\n\n");

    // Determinant check
    product = 1.0;
    for(int i = 0; i < n; i++) {
        product *= a[i][i];
    }

    printf("------------------DETERMINANT------------------\n");
    printf("The determinant of the matrix is: %lf\n\n", product);

    if(fabs(product) < 1e-10) {
        printf("Error: determinant is zero, no unique solution.\n");
        exit(1);
    }

    // Back substitution
    X[n-1] = b[n-1] / a[n-1][n-1];
    for(int i = n-2; i >= 0; i--) {
        sum = b[i];
        for(int j = i+1; j < n; j++) {
            sum -= a[i][j] * X[j];
        }
        X[i] = sum / a[i][i];
    }
}