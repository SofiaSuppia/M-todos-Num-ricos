#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAX_POINTS 50
#define MAX_SIZE 100 
#include "gauss.h"

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

/* We can see that the efficiently extrapolation occurs only in certain functions because there are some cases where
the polinomial fitting curve don't fit well outside the boundaries or intervals of the measurements
You can see it in the graphs of the three exercises (2.1 2.2 and 2.3) */

int main(int argc, char const *argv[]) {
    double X_hat, sum, product, error, fx;
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

double calculate_error(double fx, double Pn) {
    return fabs(fx - Pn);
}

double valuate_function(double x) {
    return pow(x, x);
}