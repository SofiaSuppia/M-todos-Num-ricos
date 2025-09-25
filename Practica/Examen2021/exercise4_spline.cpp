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
 * Function to print the cubic spline coefficients
 * @param X Array of X values (data points)
 * @param solution Array containing the solution coefficients
 * @param n Number of data points
 */
void print_cubic_splines(double X[], double solution[], int n);

/**
 * Function to evaluate the cubic spline at a given point x
 * @param X Array of X values (data points)
 * @param solution Array of spline coefficients
 * @param n Number of data points
 * @param x The x value to evaluate
 * @return The interpolated y value
 */
double evaluate_spline(double X[], double solution[], int n, double x);

/**
 * Function to compute and print the linear splines
 * @param X Array of X values (data points)
 * @param Y Array of Y values (data points)
 * @param n Number of data points
 */
void linear_spline(double X[], double Y[], int n);

/* X_hat = 2.5
  Spline 2 (from X[1]=2.100 to X[2]=2.800):
  S2(x) = -1.453470x³ + 10.251750x² + -22.025026x + 15.802920
  S2(2.5) = 2.103
*/

int main(int argc, char const *argv[]) {
    // Arrays for data points to read from text file
    double X[MAX_POINTS], Y[MAX_POINTS];
    // Arrays for polynomial coefficients calculation
    double A[MAX_SIZE+1][MAX_SIZE+1], b[MAX_SIZE+1], solution[MAX_SIZE+1];
    // Solution of Least Squares polynomial
    double a[MAX_SIZE+1];
    // n = number of data points and degree is the polynomial degree
    int n, option, row1, row2;
    

    // Read data points from file
    if (!read_data_points("data4.txt", X, Y, &n)) {
        printf("Failed to read data from file. Exiting.\n");
        return 1;
    }
    
    // Print the data points
    print_data_points(X, Y, n);

    // Initialize A and b to zero
    for (int i = 0; i < 4*(n-1); i++) {
        b[i] = 0.0;
        for (int j = 0; j < 4*(n-1); j++) {
            A[i][j] = 0.0;
        }
    }

    printf("Choose an option for spline interpolation:\n");
    printf("1. Linear Spline\n");
    printf("2. Cubic Spline\n");
    scanf("%d", &option);

    switch(option) {
        case 1:
            linear_spline(X, Y, n);
            // We evalute points x_hat to verify if the spline works correctly
            printf("Do you want to evaluate a point x_hat? (1 for yes, 0 for no): ");
            scanf("%d", &option);
            while(option) {
                double x_hat;
                printf("Enter the value of x_hat: ");
                scanf("%lf", &x_hat);
                // Find the right interval for x_hat
                if(x_hat < X[0] || x_hat > X[n-1]) {
                    printf("x_hat is out of bounds [%lf, %lf]. Please enter a value within the range.\n", X[0], X[n-1]);
                } else {
                    // Locate the interval
                    int k = 0;
                    while(k < n-1 && x_hat > X[k+1]) {
                        k++;
                    }
                    // Linear interpolation formula
                    double y_hat = Y[k] + (Y[k+1] - Y[k]) * (x_hat - X[k]) / (X[k+1] - X[k]);
                    printf("The interpolated value at x_hat = %.3f is y_hat = %.3f\n", x_hat, y_hat);
                }
                printf("Do you want to evaluate another point x_hat? (1 for yes, 0 for no): ");
                scanf("%d", &option);
            }
        break;
        case 2:
            // Initialize A and b
            for (int i = 0; i < 4*(n-1); i++) {
                b[i] = 0.0;
                for (int j = 0; j < 4*(n-1); j++) {
                    A[i][j] = 0.0;
                }
            }

            // We calculate A[4(n-1)][4(n-1)] and b[4(n-1)]
            // First 2(n-1) equations: Each spline passes through its two endpoints
            for(int k = 0; k < n-1; k++) {
                // Spline k passes through point (X[k], Y[k])
                for(int j = 0; j <= 3; j++) {
                    A[2*k][4*k+j] = pow(X[k], 3-j);
                }
                b[2*k] = Y[k];

                // Spline k passes through point (X[k+1], Y[k+1])
                for(int j = 0; j <= 3; j++) {
                    A[2*k+1][4*k+j] = pow(X[k+1], 3-j);
                }
                b[2*k+1] = Y[k+1];
            }
        
            // Next n-2 equations: Continuity of first derivatives
            for(int k = 0; k < n-2; k++) {
                int row = 2*(n-1) + k;
                // First derivative of spline k at X[k+1]
                for(int j = 0; j <= 2; j++) {
                    A[row][4*k+j] = (3-j) * pow(X[k+1], 2-j);
                }
                // First derivative of spline k+1 at X[k+1]
                for(int j = 0; j <= 2; j++) {
                    A[row][4*(k+1)+j] = -(3-j) * pow(X[k+1], 2-j);
                }
                b[row] = 0.0;
            }
        
            // Next n-2 equations: Continuity of second derivatives
            for(int k = 0; k < n-2; k++) {
                int row = 2*(n-1) + (n-2) + k;
                // Second derivative of spline k at X[k+1]
                A[row][4*k] = 6 * X[k+1];
                A[row][4*k+1] = 2;
                // Second derivative of spline k+1 at X[k+1]
                A[row][4*(k+1)] = -6 * X[k+1];
                A[row][4*(k+1)+1] = -2;
                b[row] = 0.0;
            }
        
            // Two boundary conditions: Natural spline (second derivatives = 0 at endpoints)
            row1 = 4*(n-1) - 2;
            row2 = 4*(n-1) - 1;

            // Second derivative = 0 at X[0] (first spline)
            A[row1][0] = 6 * X[0];
            A[row1][1] = 2;
            b[row1] = 0.0;

            // Second derivative = 0 at X[n-1] (last spline)
            A[row2][4*(n-2)] = 6 * X[n-1];
            A[row2][4*(n-2)+1] = 2;
            b[row2] = 0.0;
        
            // We use the function from gauss.h to solve the system with Gaussian elimination
            gauss_elimination(4*(n-1), A, b, solution);
        
            // Print the cubic spline coefficients
            print_cubic_splines(X, solution, n);
        
            // We evalute points x_hat to verify if the spline works correctly
            printf("Do you want to evaluate a point x_hat? (1 for yes, 0 for no): ");
            scanf("%d", &option);
            while(option) {
                double x_hat;
                printf("Enter the value of x_hat: ");
                scanf("%lf", &x_hat);
                double y_hat = evaluate_spline(X, solution, n, x_hat);
                printf("The interpolated value at x_hat = %.3f is y_hat = %.3f\n", x_hat, y_hat);
                printf("Do you want to evaluate another point x_hat? (1 for yes, 0 for no): ");
                scanf("%d", &option);
            }
            break;
        default:
            printf("Invalid option. Please try again.\n");
            break;
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

void print_cubic_splines(double X[], double solution[], int n) {
    printf("------------------SOLUTION------------------\n");
    printf("Cubic Spline Coefficients:\n");
    for(int k = 0; k < n-1; k++) {
        printf("Spline %d (from X[%d]=%.3f to X[%d]=%.3f):\n", k+1, k, X[k], k+1, X[k+1]);
        printf("  S%d(x) = %.6fx³ + %.6fx² + %.6fx + %.6f\n", 
               k+1, solution[4*k], solution[4*k+1], solution[4*k+2], solution[4*k+3]);
    }
    printf("\n");
}

/**
 * Function to evaluate the cubic spline at a given point x
 * @param X Array of X values (data points)
 * @param solution Array of spline coefficients
 * @param n Number of data points
 * @param x The x value to evaluate
 * @return The interpolated y value
 */
double evaluate_spline(double X[], double solution[], int n, double x) {
    int k;

    // Find the correct spline interval
    if (x <= X[0]) {
        k = 0;  // Use first spline for x smaller than first point
    } else if (x >= X[n-1]) {
        k = n - 2; // Use last spline for x larger than last point
    } else {
        for (k = 0; k < n-1; k++) {
            if (x >= X[k] && x <= X[k+1]) {
                break;
            }
        }
    }

    // Evaluate S_k(x) = a*x^3 + b*x^2 + c*x + d
    double y = solution[4*k] * pow(x, 3) +
               solution[4*k+1] * pow(x, 2) +
               solution[4*k+2] * x +
               solution[4*k+3];

    return y;
}

/**
 * Function to compute and print the linear splines
 * @param X Array of X values (data points)
 * @param Y Array of Y values (data points)
 * @param n Number of data points
 */
void linear_spline(double X[], double Y[], int n) {
    printf("------------------LINEAR SPLINES------------------\n");
    for (int k = 0; k < n - 1; k++) {
        double mk = (Y[k+1] - Y[k]) / (X[k+1] - X[k]);
        printf("Spline %d (from X[%d]=%.3f to X[%d]=%.3f):\n", k+1, k, X[k], k+1, X[k+1]);
        printf("  f_%d(X) = %.6f + %.6f * (X - %.6f)\n\n", k+1, Y[k], mk, X[k]);
    }
}
