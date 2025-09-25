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
 * Function to evaluate the cubic spline at a given x
 * @param x The x value to evaluate
 * @param X Array of X values (data points)
 * @param solution Array containing the solution coefficients
 * @param n Number of data points
 * @return The evaluated spline value at x
 */
double evaluate_spline(double x, double X[], double solution[], int n);

/* Results for Lineal Scale:
  Cubic Spline Coefficients:
  Spline 1 (from X[0]=0.200 to X[1]=2.000):
    S1(x) = 0.700110x³ + -0.420066x² + -51.684344x + 113.348071
  Spline 2 (from X[1]=2.000 to X[2]=20.000):
    S2(x) = -0.073175x³ + 4.219642x² + -60.963761x + 119.534349
  Spline 3 (from X[2]=20.000 to X[3]=200.000):
    S3(x) = 0.000331x³ + -0.190691x² + 27.242909x + -468.510122
  Spline 4 (from X[3]=200.000 to X[4]=2000.000):
    S4(x) = -0.000002x³ + 0.008859x² + -12.667063x + 2192.154679
  Spline 5 (from X[4]=2000.000 to X[5]=20000.000):
    S5(x) = 0.000000x³ + -0.000401x² + 5.852124x + -10153.969597

  f(5) = S2(5) = -88.94 = CD1
  f(50) = S3(50) = 458.283 = CD2
  f(500) = S4(500) = -2176.6268 = CD3
  f(5000) = S5(5000) = -24980893.35 = CD4

Results for Logaritmic Scale:
  Cubic Spline Coefficients:
  Spline 1 (from X[0]=-0.699 to X[1]=0.301):
    S1(x) = 20.240541x³ + 42.442592x² + -79.674442x + 33.486146
  Spline 2 (from X[1]=0.301 to X[2]=1.301):
    S2(x) = -23.282703x³ + 81.747998x² + -91.506548x + 34.673419
  Spline 3 (from X[2]=1.301 to X[3]=2.301):
    S3(x) = 4.230273x³ + -25.637623x² + 48.205367x + -25.916378
  Spline 4 (from X[3]=2.301 to X[4]=3.301):
    S4(x) = -1.377388x³ + 13.072560x² + -40.867927x + 42.403729
  Spline 5 (from X[4]=3.301 to X[5]=4.301):
    S5(x) = 0.189278x³ + -2.442265x² + 10.346977x + -13.950249

  C_D(Re=5) = 2.70105 = f(0.699) = S2(0.699) 
  C_D(Re=50) = 2.72569 = f(1.6969) = S3(1.6969)
  C_D(Re=500) = 0.248586 = f(2.6989) = S4(2.6989)
  C_D(Re=5000) = 0.486384 = f(3.6989) = S5(3.6989)
*/

int main(int argc, char const *argv[]) {
    // Arrays for data points to read from text file
    double X[MAX_POINTS], Y[MAX_POINTS];
    // Arrays for polynomial coefficients calculation
    double A[MAX_SIZE+1][MAX_SIZE+1], b[MAX_SIZE+1], solution[MAX_SIZE+1];
    // Solution of Least Squares polynomial
    double a[MAX_SIZE+1];
    // n = number of data points and degree is the polynomial degree
    int n, degree;
    

    // Read data points from file
    if (!read_data_points("exercise2_log.txt", X, Y, &n)) {
        printf("Failed to read data from file. Exiting.\n");
        return 1;
    }
    
    // Print the data points
    print_data_points(X, Y, n);

    /* Initialize A and b to cero (to avoid trash data) */
    for (int i = 0; i < 4*(n-1); i++) {
        b[i] = 0.0;
        for (int j = 0; j < 4*(n-1); j++) {
            A[i][j] = 0.0;
        }
    }

    // We convert read data to logaritmic scale
    for (int i = 0; i < n; i++) {
      X[i] = log10(X[i]);
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
    int row1 = 4*(n-1) - 2;
    int row2 = 4*(n-1) - 1;
    
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

    // We evaluate the new reynolds of the exercise
    double Re[] = {5, 50, 500, 5000};
    for (int i = 0; i < 4; i++) {
        double logRe = log10(Re[i]); // Convert to logaritmic scale
        double Cd = evaluate_spline(logRe, X, solution, n);
        printf("C_D(Re=%g) = %g\n", Re[i], Cd);
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

double evaluate_spline(double x, double X[], double solution[], int n) {
    int k;

    // Search for the correspondient interval
    for (k = 0; k < n - 1; k++) {
        if (x >= X[k] && x <= X[k+1]) {
            break;
        }
    }

    if (k == n-1) {
        // If x is out of range, return 0 (or we can extrapolate)
        printf("Warning: x = %f fuera del rango de la spline\n", x);
        return 0.0;
    }

    // Evaluate correspondient cubic polynomial: S_k(x) = a*x^3 + b*x^2 + c*x + d
    double a = solution[4*k];
    double b = solution[4*k + 1];
    double c = solution[4*k + 2];
    double d = solution[4*k + 3];

    return a*pow(x, 3) + b*pow(x, 2) + c*x + d;
}