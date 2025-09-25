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

/* Results for NACA 65:
Lower Surface:
Cubic Spline Coefficients:
    Spline 1 (from X[0]=0.000 to X[1]=0.658):
      S1(x) = 0.548527x³ + 0.000000x² + -1.468496x + -0.000000
    Spline 2 (from X[1]=0.658 to X[2]=0.920):
      S2(x) = -1.237194x³ + 3.525014x² + -3.787955x + 0.508735
    Spline 3 (from X[2]=0.920 to X[3]=1.441):
      S3(x) = -0.020776x³ + 0.167702x² + -0.699228x + -0.438475
    Spline 4 (from X[3]=1.441 to X[4]=2.717):
      S4(x) = -0.007779x³ + 0.111512x² + -0.618258x + -0.477367
    Spline 5 (from X[4]=2.717 to X[5]=5.243):
      S5(x) = -0.019466x³ + 0.206774x² + -0.877085x + -0.242956
    Spline 6 (from X[5]=5.243 to X[6]=7.753):
      S6(x) = 0.035650x³ + -0.660145x² + 3.668174x + -8.186555
    Spline 7 (from X[6]=7.753 to X[7]=10.254):
      S7(x) = -0.031888x³ + 0.910734x² + -8.510851x + 23.288105
    Spline 8 (from X[7]=10.254 to X[8]=15.243):
      S8(x) = 0.006159x³ + -0.259672x² + 3.490487x + -17.732467
    Spline 9 (from X[8]=15.243 to X[9]=20.219):
      S9(x) = -0.001716x³ + 0.100424x² + -1.998455x + 10.156846
    Spline 10 (from X[9]=20.219 to X[10]=25.189):
      S10(x) = 0.000472x³ + -0.032275x² + 0.684584x + -7.925942
    Spline 11 (from X[10]=25.189 to X[11]=30.154):
      S11(x) = -0.000174x³ + 0.016561x² + -0.545540x + 2.402587
    Spline 12 (from X[11]=30.154 to X[12]=35.116):
      S12(x) = 0.000088x³ + -0.007187x² + 0.170560x + -4.795169
    Spline 13 (from X[12]=35.116 to X[13]=40.077):
      S13(x) = -0.000128x³ + 0.015593x² + -0.629397x + 4.568587
    Spline 14 (from X[13]=40.077 to X[14]=46.088):
      S14(x) = 0.000318x³ + -0.038059x² + 1.520836x + -24.156365
    Spline 15 (from X[14]=46.088 to X[15]=50.000):
      S15(x) = -0.000664x³ + 0.097753x² + -4.738461x + 72.003130
    Spline 16 (from X[15]=50.000 to X[16]=54.965):
      S16(x) = 0.000790x³ + -0.120294x² + 6.163869x + -109.702380
    Spline 17 (from X[16]=54.965 to X[17]=59.986):
      S17(x) = -0.002451x³ + 0.414070x² + -23.207475x + 428.429596
    Spline 18 (from X[17]=59.986 to X[18]=64.914):
      S18(x) = 0.004915x³ + -0.911462x² + 56.305888x + -1161.466593
    Spline 19 (from X[18]=64.914 to X[19]=69.899):
      S19(x) = -0.004997x³ + 1.018667x² + -68.986474x + 1549.609540
    Spline 20 (from X[19]=69.899 to X[20]=74.898):
      S20(x) = 0.002539x³ + -0.561464x² + 41.463109x + -1023.828933
    Spline 21 (from X[20]=74.898 to X[21]=79.897):
      S21(x) = -0.000857x³ + 0.201496x² + -15.681113x + 402.833706
    Spline 22 (from X[21]=79.897 to X[22]=84.910):
      S22(x) = 0.000143x³ + -0.038062x² + 3.458913x + -106.909857
    Spline 23 (from X[22]=84.910 to X[23]=89.984):
      S23(x) = 0.000018x³ + -0.006187x² + 0.752386x + -30.306111
    Spline 24 (from X[23]=89.984 to X[24]=94.967):
      S24(x) = -0.000812x³ + 0.217844x² + -19.406871x + 574.364091
    Spline 25 (from X[24]=94.967 to X[25]=100.000):
      S25(x) = 0.000900x³ + -0.270059x² + 26.927824x + -892.391567

Upper Surface:
    Spline 1 (from X[0]=0.000 to X[1]=0.347):
      S1(x) = -4.196272x³ + 0.000000x² + 3.415932x + -0.000000
    Spline 2 (from X[1]=0.347 to X[2]=0.580):
      S2(x) = 5.567296x³ + -10.163875x² + 6.942796x + -0.407941
    Spline 3 (from X[2]=0.580 to X[3]=1.059):
      S3(x) = 0.381976x³ + -1.141417x² + 1.709771x + 0.603778
    Spline 4 (from X[3]=1.059 to X[4]=2.233):
      S4(x) = -0.050213x³ + 0.231646x² + 0.255697x + 1.117066
    Spline 5 (from X[4]=2.233 to X[5]=4.757):
      S5(x) = 0.015518x³ + -0.208685x² + 1.238956x + 0.385193
    Spline 6 (from X[5]=4.757 to X[6]=7.247):
      S6(x) = -0.003923x³ + 0.068755x² + -0.080824x + 2.477924
    Spline 7 (from X[6]=7.247 to X[7]=9.746):
      S7(x) = 0.001381x³ + -0.046553x² + 0.754813x + 0.459304
    Spline 8 (from X[7]=9.746 to X[8]=14.757):
      S8(x) = 0.000139x³ + -0.010255x² + 0.401055x + 1.608545
    Spline 9 (from X[8]=14.757 to X[9]=19.781):
      S9(x) = -0.000212x³ + 0.005302x² + 0.171480x + 2.737823
    Spline 10 (from X[9]=19.781 to X[10]=24.811):
      S10(x) = 0.001185x³ + -0.077586x² + 1.811072x + -8.073100
    Spline 11 (from X[10]=24.811 to X[11]=29.846):
      S11(x) = -0.001874x³ + 0.150037x² + -3.836483x + 38.634069
    Spline 12 (from X[11]=29.846 to X[12]=39.028):
      S12(x) = 0.000808x³ + -0.090072x² + 3.329820x + -32.661099
    Spline 13 (from X[12]=39.028 to X[13]=44.962):
      S13(x) = -0.000536x³ + 0.067317x² + -2.812761x + 47.249788
    Spline 14 (from X[13]=44.962 to X[14]=50.000):
      S14(x) = 0.000080x³ + -0.015742x² + 0.921729x + -8.720254
    Spline 15 (from X[14]=50.000 to X[15]=55.035):
      S15(x) = 0.000059x³ + -0.012640x² + 0.766640x + -6.135441
    Spline 16 (from X[15]=55.035 to X[16]=60.064):
      S16(x) = 0.000024x³ + -0.006863x² + 0.448705x + -0.302925
    Spline 17 (from X[16]=60.064 to X[17]=65.086):
      S17(x) = 0.000021x³ + -0.006296x² + 0.414675x + 0.378403
    Spline 18 (from X[17]=65.086 to X[18]=70.101):
      S18(x) = 0.000039x³ + -0.009801x² + 0.642752x + -4.569803
    Spline 19 (from X[18]=70.101 to X[19]=75.107):
      S19(x) = 0.000008x³ + -0.003380x² + 0.192674x + 5.947171
    Spline 20 (from X[19]=75.107 to X[20]=80.108):
      S20(x) = 0.000026x³ + -0.007470x² + 0.499816x + -1.742335
    Spline 21 (from X[20]=80.108 to X[21]=85.090):
      S21(x) = 0.000055x³ + -0.014262x² + 1.043922x + -16.271428
    Spline 22 (from X[21]=85.090 to X[22]=90.060):
      S22(x) = -0.000117x³ + 0.029542x² + -2.683335x + 89.446027
    Spline 23 (from X[22]=90.060 to X[23]=95.088):
      S23(x) = 0.000363x³ + -0.100248x² + 9.005543x + -261.454108
    Spline 24 (from X[23]=95.088 to X[24]=100.000):
      S24(x) = -0.000233x³ + 0.069852x² + -7.168957x + 251.212841

Results for NACA 66:
Lower Surface:
Cubic Spline Coefficients:
    Spline 1 (from X[0]=0.000 to X[1]=0.611):
      S1(x) = 1.156580x³ + 0.000000x² + -2.507062x + -0.000000
    Spline 2 (from X[1]=0.611 to X[2]=0.872):
      S2(x) = -2.820897x³ + 7.290716x² + -6.961690x + 0.907259
    Spline 3 (from X[2]=0.872 to X[3]=1.385):
      S3(x) = 0.170404x³ + -0.534527x² + -0.138077x + -1.076137
    Spline 4 (from X[3]=1.385 to X[4]=2.654):
      S4(x) = -0.041911x³ + 0.347638x² + -1.359876x + -0.512074
    Spline 5 (from X[4]=2.654 to X[5]=5.178):
      S5(x) = 0.000457x³ + 0.010306x² + -0.464598x + -1.304097
    Spline 6 (from X[5]=5.178 to X[6]=7.680):
      S6(x) = -0.001313x³ + 0.037803x² + -0.606975x + -1.058353
    Spline 7 (from X[6]=7.680 to X[7]=10.182):
      S7(x) = -0.000325x³ + 0.015049x² + -0.432228x + -1.505707
    Spline 8 (from X[7]=10.182 to X[8]=15.175):
      S8(x) = 0.000306x³ + -0.004226x² + -0.235970x + -2.171807
    Spline 9 (from X[8]=15.175 to X[9]=20.159):
      S9(x) = -0.000964x³ + 0.053564x² + -1.112926x + 2.264128
    Spline 10 (from X[9]=20.159 to X[10]=25.137):
      S10(x) = 0.000847x³ + -0.055950x² + 1.094772x + -12.570862
    Spline 11 (from X[10]=25.137 to X[11]=30.113):
      S11(x) = -0.000464x³ + 0.042915x² + -1.390410x + 8.252478
    Spline 12 (from X[11]=30.113 to X[12]=35.086):
      S12(x) = 0.000134x³ + -0.011069x² + 0.235227x + -8.065126
    Spline 13 (from X[12]=35.086 to X[13]=40.058):
      S13(x) = -0.000076x³ + 0.010952x² + -0.537431x + 0.971360
    Spline 14 (from X[13]=40.058 to X[14]=45.029):
      S14(x) = 0.000113x³ + -0.011700x² + 0.369997x + -11.145215
    Spline 15 (from X[14]=45.029 to X[15]=50.000):
      S15(x) = -0.000253x³ + 0.037728x² + -1.855720x + 22.262052
    Spline 16 (from X[15]=50.000 to X[16]=54.972):
      S16(x) = 0.000508x³ + -0.076481x² + 3.854756x + -72.912543
    Spline 17 (from X[16]=54.972 to X[17]=59.964):
      S17(x) = -0.000148x³ + 0.031677x² + -2.090919x + 36.035992
    Spline 18 (from X[17]=59.964 to X[18]=64.925):
      S18(x) = -0.000027x³ + 0.009918x² + -0.786175x + 9.956770
    Spline 19 (from X[18]=64.925 to X[19]=69.911):
      S19(x) = -0.000145x³ + 0.033065x² + -2.288994x + 42.480285
    Spline 20 (from X[19]=69.911 to X[20]=74.905):
      S20(x) = -0.000036x³ + 0.010178x² + -0.688928x + 5.192891
    Spline 21 (from X[20]=74.905 to X[21]=79.907):
      S21(x) = -0.000171x³ + 0.040538x² + -2.963036x + 61.973562
    Spline 22 (from X[21]=79.907 to X[22]=84.919):
      S22(x) = -0.000090x³ + 0.020996x² + -1.401481x + 20.380521
    Spline 23 (from X[22]=84.919 to X[23]=89.940):
      S23(x) = 0.000034x³ + -0.010506x² + 1.273639x + -55.342331
    Spline 24 (from X[23]=89.940 to X[24]=94.970):
      S24(x) = -0.000946x³ + 0.253765x² + -22.494971x + 657.240598
    Spline 25 (from X[24]=94.970 to X[25]=100.000):
      S25(x) = 0.001038x³ + -0.311329x² + 31.172038x + -1041.678014

Upper Surface:
    Cubic Spline Coefficients:
    Spline 1 (from X[0]=0.000 to X[1]=0.389):
      S1(x) = -12.898520x³ + 0.000000x² + 6.753874x + -0.000000
    Spline 2 (from X[1]=0.389 to X[2]=0.628):
      S2(x) = 30.259118x³ + -50.364964x² + 26.345845x + -2.540426
    Spline 3 (from X[2]=0.628 to X[3]=1.115):
      S3(x) = -5.335922x³ + 16.696092x² + -15.768499x + 6.275510
    Spline 4 (from X[3]=1.115 to X[4]=2.846):
      S4(x) = 0.293590x³ + -2.134627x² + 5.227753x + -1.528097
    Spline 5 (from X[4]=2.846 to X[5]=4.827):
      S5(x) = -0.090677x³ + 1.146250x² + -4.109623x + 7.329961
    Spline 6 (from X[5]=4.827 to X[6]=7.320):
      S6(x) = 0.030811x³ + -0.613021x² + 4.382375x + -6.333663
    Spline 7 (from X[6]=7.320 to X[7]=0.818):
      S7(x) = -0.001444x³ + 0.095299x² + -0.802523x + 6.317488
    Spline 8 (from X[7]=0.818 to X[8]=14.825):
      S8(x) = -0.002774x³ + 0.098564x² + -0.805194x + 6.318216
    Spline 9 (from X[8]=14.825 to X[9]=19.841):
      S9(x) = 0.001254x³ + -0.080593x² + 1.850809x + -6.806864
    Spline 10 (from X[9]=19.841 to X[10]=24.868):
      S10(x) = 0.001883x³ + -0.118057x² + 2.594125x + -11.722908
    Spline 11 (from X[10]=24.868 to X[11]=29.887):
      S11(x) = -0.004560x³ + 0.362663x² + -9.360426x + 87.372348
    Spline 12 (from X[11]=29.887 to X[12]=34.914):
      S12(x) = 0.004676x³ + -0.465480x² + 15.390296x + -159.202599
    Spline 13 (from X[12]=34.914 to X[13]=39.942):
      S13(x) = -0.002284x³ + 0.263493x² + -10.061065x + 137.000350
    Spline 14 (from X[13]=39.942 to X[14]=44.971):
      S14(x) = 0.000581x³ + -0.079750x² + 3.648718x + -45.531713
    Spline 15 (from X[14]=44.971 to X[15]=50.000):
      S15(x) = -0.000142x³ + 0.017783x² + -0.737398x + 20.217640
    Spline 16 (from X[15]=50.000 to X[16]=55.028):
      S16(x) = -0.000029x³ + 0.000788x² + 0.112333x + 6.055452
    Spline 17 (from X[16]=55.028 to X[17]=60.054):
      S17(x) = -0.000162x³ + 0.022712x² + -1.094105x + 28.184741
    Spline 18 (from X[17]=60.054 to X[18]=65.075):
      S18(x) = -0.000380x³ + 0.062109x² + -3.460041x + 75.546054
    Spline 19 (from X[18]=65.075 to X[19]=70.089):
      S19(x) = 0.002404x³ + -0.481391x² + 31.908187x + -691.649756
    Spline 20 (from X[19]=70.089 to X[20]=75.095):
      S20(x) = -0.004667x³ + 1.005384x² + -72.298392x + 1742.928546
    Spline 21 (from X[20]=75.095 to X[21]=80.098):
      S21(x) = 0.004791x³ + -1.125411x² + 87.713706x + -2262.440964
    Spline 22 (from X[21]=80.098 to X[22]=85.081):
      S22(x) = -0.001996x³ + 0.505475x² + -42.917049x + 1225.313117
    Spline 23 (from X[22]=85.081 to X[23]=90.060):
      S23(x) = 0.000155x³ + -0.043458x² + 3.786721x + -99.221362
    Spline 24 (from X[23]=90.060 to X[24]=95.030):
      S24(x) = 0.000870x³ + -0.236747x² + 21.194352x + -621.798459
    Spline 25 (from X[24]=95.030 to X[25]=100.000):
      S25(x) = -0.000757x³ + 0.227173x² + -22.891907x + 774.707279
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
    if (!read_data_points("upper_surface_NACA66.txt", X, Y, &n)) {
        printf("Failed to read data from file. Exiting.\n");
        return 1;
    }
    
    // Print the data points
    print_data_points(X, Y, n);

    /* Inicializar A y b a cero (importante) */
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
