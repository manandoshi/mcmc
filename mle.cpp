#include <iostream>
#include <cstring>
#include <stdlib.h>
#include <lapacke.h>
#include <cblas.h>
#include <math.h>
#include <fstream>
#include <string.h>
#include <vector>

using namespace std;

void AR(int n, double* Id, double* B, double* result, float phi_1)
{
    int info;
    memcpy(result, B, n*n*sizeof(double));
    info = LAPACKE_dlascl(LAPACK_ROW_MAJOR, 'L', 1, 1, -1, phi_1, n, n, result, n);
    if(info != 0)
    {
        cout<<"Error in scaling of B";
        exit;
    }
    cblas_daxpy(n*n, 1, Id, 1, result, 1);
    return;
}

void differential(int n, double *Id, double *B, double *result, float d)
{
    int iter = 0;
    if(ceilf(d) == d)
    {
        iter = fmin(d,n);
        cout<<"Integer hai re"<<endl;
    }
    else
    {
        iter = n;
        cout<<"Not integer re"<<endl;
    }
    double *B_pow = new double[n*n];
    memcpy(B_pow, B, n*n*sizeof(double));
    memcpy(result, Id, n*n*sizeof(double));
    double coeff = 1.0;
    for(int j=0; j<iter; j++)
    {
        coeff = -1*coeff*(d - j)/(j + 1);
        cblas_daxpy(n*n, coeff, B_pow, 1, result, 1);
        if(j == iter - 1) break;
        cblas_dtrmm(CblasRowMajor, CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, n, n, 1, B, n, B_pow, n);
    }
    delete B_pow;
    return;
}

double calc_MLE(int n, float var, float d, float phi_1, string fname)
{
    //Getting the data
    ifstream ifile(fname);
    string value;
    double y[n];
    int i = 0;
    while(getline(ifile, value))
    {
        y[i] = stod(value);
        i++;
    }

    double* Id = new double[n*n];
    double* B = new double[n*n];
    for(int i = 0; i < n*n; i++)
    {
        Id[i] = 0;
        B[i] = 0;
    }
    Id[0] = 1;
    for(int i = 1; i < n; i++)
    {
        Id[i*n + i] = 1;
        B[i*n + i - 1] = 1;
    }
    double *temp = new double[n*n];
    double *result = new double[n*n];
    double *z = new double[n];
    //Load differential to result
    differential(n, Id, B, result, d);
    //Load AR to temp
    AR(n, Id, B, temp, phi_1);
    //Compute AR*differential
    cblas_dtrmm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasNonUnit, n, n, 1, temp, n, result, n);

    //Multiply the resulting operator with y to get z
    memcpy(z, y, n*sizeof(double));
    cblas_dtrmv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit, n, result, n, z, 1);

    double nll = 0;
    for(int i = 0; i < n; i++)
    {
        double factor = z[i]*z[i]/var;
        cout<<z[i]<<endl;
        nll += factor;
    }
    delete temp;
    delete result;
    delete z;

    return nll;
}
