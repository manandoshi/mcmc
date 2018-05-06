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
    memcpy(result, B, n*n);
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
    memcpy(B_pow, B, n*n);
    memcpy(result, Id, n*n);
    double coeff = 1.0;
    for(int j=0; j<iter; j++)
    {
        coeff = -1*coeff*(d - j)/(j + 1);
        cblas_daxpy(n, coeff, B_pow, 1, result, 1);
        if(j == iter - 1) break;
        cblas_dtrmm(CblasRowMajor, CblasRight, CblasLower, CblasNoTrans, CblasNonUnit, n, n, 1, B, n, B_pow, n);
    }
    delete B_pow;
    return;
}

double calc_MLE(double *y, int n, float var, float d, float phi_1)
{
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
    memcpy(z, y, n);
    cblas_dtrmv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit, n, result, n, z, 1);

    double mle = 1;
    for(int i = 0; i < n; i++)
    {
        mle = mle * 1/sqrt(2*M_PI*var)*exp(-1*z[i]*z[i]/var);
    }
    
    return mle;
}

int main()
{
    //float c,d,e;
    //d = 1.0;
    //calc_MLE(a,3,c,d,e);
    //double a[] = {1,0,4,1};
    //double b[] = {1,0,8,4};
    //double x[] = {1,2}; 
    //double *c = new double[4];
    //double *result = new double[4];
    //for(int i = 0; i < 4; i++) result[i] = 0;
    //int info_a, info_b, info; info_a = info_b = 2;
    //cblas_daxpy(4, 1, a, 1, result, 1);
    //cblas_daxpy(4, 1, b, 1, result, 1);
    //cblas_daxpy(4, 1, b, 1, c, 1);
    //cblas_dtrmm(CblasRowMajor, CblasLeft, CblasLower, CblasNoTrans, CblasUnit, 2, 2, 1, a, info_a, c, info_b);
    //for(int i = 0; i < 4; i++)
    //    cout<<result[i]<<endl;
    //info = LAPACKE_dlascl(LAPACK_ROW_MAJOR, 'L', 1, 1, 1, 2.3, 2, 2, c, 2);
    //cout<<"Info:"<<info;
    //for(int i = 0; i < 2; i++)
    //{
    //    for(int j = 0; j < 2; j++)
    //    {
    //        cout<<c[j + 2*i]<<" ";
    //    }
    //    cout<<endl;
    //}
    //cblas_dtrmv(CblasRowMajor, CblasLower, CblasNoTrans, CblasNonUnit, 2, a, 2, x, 1);
    //for(int i = 0; i < 2; i++)
    //{
    //    cout<<x[i]<<endl;
    //}

    ifstream ifile("data.csv");
    string value;
    while(ifile.good())
    {
        getline(ifile, value, ',');
        cout<<value<<endl;
    }    
    return 0;
}
