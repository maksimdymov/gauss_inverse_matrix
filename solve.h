#ifndef SOLVE_H_
#define SOLVE_H_

#include <cstring>
#include <new>

int Ind (int ibl, int jbl, int iel, int jel, int n, int m, int *ind);
void GetBlock (double *src, double *dest, int ibl, int jbl, int n, int m, int *ind);
void SetBlock (double *dest, double *src, int ibl, int jbl, int n, int m, int *ind);
double Residual (double* a, double* x, int n, int m, int pos_mul);
double Norm (double *a, int v, int h);
void Multi (double *x, double *y, double *res, int vx, int hx, int hy);
void Plus (double *x, double *y, double *res, int v, int h);
void Minus (double *x, double *y, double *res, int v, int h);
int Inverse (double *matrix, double *inversed, int matrix_n, int matrix_norm, int *ind);
int Solve (double *matrix, double *inversed_matrix, int matrix_n, int block_m, int matrix_norm, int *ind);
void ReplaceLines (double *matrix, int matrix_n, int block_m, int *ind);

#endif