#ifndef SOLVE_H_
#define SOLVE_H_

#include <cstring>
#include <new>

int Ind (int ibl, int jbl, int iel, int jel, int n, int m);
void GetBlock (double *a, double *x, int ibl, int jbl, int n, int m);
void SetBlock (double *a, double *x, int ibl, int jbl, int n, int m);
double Residual (double* a, double* x, int n, int m, int pos_mul);
double Norm (double *a, int v, int h);
void Multi (double *x, double *y, double *res, int vx, int hx, int hy);
void Plus (double *x, double *y, double *res, int v, int h);
void Minus (double *x, double *y, double *res, int v, int h);

#endif