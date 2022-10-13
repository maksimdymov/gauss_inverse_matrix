#ifndef MATRIX_H_
#define MATRIX_H_

#include <cmath>
#include <cstdio>

int InitMatr (double *a, int n, int s, FILE *f);
double Aij (int s, int n, int i, int j);
void PrintMatr (const double *const a, int l, int k, int n);

#endif