#include "matrix.h"
#include "my_errors.h"
#include "solve.h"

#define EPS 1.110223e-16

int
Ind (int ibl, int jbl, int iel, int jel, int n, int m, int *ind)
{
  if (ind)
    {
      return n * ind[m * ibl + iel] + m * jbl + jel;
    }
  return n * (m * ibl + iel) + m * jbl + jel;
}

void
GetBlock (double *a, double *x, int ibl, int jbl, int n, int m, int *ind)
{
  int v = 0; //Number of lines
  int h = 0; //Number of columns
  int l = n % m, k = n / m;
  v = (ibl < k ? m : l);
  h = (jbl < k ? m : l);
  for (int i = 0; i < v; i++)
    {
      memcpy (x + i * h, a + Ind (ibl, jbl, i, 0, n, m, ind), h * sizeof (double));
    }
}

void
SetBlock (double *a, double *x, int ibl, int jbl, int n, int m, int *ind)
{
  int v = 0, h = 0, l = n % m, k = n / m;
  v = (ibl < k ? m : l);
  h = (jbl < k ? m : l);
  for (int i = 0; i < v; i++)
    {
      memcpy (a + Ind (ibl, jbl, i, 0, n, m, ind), x + i * h, h * sizeof (double));
    }
}

void
Plus (double *x, double *y, double *res, int v, int h)
{
  for (int i = 0; i < v; i++)
    {
      for (int j = 0; j < h; j++)
        {
          res[i * h + j] = x[i * h + j] + y[i * h + j];
        }
    }
}

void
Minus (double *x, double *y, double *res, int v, int h)
{
  for (int i = 0; i < v; i++)
    {
      for (int j = 0; j < h; j++)
        {
          res[i * h + j] = x[i * h + j] - y[i * h + j];
        }
    }
}

void
Multi (double *x, double *y, double *res, int vx, int hx, int hy)
{
  int v3 = vx % 3, h3 = hy % 3, i = 0, j = 0;
  double sum = 0;
  for (i = 0; i < v3; i++)
    {
      for (j = 0; j < h3; j++)
        {
          sum = 0;
          for (int s = 0; s < hx; s++)
            {
              sum += x[i * hx + s] * y[s * hy + j];
            }
          res[i * hy + j] = sum;
        }
      for (; j < hy; j += 3)
        {
          double s00 = 0, s01 = 0, s02 = 0;
          for (int s = 0; s < hx; s++)
            {
              s00 += x[i * hx + s] * y[s * hy + j];
              s01 += x[i * hx + s] * y[s * hy + j + 1];
              s02 += x[i * hx + s] * y[s * hy + j + 2];
            }
          res[i * hy + j] = s00;
          res[i * hy + j + 1] = s01;
          res[i * hy + j + 2] = s02;
        }
    }
  for (; i < vx; i += 3)
    {
      for (j = 0; j < h3; j++)
        {
          double s00 = 0, s10 = 0, s20 = 0;
          for (int s = 0; s < hx; s++)
            {
              s00 += x[i * hx + s] * y[s * hy + j];
              s10 += x[(i + 1) * hx + s] * y[s * hy + j];
              s20 += x[(i + 2) * hx + s] * y[s * hy + j];
            }
          res[i * hy + j] = s00;
          res[(i + 1) * hy + j] = s10;
          res[(i + 2) * hy + j] = s20;
        }
      for (; j < hy; j += 3)
        {
          double s00 = 0, s01 = 0, s02 = 0, s10 = 0, s11 = 0, s12 = 0, s20 = 0, s21 = 0, s22 = 0;
          for (int s = 0; s < hx; s++)
            {
              s00 += x[i * hx + s] * y[s * hy + j];
              s01 += x[i * hx + s] * y[s * hy + j + 1];
              s02 += x[i * hx + s] * y[s * hy + j + 2];
              s10 += x[(i + 1) * hx + s] * y[s * hy + j];
              s11 += x[(i + 1) * hx + s] * y[s * hy + j + 1];
              s12 += x[(i + 1) * hx + s] * y[s * hy + j + 2];
              s20 += x[(i + 2) * hx + s] * y[s * hy + j];
              s21 += x[(i + 2) * hx + s] * y[s * hy + j + 1];
              s22 += x[(i + 2) * hx + s] * y[s * hy + j + 2];
            }
          res[i * hy + j] = s00;
          res[i * hy + j + 1] = s01;
          res[i * hy + j + 2] = s02;
          res[(i + 1) * hy + j] = s10;
          res[(i + 1) * hy + j + 1] = s11;
          res[(i + 1) * hy + j + 2] = s12;
          res[(i + 2) * hy + j] = s20;
          res[(i + 2) * hy + j + 1] = s21;
          res[(i + 2) * hy + j + 2] = s22;
        }
    }
}


double
Norm (double *a, int v, int h)
{
  double norm = 0;
  for (int i = 0; i < v; i++)
    {
      double sum = 0;
      for (int j = 0; j < h; j++)
        {
          sum += fabs (a[i * h + j]);
        }
      if (sum > norm)
        {
          norm = sum;
        }
    }
  return norm;
}

double
Residual (double* a, double* x, int n, int m, int pos_mul)
{
  double res = 0; //Residual
  int l = n % m;
  int k = n / m;
  double* mul = a + pos_mul * n * n; //Product
  double* d1 = nullptr;
  double* d2 = nullptr;
  double* d3 = nullptr;
  double norm = 0; //Norm of block line of product
  int *ind = nullptr;

  memset(mul, 0, (n + 3 * m) * m * sizeof (double));
  d1 = mul + n * m;
  d2 = d1 + m * m;
  d3 = d2 + m * m;
  
  for (int i = 0; i < (l ? k + 1 : k); i++)
    {
      int vx = (i < k) ? m : l; //Number of lines in first multiplier
      for (int j = 0; j < (l ? k + 1 : k); j++)
        {
          for (int s = 0; s < (l ? k + 1 : k); s++)
            {
              int hx = (s < k) ? m : l; //Number of columns in first multiplier
              int hy = (j < k) ? m : l; //Number of columns in second multiplier
              GetBlock (a, d1, i, s, n, m, ind);
              GetBlock (x, d2, s, j, n, m, ind);
              Multi (d1, d2, d3, vx, hx, hy);
              GetBlock (mul, d1, 0, j, n, m, ind);
              Plus (d1, d3, d2, vx, hy);
              SetBlock (mul, d2, 0, j, n, m, ind);
            }
        }
      for (int j = 0; j < vx; j++)
        {
          mul[Ind (0, i, j, j, n, m, ind)]--;
        }
      norm = Norm (mul, vx, n);
      if (norm > res)
        {
          res = norm;
        }
      memset (mul, 0, n * m * sizeof (double));
    }
  return res;
}

//Перед обращением надо правильно заполнить массив перестановки индексов
int
Inverse (double *matrix, double *inversed, int matrix_n, int matrix_norm, int *ind)
{
  memset (inversed, 0, matrix_n * matrix_n * sizeof (double));
  for (int i = 0; i < matrix_n; i++)
    {
      inversed[i * matrix_n + i] = 1.;
    }
  //Прямой ход
  for (int s = 0; s < matrix_n; s++)
    {
      //Нахождение главного
      double max = fabs (matrix[ind[s] * matrix_n + s]);
      double reverse = 0;
      int n_max = s;
      int swap = 0;
      int i = 0;
      for (i = s + 1; i < matrix_n; i++)
        {
          if (fabs (matrix[ind[i] * matrix_n + s]) > max)
            {
              max = fabs (matrix[ind[i] * matrix_n + s]);
              n_max = i;
            }
        }
      swap = ind[s];
      ind[s] = ind[n_max];
      ind[n_max] = swap;
      
      //Домножаем строку на обратный к первому элементу
      if (fabs (matrix[ind[s] * matrix_n + s]) < matrix_norm * EPS)
        {
          return -1;
        }
      reverse = 1. / matrix[ind[s] * matrix_n + s];
      matrix[ind[s] * matrix_n + s] = 1.;
      for (i = 0; i < s + 1; i++)
        {
          inversed[ind[s] * matrix_n + i] *= reverse;
        }
      for (int i = s + 1; i < matrix_n; i++)
        {
          inversed[ind[s] * matrix_n + i] *= reverse;
          matrix[ind[s] * matrix_n + i] *= reverse;
        }
      
      //Из каждой нижележащей строки вычитаем первую, домноженную на первый элемент строки
      for (i = s + 1; i < matrix_n; i++)
        {
          int j = 0;
          
          for (j = 0; j < s + 1; j++)
            {
              inversed[ind[i] * matrix_n + j] -= matrix[ind[i] * matrix_n + s] * inversed[ind[s] * matrix_n + j];
            }
          for (; j < matrix_n; j++)
            {
              inversed[ind[i] * matrix_n + j] -= matrix[ind[i] * matrix_n + s] * inversed[ind[s] * matrix_n + j];
              matrix[ind[i] * matrix_n + j] -= matrix[ind[i] * matrix_n + s] * matrix[ind[s] * matrix_n + j];
            }
          //matrix[ind[i] * matrix_n + s] = 0.;
        }
    }
  for (int s = matrix_n - 1; s >= 0; s--)
    {
      //Для всех вышележащих вычитаем последнюю строку, домноженную на s-ый элемент строки
      for (int i = s - 1; i >= 0; i--)
        {
          for (int j = 0; j < matrix_n; j++)
            {
              inversed[ind[i] * matrix_n + j] -= inversed[ind[s] * matrix_n + j] * matrix[ind[i] * matrix_n + s];
            }
          //matrix[ind[i] * matrix_n + s] = 0.;
        }
    }
  return 0;
}