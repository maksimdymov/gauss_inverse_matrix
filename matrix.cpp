#include "matrix.h"
#include "my_errors.h"

int
InitMatr (double *a, int n, int s, FILE *f)
{
  if (s == 0)
    {
      int i = 0, j = 0, k = 0;
      for (i = 0; i < n; i++)
        {
          for (j = 0; j<n; j++)
            {
              k = fscanf (f, "%lf", a + n*i+j);
              if (k != 1)
                {
                  if (k == EOF)
                    {
                      return END_OF_FILE;
                    }
                  else
                    {
                      return FILE_READING_ERROR;
                    }
                }
            }
        }
      if (fscanf (f, "%d", &k) != EOF)
        {
          return TOO_MANY_INFO;
        }
      return SUCCESS;
    }
  for (int i = 0; i < n; i++)
    {
      for (int j = 0; j < n; j++)
        {
          a[n * i + j] = Aij (s, n, i + 1, j + 1);
        }
    }
  return SUCCESS;
}

double
Aij (int s, int n, int i, int j)
{
  double aij = 0; //(i, j) - element of matrix
  switch (s)
    {
      case 1:
        aij = n + 1. - (i < j ? j : i);
        break;

      case 2:
        aij = i < j ? j : i;
        break;

      case 3:
        aij = fabs (i - j);
        break;

      case 4:
        aij = 1. / (i + j - 1.);
        break;

      default:
        aij = (i == j);
    }
  return aij;
}

void
PrintMatr (const double *const a, int l, int k, int n)
{
  int min1 = ((l < n) ? l : n);
  int min2 = ((k < n) ? k : n);
  for (int i = 0; i < min1; i++)
    {
      for (int j = 0; j < min2; j++)
        {
          fprintf(stdout, " %10.3e", a[i * n + j]);
        }
      fprintf(stdout, "\n");
    }
  fprintf(stdout, "\n");
}