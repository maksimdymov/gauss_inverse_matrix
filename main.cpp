#include <ctime>

#include "matrix.h"
#include "my_errors.h"
#include "solve.h"

int main(int argc, char * argv[])
{
  int n = 0; //Size of matrix
  int m = 0; //Size of block
  int r = 0; //Size of printed matrix
  int s = 0; //Number of formulae
  int err = 0;
  FILE *fi = nullptr; //Input file
  double* a = nullptr; //Pointer of matrix
  double* x = nullptr; //Pointer of attached matrix
  double time1 = 0, time2 = 0, time3 = 0, time4 = 0;
  double elapsed = 0; //Time of calculating inversed matrix
  double residual_time = 0; //Time of calculating residual
  double r1 = 0; //Residual 1
  double r2 = 0; //Residual 2

  if ((argc != 5) && (argc != 6))
    {
      return PrintErrorMsgByCode (MAIN_ARGS_ERROR, argv[0]);
    }
  if ((!sscanf (argv[1], "%d", &n)) || (!sscanf (argv[2], "%d", &m)) || (!sscanf (argv[3], "%d", &r)) || (!sscanf (argv[4], "%d", &s)))
    {
      return PrintErrorMsgByCode (MAIN_ARGS_ERROR, argv[0]);
    }
  if (((argc == 6) && (s != 0)) || ((s == 0) && (argc == 5)))
    {
      return PrintErrorMsgByCode (MAIN_ARGS_ERROR, argv[0]);
    }
  if ((s < 0) || (s > 4))
    {
      return PrintErrorMsgByCode (MAIN_ARGS_ERROR, argv[0]);
    }
  if (s == 0)
    {
      fi = fopen (argv[5], "r");
      if (fi == NULL)
        {
          return PrintErrorMsgByCode (FILE_OPENING_ERROR, argv[5]);
        }
    }
  try
    {
      a = new double[2 * n * n + (n + 3 * m) * m];
    }
  catch (std::bad_alloc& e)
    {
      char msg[100];
      sprintf (msg, "%d", 2 * n * n + (n + 3 * m) * m);
      return PrintErrorMsgByCode (ALLOCATE_MEMORY_ERROR, msg);
    }
  if (a == nullptr)
    {
      char msg[100];
      sprintf (msg, "%d", 2 * n * n + (n + 3 * m) * m);
      return PrintErrorMsgByCode (ALLOCATE_MEMORY_ERROR, msg);
    }
  x = a + n * n;
  err = InitMatr (a, n, s, fi);
  if (err != SUCCESS) {
      delete[] a;
      return PrintErrorMsgByCode (err, argv[5]);
    }
  err = InitMatr (x, n, 10, fi);
  fprintf (stdout, "Initial matrix:\n");
  PrintMatr (a, r, r, n);
  //double bignorm = Norm (a, n, n); //Norm of the whole matrix
  time1 = clock ();
  //err = Solve(a, b, x, n, m, bignorm);
  time2 = clock ();
  /*if (err)
    {
      cerr << "Cannot solve system" << endl;
      delete[] a;
      return -1;
    }*/
  fprintf (stdout, "Inversed matrix:\n");
  PrintMatr (x, r, r, n);
  if (s == 0)
    {
      rewind (fi);
    }
  InitMatr (a, n, s, fi);
  if (s == 0) {
      fclose (fi);
    }
  if (n <= 11000)
    {
      time3 = clock ();
      r1 = Residual (a, x, n, m, 2);
      r2 = Residual (x, a, n, m, 1);
      time4 = clock ();
    }
  delete[] a;
  elapsed = (time2 - time1) / CLOCKS_PER_SEC;
  residual_time = (time4 - time3) / CLOCKS_PER_SEC;
  //printf ("%e\n", bignorm);
  printf ("%s : Task = %d Res1 = %e Res2 = %e T1 = %.2f T2 = %.2f S = %d N = %d M = %d\n", argv[0], 12, r1, r2, elapsed, residual_time, s, n, m);
  return 0;
}