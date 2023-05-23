#include <iostream>
#include <omp.h>
#include "matvec_func.hpp"

double a[N][N];
double b[N], c[N];

int main(void)
{
  double t1, sum;
  int i, j;
  
  for (i = 0; i < N; i++)
    for (j = 0; j < N; j++)
      a[i][j] = 1.0;

  for (i = 0; i < N; i++)
    b[i] = 1.0;

  t1 = omp_get_wtime();
  for (i = 0; i < NITER; i++)
    matvec(a, b, c);
  t1 = omp_get_wtime() - t1;

  sum = 0;
  for (i = 0; i < N; i++)
    sum += c[i];

  if (sum == N*N)
    std::cout << "correct, time=" << t1/NITER << std::endl;
  else
    std::cout << "FAIL\n";

  t1 = omp_get_wtime();
  for (i = 0; i < NITER; i++)
    matvec_cblas(a, b, c);
  t1 = omp_get_wtime() - t1;

  sum = 0;
  for (i = 0; i < N; i++)
    sum += c[i];

  if (sum == N*N)
    std::cout << "CBLAS: correct, time=" << t1/NITER << std::endl;
  else
    std::cout << "CBLAS: FAIL\n";
}
