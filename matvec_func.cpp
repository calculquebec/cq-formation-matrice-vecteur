#include <cblas.h>
#include "matvec_func.hpp"

void matvec(const double a[N][N], const double b[N], double c[N])
{
  int i, j;
  for (i = 0; i < N; i++)
    c[i] = 0;

  for (j = 0; j < N; j++)
    for (i = 0; i < N; i++)
      c[i] += a[i][j] * b[j];
}

void matvec_cblas(const double a[N][N], const double b[N], double c[N])
{
  for (int i = 0; i < N; i++)
    c[i] = 0;

  cblas_dgemv(CblasRowMajor, CblasNoTrans,
	      N, N, 1.0, &a[0][0], N, b, 1, 0, c, 1);
}
