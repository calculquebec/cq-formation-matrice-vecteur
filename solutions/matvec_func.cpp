#include <cblas.h>
#include "matvec_func.hpp"

#if SLOW
void matvec(const double a[N][N], const double b[N], double c[N])
{
  int i, j;
  for (i = 0; i < N; i++)
    c[i] = 0;

  for (j = 0; j < N; j++)
    for (i = 0; i < N; i++)
      c[i] += a[i][j] * b[j];
}
#endif

#if TRANSPOSE
void matvec(const double a[N][N], const double b[N], double c[N])
{
  int i, j;
  for (i = 0; i < N; i++) {
    double ci = 0;
    for (j = 0; j < N; j++)
      ci += a[i][j] * b[j];
    c[i] = ci;
  }
}
#endif

#if BLOCKED
void matvec(const double a[N][N], const double b[N], double c[N])
{
  int i, j;
  for (i = 0; i < N; i+=4) {
    double ci0 = 0, ci1 = 0, ci2 = 0, ci3 = 0;
    for (j = 0; j < N; j++) {
      ci0 += a[i][j] * b[j];
      ci1 += a[i+1][j] * b[j];
      ci2 += a[i+2][j] * b[j];
      ci3 += a[i+3][j] * b[j];
    }
    c[i] = ci0;
    c[i+1] = ci1;
    c[i+2] = ci2;
    c[i+3] = ci3;
  }
}
#endif

void matvec(const double a[N][N], const double b[N], double c[N])
{
  int i, j;
  for (i = 0; i < N; i+=4) {
    double ci0 = 0, ci1 = 0, ci2 = 0, ci3 = 0;
#pragma omp simd reduction(+:ci0,ci1,ci2,ci3)
    for (j = 0; j < N; j++) {
      ci0 += a[i][j] * b[j];
      ci1 += a[i+1][j] * b[j];
      ci2 += a[i+2][j] * b[j];
      ci3 += a[i+3][j] * b[j];
    }
    c[i] = ci0;
    c[i+1] = ci1;
    c[i+2] = ci2;
    c[i+3] = ci3;
  }
}

void matvec_cblas(const double a[N][N], const double b[N], double c[N])
{
  for (int i = 0; i < N; i++)
    c[i] = 0;

  cblas_dgemv(CblasRowMajor, CblasNoTrans,
	      N, N, 1.0, &a[0][0], N, b, 1, 0, c, 1);
}
