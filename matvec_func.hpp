#define N 10000
#define NITER 10

extern void matvec(const double a[N][N], const double b[N], double c[N]);
extern void matvec_cblas(const double a[N][N], const double b[N], double c[N]);
