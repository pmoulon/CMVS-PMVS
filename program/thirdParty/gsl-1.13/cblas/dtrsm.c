#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_dtrsm (const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
             const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
             const enum CBLAS_DIAG Diag, const int M, const int N,
             const double alpha, const double *A, const int lda, double *B,
             const int ldb)
{
#define BASE double
#include "source_trsm_r.h"
#undef BASE
}
