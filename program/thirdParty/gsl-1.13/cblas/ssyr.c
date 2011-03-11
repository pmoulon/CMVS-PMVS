#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_ssyr (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
            const int N, const float alpha, const float *X, const int incX,
            float *A, const int lda)
{
#define BASE float
#include "source_syr.h"
#undef BASE
}
