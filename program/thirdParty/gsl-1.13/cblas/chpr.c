#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_chpr (const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
            const int N, const float alpha, const void *X, const int incX,
            void *Ap)
{
#define BASE float
#include "source_hpr.h"
#undef BASE
}
