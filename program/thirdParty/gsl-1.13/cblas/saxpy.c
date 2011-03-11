#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_saxpy (const int N, const float alpha, const float *X, const int incX,
             float *Y, const int incY)
{
#define BASE float
#include "source_axpy_r.h"
#undef BASE
}
