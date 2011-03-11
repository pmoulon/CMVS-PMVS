#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_srotmg (float *d1, float *d2, float *b1, const float b2, float *P)
{
#define BASE float
#include "source_rotmg.h"
#undef BASE
}
