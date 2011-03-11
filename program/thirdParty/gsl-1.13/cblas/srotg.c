#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_srotg (float *a, float *b, float *c, float *s)
{
#define BASE float
#include "source_rotg.h"
#undef BASE
}
