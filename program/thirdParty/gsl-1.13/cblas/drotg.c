#include <gsl/gsl_math.h>
#include <gsl/gsl_cblas.h>
#include "cblas.h"

void
cblas_drotg (double *a, double *b, double *c, double *s)
{
#define BASE double
#include "source_rotg.h"
#undef BASE
}
