#include <config.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_complex.h>

#include <gsl/gsl_fft_complex.h>
#include <gsl/gsl_fft_complex_float.h>

#define BASE_DOUBLE
#include "templates_on.h"
#include "bitreverse.c"
#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
#include "bitreverse.c"
#include "templates_off.h"
#undef  BASE_FLOAT

#include "factorize.c"

#define BASE_DOUBLE
#include "templates_on.h"
#include "c_init.c"
#include "c_main.c"
#include "c_pass_2.c"
#include "c_pass_3.c"
#include "c_pass_4.c"
#include "c_pass_5.c"
#include "c_pass_6.c"
#include "c_pass_7.c"
#include "c_pass_n.c"
#include "c_radix2.c"
#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
#include "c_init.c"
#include "c_main.c"
#include "c_pass_2.c"
#include "c_pass_3.c"
#include "c_pass_4.c"
#include "c_pass_5.c"
#include "c_pass_6.c"
#include "c_pass_7.c"
#include "c_pass_n.c"
#include "c_radix2.c"
#include "templates_off.h"
#undef  BASE_FLOAT

#include <gsl/gsl_fft_halfcomplex.h>
#include <gsl/gsl_fft_halfcomplex_float.h>

#define BASE_DOUBLE
#include "templates_on.h"
#include "hc_init.c"
#include "hc_main.c"
#include "hc_pass_2.c"
#include "hc_pass_3.c"
#include "hc_pass_4.c"
#include "hc_pass_5.c"
#include "hc_pass_n.c"
#include "hc_radix2.c"
#include "hc_unpack.c"
#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
#include "hc_init.c"
#include "hc_main.c"
#include "hc_pass_2.c"
#include "hc_pass_3.c"
#include "hc_pass_4.c"
#include "hc_pass_5.c"
#include "hc_pass_n.c"
#include "hc_radix2.c"
#include "hc_unpack.c"
#include "templates_off.h"
#undef  BASE_FLOAT

#include <gsl/gsl_fft_real.h>
#include <gsl/gsl_fft_real_float.h>

#define BASE_DOUBLE
#include "templates_on.h"
#include "real_init.c"
#include "real_main.c"
#include "real_pass_2.c"
#include "real_pass_3.c"
#include "real_pass_4.c"
#include "real_pass_5.c"
#include "real_pass_n.c"
#include "real_radix2.c"
#include "real_unpack.c"
#include "templates_off.h"
#undef  BASE_DOUBLE

#define BASE_FLOAT
#include "templates_on.h"
#include "real_init.c"
#include "real_main.c"
#include "real_pass_2.c"
#include "real_pass_3.c"
#include "real_pass_4.c"
#include "real_pass_5.c"
#include "real_pass_n.c"
#include "real_radix2.c"
#include "real_unpack.c"
#include "templates_off.h"
#undef  BASE_FLOAT
