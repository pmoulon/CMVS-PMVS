/* cheb/init.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 * Copyright (C) 2009 Brian Gough
 * 
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_chebyshev.h>

/*-*-*-*-*-*-*-*-*-*-*-* Allocators *-*-*-*-*-*-*-*-*-*-*-*/

gsl_cheb_series * 
gsl_cheb_alloc(const size_t order)
{
  gsl_cheb_series * cs = (gsl_cheb_series *) malloc(sizeof(gsl_cheb_series));
  
  if(cs == 0) {
    GSL_ERROR_VAL("failed to allocate gsl_cheb_series struct", GSL_ENOMEM, 0);
  }
  
  cs->order    = order;
  cs->order_sp = order;

  cs->c = (double *) malloc((order+1) * sizeof(double));

  if(cs->c == 0) {
    GSL_ERROR_VAL("failed to allocate cheb coefficients", GSL_ENOMEM, 0);
  }

  cs->f = (double *) malloc((order+1) * sizeof(double));

  if(cs->f == 0) {
    GSL_ERROR_VAL("failed to allocate cheb function space", GSL_ENOMEM, 0);
  }

  return cs;
}


void gsl_cheb_free(gsl_cheb_series * cs)
{
  RETURN_IF_NULL (cs);
  free(cs->f);
  free(cs->c);
  free(cs);
}

/*-*-*-*-*-*-*-*-*-*-*-* Initializer *-*-*-*-*-*-*-*-*-*-*-*/

int gsl_cheb_init(gsl_cheb_series * cs, const gsl_function *func,
                  const double a, const double b)
{
  size_t k, j;

  if(a >= b) {
    GSL_ERROR_VAL("null function interval [a,b]", GSL_EDOM, 0);
  }
  cs->a = a;
  cs->b = b;
  /* cs->err = 0.0; */

  { 
    double bma = 0.5 * (cs->b - cs->a);
    double bpa = 0.5 * (cs->b + cs->a);
    double fac = 2.0/(cs->order +1.0);

    for(k = 0; k<=cs->order; k++) {
      double y = cos(M_PI * (k+0.5)/(cs->order+1));
      cs->f[k] = GSL_FN_EVAL(func, (y*bma + bpa));
    }
    
    for(j = 0; j<=cs->order; j++) {
      double sum = 0.0;
      for(k = 0; k<=cs->order; k++) 
        sum += cs->f[k]*cos(M_PI * j*(k+0.5)/(cs->order+1));
      cs->c[j] = fac * sum;
    }
    
  }
  return GSL_SUCCESS;
}

size_t
gsl_cheb_order (const gsl_cheb_series * cs)
{
  return cs->order;
}

size_t
gsl_cheb_size (const gsl_cheb_series * cs)
{
  return (cs->order + 1);
}

double *
gsl_cheb_coeffs (const gsl_cheb_series * cs)
{
  return cs->c;
}

