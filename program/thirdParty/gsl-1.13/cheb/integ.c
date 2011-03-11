/* cheb/integ.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
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

int gsl_cheb_calc_integ(gsl_cheb_series * integ, const gsl_cheb_series * f)
{
  const size_t n = f->order + 1;
  const double con = 0.25 * (f->b - f->a);

  if(integ->order != f->order) 
    {
      GSL_ERROR ("order of chebyshev series must be equal", GSL_ENOMEM);
    }

  /* set the other parameters in the chebyshev struct */

  integ->a = f->a;
  integ->b = f->b;

  /* FIXME:  should probably set integ->f[] as well */

  if(n == 1) {
    integ->c[0] = 0.;
  }
  else if(n == 2) {
    integ->c[1] = con * f->c[0];
    integ->c[0] = 2.0 * integ->c[1];
  }
  else {
    double sum = 0.0;
    double fac = 1.0;
    size_t i;
    for(i=1; i<=n-2; i++) {
      integ->c[i] = con * (f->c[i-1] - f->c[i+1])/((double)i);
      sum += fac * integ->c[i];
      fac = -fac;
    }
    integ->c[n-1] = con * f->c[n-2]/(n-1.0);
    sum += fac * integ->c[n-1];
    integ->c[0] = 2.0 * sum;
  }

  return GSL_SUCCESS;
}
