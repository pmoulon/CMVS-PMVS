/* cheb/deriv.c
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

int gsl_cheb_calc_deriv(gsl_cheb_series * deriv, const gsl_cheb_series * f)
{
  const size_t n = f->order + 1;
  const double con = 2.0 / (f->b - f->a);
  size_t i;
  
  if(deriv->order != f->order) 
    {
      GSL_ERROR ("order of chebyshev series must be equal", GSL_ENOMEM);
    }
  
  /* set the other parameters in the chebyshev struct */

  deriv->a = f->a;
  deriv->b = f->b;

  /* error in derivative is n^2 c_n */ 
  /* deriv->err = n * n * f->c[n-1];*/   

  /* FIXME:  should probably set deriv->f[] as well */
  
  deriv->c[n-1] = 0.0;
  
  if(n > 1) {
    deriv->c[n-2] = 2.0 *(n-1.0) * f->c[n-1];

    for(i = n-3; i>0; i--) 
      deriv->c[i] = deriv->c[i+2] + 2.0 *(i+1.0) * f->c[i+1];

    deriv->c[0] = deriv->c[2] + 2.0 * f->c[1];

    for(i = 0  ; i<n ; i++) 
      deriv->c[i] *= con;
  }

  return GSL_SUCCESS;
}
