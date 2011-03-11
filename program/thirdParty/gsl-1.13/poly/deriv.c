/* poly/eval.c
 * 
 * Copyright (C) 2009 Marc JOURDAIN
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
#include <gsl/gsl_errno.h>
#include <gsl/gsl_poly.h>

int
gsl_poly_eval_derivs (const double c[], const size_t lenc, const double x,
                      double res[], const size_t lenres)
{
  size_t i, n, nmax;
  size_t k, l, lmax;

  for (i = 0, n = 0, nmax = 0; i < lenres; i++)
    {
      if (n < lenc)
	{
	  res[i] = c[lenc - 1];
	  nmax = n;
	  n++;
	}
      else
	res[i] = 0.0;
    }

  for (i = 0; i < lenc - 1; i++)
    {
      k = (lenc - 1) - i;
      res[0] = ((x * res[0]) + c[k - 1]);
      lmax = (nmax < k) ? nmax : k - 1;
      for (l = 1; l <= lmax; l++)
	{
	  res[l] = ((x * res[l]) + res[l - 1]);
	}
    }

  {
    double f = 1.0;
    for (i = 2; i <= nmax; i++)
      {
        f *= i;
        res[i] *= f;
      }
  }

  return GSL_SUCCESS;
}
