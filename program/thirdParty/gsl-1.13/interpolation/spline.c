/* interpolation/spline.c
 * 
 * Copyright (C) 2001, 2007 Brian Gough
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
#include <string.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

gsl_spline *
gsl_spline_alloc (const gsl_interp_type * T, size_t size)
{
  gsl_spline * spline = (gsl_spline *) malloc (sizeof(gsl_spline));
  
  if (spline == NULL)
    {
      GSL_ERROR_NULL ("failed to allocate space for spline struct", 
                      GSL_ENOMEM);
    }
  
  spline->interp = gsl_interp_alloc (T, size);
  
  if (spline->interp == NULL)
    {
      free (spline);          
      GSL_ERROR_NULL ("failed to allocate space for interp", GSL_ENOMEM);
    };
    
  spline->x = (double *) malloc (size * sizeof(double));

  if (spline->x == NULL)
    {
      gsl_interp_free(spline->interp);
      free(spline);
      GSL_ERROR_NULL ("failed to allocate space for x", GSL_ENOMEM);
    }

  spline->y = (double *) malloc (size * sizeof(double));

  if (spline->y == NULL)
    {
      free(spline->x);
      gsl_interp_free(spline->interp);
      free(spline);
      GSL_ERROR_NULL ("failed to allocate space for y", GSL_ENOMEM);
    }
  
  spline->size = size;

  return spline;
}

int
gsl_spline_init (gsl_spline * spline, const double x_array[], const double y_array[], size_t size)
{
  if (size != spline->size)
    {
      GSL_ERROR ("data must match size of spline object", GSL_EINVAL);
    }
  
  memcpy (spline->x, x_array, size * sizeof(double));
  memcpy (spline->y, y_array, size * sizeof(double));

  {
    int status = gsl_interp_init (spline->interp, x_array, y_array, size);
    return status;
  }
}

const char *
gsl_spline_name(const gsl_spline * spline)
{
  return gsl_interp_name(spline->interp);
}

unsigned int
gsl_spline_min_size(const gsl_spline * spline)
{
  return gsl_interp_min_size(spline->interp);
}

void
gsl_spline_free (gsl_spline * spline)
{
  RETURN_IF_NULL (spline);
  gsl_interp_free (spline->interp);
  free (spline->x);
  free (spline->y);
  free (spline);
}

int
gsl_spline_eval_e (const gsl_spline * spline, 
                   double x,
                   gsl_interp_accel * a, double *y)
{
  return gsl_interp_eval_e (spline->interp, 
                            spline->x, spline->y,
                            x, a, y);
}

double
gsl_spline_eval (const gsl_spline * spline,
                 double x,
                 gsl_interp_accel * a)
{
  return gsl_interp_eval (spline->interp, 
                          spline->x, spline->y,
                          x, a);
}


int
gsl_spline_eval_deriv_e (const gsl_spline * spline,
                         double x,
                         gsl_interp_accel * a,
                         double *dydx)
{
  return gsl_interp_eval_deriv_e (spline->interp, 
                                  spline->x, spline->y,
                                  x, a, dydx);
}

double
gsl_spline_eval_deriv (const gsl_spline * spline,
                       double x,
                       gsl_interp_accel * a)
{
  return gsl_interp_eval_deriv (spline->interp, 
                                spline->x, spline->y,
                                x, a);
}


int
gsl_spline_eval_deriv2_e (const gsl_spline * spline,
                          double x,
                          gsl_interp_accel * a,
                          double * d2)
{
  return gsl_interp_eval_deriv2_e (spline->interp, 
                                   spline->x, spline->y,
                                   x, a, d2);
}

double
gsl_spline_eval_deriv2 (const gsl_spline * spline,
                        double x,
                        gsl_interp_accel * a)
{
  return gsl_interp_eval_deriv2 (spline->interp, 
                                 spline->x, spline->y,
                                 x, a);
}


int
gsl_spline_eval_integ_e (const gsl_spline * spline,
                         double a, double b,
                         gsl_interp_accel * acc,
                         double * result)
{
  return gsl_interp_eval_integ_e (spline->interp, 
                                  spline->x, spline->y,
                                  a, b, acc, result);
}


double
gsl_spline_eval_integ (const gsl_spline * spline,
                       double a, double b,
                       gsl_interp_accel * acc)
{
  return gsl_interp_eval_integ (spline->interp, 
                                spline->x, spline->y,
                                a, b, acc);
}

