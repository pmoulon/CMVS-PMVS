/*----------------------------------------------------------------------------*/
/*                                                                            */
/*  quad_golden.c                                                             */
/*                                                                            */
/*  Copyright (C) 2007 James Howse                                            */
/*  Copyright (C) 2009 Brian Gough                                             */
/*                                                                            */
/*  This program is free software; you can redistribute it and/or modify      */
/*  it under the terms of the GNU General Public License as published by      */
/*  the Free Software Foundation; either version 3 of the License, or (at     */
/*  your option) any later version.                                           */
/*                                                                            */
/*  This program is distributed in the hope that it will be useful, but       */
/*  WITHOUT ANY WARRANTY; without even the implied warranty of                */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU         */
/*  General Public License for more details.                                  */
/*                                                                            */
/*  You should have received a copy of the GNU General Public License         */
/*  along with this program; if not, write to the Free Software               */
/*  Foundation, Inc., 51 Franklin Street, Fifth Floor,                        */
/*  Boston, MA 02110-1301, USA.                                               */
/*                                                                            */
/*  ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::  */
/*                                                                            */
/*  This algorithm performs univariate minimization (i.e., line search).      */
/*  It requires only objective function values g(x) to compute the minimum.   */
/*  The algorithm maintains an interval of uncertainty [a,b] and a point x    */
/*  in the interval [a,b] such that a < x < b, and g(a) > g(x) and            */
/*  g(x) < g(b).  The algorithm also maintains the three points with the      */
/*  smallest objective values x, v and w such that g(x) < g(v) < g(w).  The   */
/*  algorithm terminates when max( x - a, b - x ) <  2(r |x| + t) where r     */
/*  and t are small positive reals.  At a given iteration, the algorithm      */
/*  first fits a quadratic through the three points (x, g(x)), (v, g(v))      */
/*  and (w, g(w)) and computes the location of the minimum u of the           */
/*  resulting quadratic.  If u is in the interval [a,b] then g(u) is          */
/*  computed.  If u is not in the interval [a,b], and either v < x and        */
/*  w < x, or v > x and w > x (i.e., the quadratic is extrapolating), then    */
/*  a point u' is computed using a safeguarding procedure and g(u') is        */
/*  computed.  If u is not in the interval [a,b], and the quadratic is not    */
/*  extrapolating, then a point u'' is computed using approximate golden      */
/*  section and g(u'') is computed.  After evaluating g() at the              */
/*  appropriate new point, a, b, x, v, and w are updated and the next         */
/*  iteration is performed.  The algorithm is based on work presented in      */
/*  the following references.                                                 */
/*                                                                            */
/*  Algorithms for Minimization without derivatives                           */
/*  Richard Brent                                                             */
/*  Prentice-Hall Inc., Englewood Cliffs, NJ, 1973                            */
/*                                                                            */
/*  Safeguarded Steplength Algorithms for Optimization using Descent Methods  */
/*  Philip E. Gill and Walter Murray                                          */
/*  Division of Numerical Analysis and Computing                              */
/*  National Physical Laboratory, Teddington, United Kingdom                  */
/*  NPL Report NAC 37, August 1974                                            */
/*                                                                            */
/*----------------------------------------------------------------------------*/
#include <config.h>

#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <gsl/gsl_math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_min.h>

#include "min.h"

#define REL_ERR_VAL   1.0e-06
#define ABS_ERR_VAL   1.0e-10
#define GOLDEN_MEAN   0.3819660112501052	/* (3 - sqrt(5))/2 */
#define GOLDEN_RATIO  1.6180339887498950	/* (1 + sqrt(5))/2 */

#define DEBUG_PRINTF(x) /* do nothing */

typedef struct
{
  double step_size, stored_step, prev_stored_step;
  double x_prev_small, f_prev_small, x_small, f_small;
  unsigned int num_iter;
}
quad_golden_state_t;

static int
quad_golden_init (void *vstate, gsl_function * f, double x_minimum,
		  double f_minimum, double x_lower, double f_lower,
		  double x_upper, double f_upper)
{
  quad_golden_state_t *state = (quad_golden_state_t *) vstate;

  /* For the original behavior, the first value for x_minimum_minimum
     passed in by the user should be a golden section step but we
     don't enforce this here. */

  state->x_prev_small = x_minimum;
  state->x_small = x_minimum;

  state->f_prev_small = f_minimum;
  state->f_small = f_minimum;

  state->step_size = 0.0;
  state->stored_step = 0.0;
  state->prev_stored_step = 0.0;
  state->num_iter = 0;

  x_lower = 0 ;  /* avoid warnings about unused variables */
  x_upper = 0 ;
  f_lower = 0 ;
  f_upper = 0 ;
  f = 0;

  return GSL_SUCCESS;
}

static int
quad_golden_iterate (void *vstate, gsl_function * f, double *x_minimum,
		     double *f_minimum, double *x_lower, double *f_lower,
		     double *x_upper, double *f_upper)
{
  quad_golden_state_t *state = (quad_golden_state_t *) vstate;

  const double x_m = *x_minimum;
  const double f_m = *f_minimum;

  const double x_l = *x_lower;
  const double x_u = *x_upper;

  const double x_small = state->x_small;
  const double f_small = state->f_small;

  const double x_prev_small = state->x_prev_small;
  const double f_prev_small = state->f_prev_small;
  
  double stored_step = state->stored_step; /* update on exit */
  double prev_stored_step = state->prev_stored_step; /* update on exit */
  double step_size = state->step_size; /* update on exit */

  double quad_step_size = prev_stored_step;
  
  double x_trial;
  double x_eval, f_eval;

  double x_midpoint = 0.5 * (x_l + x_u);
  double tol = REL_ERR_VAL * fabs (x_m) + ABS_ERR_VAL; /* total error tolerance */

  if (fabs (stored_step) - tol > -2.0 * GSL_DBL_EPSILON)
    {
      /* Fit quadratic */
      double c3 = (x_m - x_small) * (f_m - f_prev_small);
      double c2 = (x_m - x_prev_small) * (f_m - f_small);
      double c1 = (x_m - x_prev_small) * c2 - (x_m - x_small) * c3;

      c2 = 2.0 * (c2 - c3);

      if (fabs (c2) > GSL_DBL_EPSILON)	/* if( c2 != 0 ) */
	{
	  if (c2 > 0.0)
	    c1 = -c1;

	  c2 = fabs (c2);

	  quad_step_size = c1 / c2;
	}
      else
	{
	  /* Handle case where c2 ~=~ 0  */
	  /* Insure that the line search will NOT take a quadratic
	     interpolation step in this iteration */
	  quad_step_size = stored_step;
	}

      prev_stored_step = stored_step;
      stored_step = step_size;
    }

  x_trial = x_m + quad_step_size;

  if (fabs (quad_step_size) < fabs (0.5 * prev_stored_step) && x_trial > x_l && x_trial < x_u)
    {
      /* Take quadratic interpolation step */
      step_size = quad_step_size;

      /* Do not evaluate function too close to x_l or x_u */
      if ((x_trial - x_l) < 2.0 * tol || (x_u - x_trial) < 2.0 * tol)
        {
          step_size = (x_midpoint >= x_m ? +1.0 : -1.0) * fabs(tol);
        }

      DEBUG_PRINTF(("quadratic step: %g\n", step_size));
    }
  else if ((x_small != x_prev_small && x_small < x_m && x_prev_small < x_m) ||
           (x_small != x_prev_small && x_small > x_m && x_prev_small > x_m))
    {
      /* Take safeguarded function comparison step */
      double outside_interval, inside_interval;

      if (x_small < x_m)
	{
	  outside_interval = x_l - x_m;
	  inside_interval = x_u - x_m;
	}
      else
	{
	  outside_interval = x_u - x_m;
	  inside_interval = x_l - x_m;
	}

      if (fabs (inside_interval) <= tol)
	{
          /* Swap inside and outside intervals */
          double tmp = outside_interval;
	  outside_interval = inside_interval;
	  inside_interval = tmp;
	}

      {
        double step = inside_interval;
        double scale_factor;

        if (fabs (outside_interval) < fabs (inside_interval))
          {
            scale_factor = 0.5 * sqrt (-outside_interval / inside_interval);
          }
        else
          {
            scale_factor = (5.0 / 11.0) * (0.1 - inside_interval / outside_interval);
          }

        state->stored_step = step;
        step_size = scale_factor * step;
      }

      DEBUG_PRINTF(("safeguard step: %g\n", step_size));
    }
  else
    {
      /* Take golden section step */
      double step;

      if (x_m < x_midpoint)
        {
          step = x_u - x_m;
        }
      else
        {
          step = x_l - x_m;
        }

      state->stored_step = step;
      step_size = GOLDEN_MEAN * step;

      DEBUG_PRINTF(("golden step: %g\n", step_size));
    }

  /* Do not evaluate function too close to x_minimum */
  if (fabs (step_size) > tol)
    {
      x_eval = x_m + step_size;
    }
  else
    {
      x_eval = x_m + (step_size >= 0 ? +1.0 : -1.0) * fabs(tol);
    }

  /* Evaluate function at the new point x_eval */
  SAFE_FUNC_CALL(f, x_eval, &f_eval);

  /* Update {x,f}_lower, {x,f}_upper, {x,f}_prev_small, {x,f}_small, and {x,f}_minimum */
  if (f_eval <= f_m)
    {
      if (x_eval < x_m)
	{
          *x_upper = x_m;
          *f_upper = f_m;
        }
      else
      	{
          *x_lower = x_m;
          *f_upper = f_m;
        }

      state->x_prev_small = x_small;
      state->f_prev_small = f_small;

      state->x_small = x_m;
      state->f_small = f_m;

      *x_minimum = x_eval;
      *f_minimum = f_eval;
    }
  else
    {
      if (x_eval < x_m)
        {
          *x_lower = x_eval;
          *f_lower = f_eval;
        }
      else
        {    
          *x_upper = x_eval;
          *f_upper = f_eval;
        }

      if (f_eval <= f_small || fabs (x_small - x_m) < 2.0 * GSL_DBL_EPSILON)
	{
	  state->x_prev_small = x_small;
	  state->f_prev_small = f_small;

	  state->x_small = x_eval;
	  state->f_small = f_eval;
	}
      else if (f_eval <= f_prev_small ||
	       fabs (x_prev_small - x_m) < 2.0 * GSL_DBL_EPSILON ||
	       fabs (x_prev_small - x_small) < 2.0 * GSL_DBL_EPSILON)
	{
	  state->x_prev_small = x_eval;
	  state->f_prev_small = f_eval;
	}
    }

  /* Update stored values for next iteration */

  state->stored_step = stored_step;
  state->prev_stored_step = prev_stored_step;
  state->step_size = step_size;
  state->num_iter++;

  DEBUG_PRINTF(("[%d] Final State: %g  %g  %g\n", state->num_iter, x_l, x_m, x_u));

  return GSL_SUCCESS;
}


static const gsl_min_fminimizer_type quad_golden_type = { "quad-golden",	/* name */
  sizeof (quad_golden_state_t),
  &quad_golden_init,
  &quad_golden_iterate
};

const gsl_min_fminimizer_type *gsl_min_fminimizer_quad_golden =
  &quad_golden_type;
