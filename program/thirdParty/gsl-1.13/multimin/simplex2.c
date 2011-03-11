/* multimin/simplex2.c
 * 
 * Copyright (C) 2007, 2008, 2009 Brian Gough
 * Copyright (C) 2002 Tuomo Keskitalo, Ivo Alxneit
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

/*
   - Originally written by Tuomo Keskitalo <tuomo.keskitalo@iki.fi>
   - Corrections to nmsimplex_iterate and other functions 
     by Ivo Alxneit <ivo.alxneit@psi.ch>
   - Additional help by Brian Gough <bjg@network-theory.co.uk>
   - Optimisations added by Brian Gough <bjg@network-theory.co.uk>
         + use BLAS for frequently-called functions
         + keep track of the center to avoid unnecessary computation
         + compute size as RMS value, allowing linear update on each step
           instead of recomputing from all N+1 vectors.
*/

/* The Simplex method of Nelder and Mead, also known as the polytope
   search alogorithm.  Ref: Nelder, J.A., Mead, R., Computer Journal 7
   (1965) pp. 308-313.

   This implementation uses n+1 corner points in the simplex.
*/

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_matrix_double.h>

typedef struct
{
  gsl_matrix *x1;		/* simplex corner points */
  gsl_vector *y1;		/* function value at corner points */
  gsl_vector *ws1;		/* workspace 1 for algorithm */
  gsl_vector *ws2;		/* workspace 2 for algorithm */
  gsl_vector *center;		/* center of all points */
  gsl_vector *delta;		/* current step */
  gsl_vector *xmc;		/* x - center (workspace) */
  double S2;
  unsigned long count;
}
nmsimplex_state_t;

static int
compute_center (const nmsimplex_state_t * state, gsl_vector * center);
static double
compute_size (nmsimplex_state_t * state, const gsl_vector * center);

static double
try_corner_move (const double coeff,
		 const nmsimplex_state_t * state,
		 size_t corner,
		 gsl_vector * xc, const gsl_multimin_function * f)
{
  /* moves a simplex corner scaled by coeff (negative value represents 
     mirroring by the middle point of the "other" corner points)
     and gives new corner in xc and function value at xc as a 
     return value 
   */

  gsl_matrix *x1 = state->x1;
  const size_t P = x1->size1;
  double newval;

  /* xc = (1-coeff)*(P/(P-1)) * center(all) + ((P*coeff-1)/(P-1))*x_corner */
  {
    double alpha = (1 - coeff) * P / (P - 1.0);
    double beta = (P * coeff - 1.0) / (P - 1.0);
    gsl_vector_const_view row = gsl_matrix_const_row (x1, corner);

    gsl_vector_memcpy (xc, state->center);
    gsl_blas_dscal (alpha, xc);
    gsl_blas_daxpy (beta, &row.vector, xc);
  }

  newval = GSL_MULTIMIN_FN_EVAL (f, xc);

  return newval;
}


static void
update_point (nmsimplex_state_t * state, size_t i,
	      const gsl_vector * x, double val)
{
  gsl_vector_const_view x_orig = gsl_matrix_const_row (state->x1, i);
  const size_t P = state->x1->size1;

  /* Compute delta = x - x_orig */
  gsl_vector_memcpy (state->delta, x);
  gsl_blas_daxpy (-1.0, &x_orig.vector, state->delta);

  /* Compute xmc = x_orig - c */
  gsl_vector_memcpy (state->xmc, &x_orig.vector);
  gsl_blas_daxpy (-1.0, state->center, state->xmc);

  /* Update size: S2' = S2 + (2/P) * (x_orig - c).delta + (P-1)*(delta/P)^2 */
  {
    double d = gsl_blas_dnrm2 (state->delta);
    double xmcd;
    gsl_blas_ddot (state->xmc, state->delta, &xmcd);
    state->S2 += (2.0 / P) * xmcd + ((P - 1.0) / P) * (d * d / P);
  }

  /* Update center:  c' = c + (x - x_orig) / P */

  {
    double alpha = 1.0 / P;
    gsl_blas_daxpy (-alpha, &x_orig.vector, state->center);
    gsl_blas_daxpy (alpha, x, state->center);
  }

  gsl_matrix_set_row (state->x1, i, x);
  gsl_vector_set (state->y1, i, val);
}

static int
contract_by_best (nmsimplex_state_t * state, size_t best,
		  gsl_vector * xc, gsl_multimin_function * f)
{

  /* Function contracts the simplex in respect to best valued
     corner. That is, all corners besides the best corner are moved.
     (This function is rarely called in practice, since it is the last
     choice, hence not optimised - BJG)  */

  /* the xc vector is simply work space here */

  gsl_matrix *x1 = state->x1;
  gsl_vector *y1 = state->y1;

  size_t i, j;
  double newval;

  int status = GSL_SUCCESS;

  for (i = 0; i < x1->size1; i++)
    {
      if (i != best)
	{
	  for (j = 0; j < x1->size2; j++)
	    {
	      newval = 0.5 * (gsl_matrix_get (x1, i, j)
			      + gsl_matrix_get (x1, best, j));
	      gsl_matrix_set (x1, i, j, newval);
	    }

	  /* evaluate function in the new point */

	  gsl_matrix_get_row (xc, x1, i);
	  newval = GSL_MULTIMIN_FN_EVAL (f, xc);
	  gsl_vector_set (y1, i, newval);

	  /* notify caller that we found at least one bad function value.
	     we finish the contraction (and do not abort) to allow the user
	     to handle the situation */

	  if (!gsl_finite (newval))
	    {
	      status = GSL_EBADFUNC;
	    }
	}
    }

  /* We need to update the centre and size as well */
  compute_center (state, state->center);
  compute_size (state, state->center);

  return status;
}

static int
compute_center (const nmsimplex_state_t * state, gsl_vector * center)
{
  /* calculates the center of the simplex and stores in center */

  gsl_matrix *x1 = state->x1;
  const size_t P = x1->size1;
  size_t i;

  gsl_vector_set_zero (center);

  for (i = 0; i < P; i++)
    {
      gsl_vector_const_view row = gsl_matrix_const_row (x1, i);
      gsl_blas_daxpy (1.0, &row.vector, center);
    }

  {
    const double alpha = 1.0 / P;
    gsl_blas_dscal (alpha, center);
  }

  return GSL_SUCCESS;
}

static double
compute_size (nmsimplex_state_t * state, const gsl_vector * center)
{
  /* calculates simplex size as rms sum of length of vectors 
     from simplex center to corner points:     

     sqrt( sum ( || y - y_middlepoint ||^2 ) / n )
   */

  gsl_vector *s = state->ws1;
  gsl_matrix *x1 = state->x1;
  const size_t P = x1->size1;
  size_t i;

  double ss = 0.0;

  for (i = 0; i < P; i++)
    {
      double t;
      gsl_matrix_get_row (s, x1, i);
      gsl_blas_daxpy (-1.0, center, s);
      t = gsl_blas_dnrm2 (s);
      ss += t * t;
    }

  /* Store squared size in the state */
  state->S2 = (ss / P);

  return sqrt (ss / P);
}

static int
nmsimplex_alloc (void *vstate, size_t n)
{
  nmsimplex_state_t *state = (nmsimplex_state_t *) vstate;

  if (n == 0)
    {
      GSL_ERROR ("invalid number of parameters specified", GSL_EINVAL);
    }

  state->x1 = gsl_matrix_alloc (n + 1, n);

  if (state->x1 == NULL)
    {
      GSL_ERROR ("failed to allocate space for x1", GSL_ENOMEM);
    }

  state->y1 = gsl_vector_alloc (n + 1);

  if (state->y1 == NULL)
    {
      gsl_matrix_free (state->x1);
      GSL_ERROR ("failed to allocate space for y", GSL_ENOMEM);
    }

  state->ws1 = gsl_vector_alloc (n);

  if (state->ws1 == NULL)
    {
      gsl_matrix_free (state->x1);
      gsl_vector_free (state->y1);
      GSL_ERROR ("failed to allocate space for ws1", GSL_ENOMEM);
    }

  state->ws2 = gsl_vector_alloc (n);

  if (state->ws2 == NULL)
    {
      gsl_matrix_free (state->x1);
      gsl_vector_free (state->y1);
      gsl_vector_free (state->ws1);
      GSL_ERROR ("failed to allocate space for ws2", GSL_ENOMEM);
    }

  state->center = gsl_vector_alloc (n);

  if (state->center == NULL)
    {
      gsl_matrix_free (state->x1);
      gsl_vector_free (state->y1);
      gsl_vector_free (state->ws1);
      gsl_vector_free (state->ws2);
      GSL_ERROR ("failed to allocate space for center", GSL_ENOMEM);
    }

  state->delta = gsl_vector_alloc (n);

  if (state->delta == NULL)
    {
      gsl_matrix_free (state->x1);
      gsl_vector_free (state->y1);
      gsl_vector_free (state->ws1);
      gsl_vector_free (state->ws2);
      gsl_vector_free (state->center);
      GSL_ERROR ("failed to allocate space for delta", GSL_ENOMEM);
    }

  state->xmc = gsl_vector_alloc (n);

  if (state->xmc == NULL)
    {
      gsl_matrix_free (state->x1);
      gsl_vector_free (state->y1);
      gsl_vector_free (state->ws1);
      gsl_vector_free (state->ws2);
      gsl_vector_free (state->center);
      gsl_vector_free (state->delta);
      GSL_ERROR ("failed to allocate space for xmc", GSL_ENOMEM);
    }

  state->count = 0;

  return GSL_SUCCESS;
}

static void
nmsimplex_free (void *vstate)
{
  nmsimplex_state_t *state = (nmsimplex_state_t *) vstate;

  gsl_matrix_free (state->x1);
  gsl_vector_free (state->y1);
  gsl_vector_free (state->ws1);
  gsl_vector_free (state->ws2);
  gsl_vector_free (state->center);
  gsl_vector_free (state->delta);
  gsl_vector_free (state->xmc);
}

static int
nmsimplex_set (void *vstate, gsl_multimin_function * f,
	       const gsl_vector * x,
	       double *size, const gsl_vector * step_size)
{
  int status;
  size_t i;
  double val;

  nmsimplex_state_t *state = (nmsimplex_state_t *) vstate;

  gsl_vector *xtemp = state->ws1;

  if (xtemp->size != x->size)
    {
      GSL_ERROR ("incompatible size of x", GSL_EINVAL);
    }

  if (xtemp->size != step_size->size)
    {
      GSL_ERROR ("incompatible size of step_size", GSL_EINVAL);
    }

  /* first point is the original x0 */

  val = GSL_MULTIMIN_FN_EVAL (f, x);

  if (!gsl_finite (val))
    {
      GSL_ERROR ("non-finite function value encountered", GSL_EBADFUNC);
    }

  gsl_matrix_set_row (state->x1, 0, x);
  gsl_vector_set (state->y1, 0, val);

  /* following points are initialized to x0 + step_size */

  for (i = 0; i < x->size; i++)
    {
      status = gsl_vector_memcpy (xtemp, x);

      if (status != 0)
	{
	  GSL_ERROR ("vector memcopy failed", GSL_EFAILED);
	}

      {
	double xi = gsl_vector_get (x, i);
	double si = gsl_vector_get (step_size, i);

	gsl_vector_set (xtemp, i, xi + si);
	val = GSL_MULTIMIN_FN_EVAL (f, xtemp);
      }

      if (!gsl_finite (val))
	{
	  GSL_ERROR ("non-finite function value encountered", GSL_EBADFUNC);
	}

      gsl_matrix_set_row (state->x1, i + 1, xtemp);
      gsl_vector_set (state->y1, i + 1, val);
    }

  compute_center (state, state->center);

  /* Initialize simplex size */
  *size = compute_size (state, state->center);

  state->count++;

  return GSL_SUCCESS;
}

static int
nmsimplex_iterate (void *vstate, gsl_multimin_function * f,
		   gsl_vector * x, double *size, double *fval)
{

  /* Simplex iteration tries to minimize function f value */
  /* Includes corrections from Ivo Alxneit <ivo.alxneit@psi.ch> */

  nmsimplex_state_t *state = (nmsimplex_state_t *) vstate;

  /* xc and xc2 vectors store tried corner point coordinates */

  gsl_vector *xc = state->ws1;
  gsl_vector *xc2 = state->ws2;
  gsl_vector *y1 = state->y1;
  gsl_matrix *x1 = state->x1;

  const size_t n = y1->size;
  size_t i;
  size_t hi, s_hi, lo;
  double dhi, ds_hi, dlo;
  int status;
  double val, val2;

  if (xc->size != x->size)
    {
      GSL_ERROR ("incompatible size of x", GSL_EINVAL);
    }

  /* get index of highest, second highest and lowest point */

  dhi = dlo = gsl_vector_get (y1, 0);
  hi = 0;
  lo = 0;

  ds_hi = gsl_vector_get (y1, 1);
  s_hi = 1;

  for (i = 1; i < n; i++)
    {
      val = (gsl_vector_get (y1, i));
      if (val < dlo)
	{
	  dlo = val;
	  lo = i;
	}
      else if (val > dhi)
	{
	  ds_hi = dhi;
	  s_hi = hi;
	  dhi = val;
	  hi = i;
	}
      else if (val > ds_hi)
	{
	  ds_hi = val;
	  s_hi = i;
	}
    }

  /* try reflecting the highest value point */

  val = try_corner_move (-1.0, state, hi, xc, f);

  if (gsl_finite (val) && val < gsl_vector_get (y1, lo))
    {
      /* reflected point is lowest, try expansion */

      val2 = try_corner_move (-2.0, state, hi, xc2, f);

      if (gsl_finite (val2) && val2 < gsl_vector_get (y1, lo))
	{
	  update_point (state, hi, xc2, val2);
	}
      else
	{
	  update_point (state, hi, xc, val);
	}
    }
  else if (!gsl_finite (val) || val > gsl_vector_get (y1, s_hi))
    {
      /* reflection does not improve things enough, or we got a
         non-finite function value */

      if (gsl_finite (val) && val <= gsl_vector_get (y1, hi))
	{
	  /* if trial point is better than highest point, replace
	     highest point */

	  update_point (state, hi, xc, val);
	}

      /* try one-dimensional contraction */

      val2 = try_corner_move (0.5, state, hi, xc2, f);

      if (gsl_finite (val2) && val2 <= gsl_vector_get (y1, hi))
	{
	  update_point (state, hi, xc2, val2);
	}
      else
	{
	  /* contract the whole simplex about the best point */

	  status = contract_by_best (state, lo, xc, f);

	  if (status != GSL_SUCCESS)
	    {
	      GSL_ERROR ("contraction failed", GSL_EFAILED);
	    }
	}
    }
  else
    {
      /* trial point is better than second highest point.  Replace
         highest point by it */

      update_point (state, hi, xc, val);
    }

  /* return lowest point of simplex as x */

  lo = gsl_vector_min_index (y1);
  gsl_matrix_get_row (x, x1, lo);
  *fval = gsl_vector_get (y1, lo);

  /* Update simplex size */
  
  {
    double S2 = state->S2;

    if (S2 > 0)
      {
	*size = sqrt (S2);
      }
    else
      {
	/* recompute if accumulated error has made size invalid */
	*size = compute_size (state, state->center);
      }
  }

  return GSL_SUCCESS;
}

static const gsl_multimin_fminimizer_type nmsimplex_type = 
{ "nmsimplex2",	/* name */
  sizeof (nmsimplex_state_t),
  &nmsimplex_alloc,
  &nmsimplex_set,
  &nmsimplex_iterate,
  &nmsimplex_free
};

const gsl_multimin_fminimizer_type
  * gsl_multimin_fminimizer_nmsimplex2 = &nmsimplex_type;


static inline double
ran_unif (unsigned long *seed)
{
  unsigned long s = *seed;
  *seed = (s * 69069 + 1) & 0xffffffffUL;
  return (*seed) / 4294967296.0;
}

static int
nmsimplex_set_rand (void *vstate, gsl_multimin_function * f,
		    const gsl_vector * x,
		    double *size, const gsl_vector * step_size)
{
  size_t i, j;
  double val;

  nmsimplex_state_t *state = (nmsimplex_state_t *) vstate;

  gsl_vector *xtemp = state->ws1;

  if (xtemp->size != x->size)
    {
      GSL_ERROR ("incompatible size of x", GSL_EINVAL);
    }

  if (xtemp->size != step_size->size)
    {
      GSL_ERROR ("incompatible size of step_size", GSL_EINVAL);
    }

  /* first point is the original x0 */

  val = GSL_MULTIMIN_FN_EVAL (f, x);

  if (!gsl_finite (val))
    {
      GSL_ERROR ("non-finite function value encountered", GSL_EBADFUNC);
    }

  gsl_matrix_set_row (state->x1, 0, x);
  gsl_vector_set (state->y1, 0, val);

  {
    gsl_matrix_view m =
      gsl_matrix_submatrix (state->x1, 1, 0, x->size, x->size);

    /* generate a random orthornomal basis  */
    unsigned long seed = state->count ^ 0x12345678;

    ran_unif (&seed);		/* warm it up */

    gsl_matrix_set_identity (&m.matrix);

    /* start with random reflections */
    for (i = 0; i < x->size; i++)
      {
        double s = ran_unif (&seed);
        if (s > 0.5) gsl_matrix_set (&m.matrix, i, i, -1.0);
      }

    /* apply random rotations */
    for (i = 0; i < x->size; i++)
      {
	for (j = i + 1; j < x->size; j++)
	  {
	    /* rotate columns i and j by a random angle */
	    double angle = 2.0 * M_PI * ran_unif (&seed);
	    double c = cos (angle), s = sin (angle);
	    gsl_vector_view c_i = gsl_matrix_column (&m.matrix, i);
	    gsl_vector_view c_j = gsl_matrix_column (&m.matrix, j);
	    gsl_blas_drot (&c_i.vector, &c_j.vector, c, s);
	  }
      }

    /* scale the orthonormal basis by the user-supplied step_size in
       each dimension, and use as an offset from the central point x */

    for (i = 0; i < x->size; i++)
      {
	double x_i = gsl_vector_get (x, i);
	double s_i = gsl_vector_get (step_size, i);
	gsl_vector_view c_i = gsl_matrix_column (&m.matrix, i);

	for (j = 0; j < x->size; j++)
	  {
	    double x_ij = gsl_vector_get (&c_i.vector, j);
	    gsl_vector_set (&c_i.vector, j, x_i + s_i * x_ij);
	  }
      }

    /* compute the function values at each offset point */

    for (i = 0; i < x->size; i++)
      {
	gsl_vector_view r_i = gsl_matrix_row (&m.matrix, i);

	val = GSL_MULTIMIN_FN_EVAL (f, &r_i.vector);

	if (!gsl_finite (val))
	  {
	    GSL_ERROR ("non-finite function value encountered", GSL_EBADFUNC);
	  }

	gsl_vector_set (state->y1, i + 1, val);
      }
  }

  compute_center (state, state->center);

  /* Initialize simplex size */
  *size = compute_size (state, state->center);

  state->count++;

  return GSL_SUCCESS;
}

static const gsl_multimin_fminimizer_type nmsimplex2rand_type = 
{ "nmsimplex2rand",	/* name */
  sizeof (nmsimplex_state_t),
  &nmsimplex_alloc,
  &nmsimplex_set_rand,
  &nmsimplex_iterate,
  &nmsimplex_free
};

const gsl_multimin_fminimizer_type
  * gsl_multimin_fminimizer_nmsimplex2rand = &nmsimplex2rand_type;
