/* bspline/bspline.c
 *
 * Copyright (C) 2006, 2007, 2008, 2009 Patrick Alken
 * Copyright (C) 2008 Rhys Ulerich
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
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_statistics.h>

/*
 * This module contains routines related to calculating B-splines.
 * The algorithms used are described in
 *
 * [1] Carl de Boor, "A Practical Guide to Splines", Springer
 *     Verlag, 1978.
 *
 * The bspline_pppack_* internal routines contain code adapted from
 *
 * [2] "PPPACK - Piecewise Polynomial Package",
 *     http://www.netlib.org/pppack/
 *
 */

#include "bspline.h"

/*
gsl_bspline_alloc()
  Allocate space for a bspline workspace. The size of the
workspace is O(5k + nbreak)

Inputs: k      - spline order (cubic = 4)
        nbreak - number of breakpoints

Return: pointer to workspace
*/

gsl_bspline_workspace *
gsl_bspline_alloc (const size_t k, const size_t nbreak)
{
  if (k == 0)
    {
      GSL_ERROR_NULL ("k must be at least 1", GSL_EINVAL);
    }
  else if (nbreak < 2)
    {
      GSL_ERROR_NULL ("nbreak must be at least 2", GSL_EINVAL);
    }
  else
    {
      gsl_bspline_workspace *w;

      w = (gsl_bspline_workspace *) malloc (sizeof (gsl_bspline_workspace));

      if (w == 0)
	{
	  GSL_ERROR_NULL ("failed to allocate space for workspace",
			  GSL_ENOMEM);
	}

      w->k = k;
      w->km1 = k - 1;
      w->nbreak = nbreak;
      w->l = nbreak - 1;
      w->n = w->l + k - 1;

      w->knots = gsl_vector_alloc (w->n + k);
      if (w->knots == 0)
	{
	  free (w);
	  GSL_ERROR_NULL ("failed to allocate space for knots vector",
			  GSL_ENOMEM);
	}

      w->deltal = gsl_vector_alloc (k);
      if (w->deltal == 0)
	{
	  gsl_vector_free (w->knots);
	  free (w);
	  GSL_ERROR_NULL ("failed to allocate space for deltal vector",
			  GSL_ENOMEM);
	}

      w->deltar = gsl_vector_alloc (k);
      if (w->deltar == 0)
	{
	  gsl_vector_free (w->deltal);
	  gsl_vector_free (w->knots);
	  free (w);
	  GSL_ERROR_NULL ("failed to allocate space for deltar vector",
			  GSL_ENOMEM);
	}


      w->B = gsl_vector_alloc (k);
      if (w->B == 0)
	{
	  gsl_vector_free (w->deltar);;
	  gsl_vector_free (w->deltal);
	  gsl_vector_free (w->knots);
	  free (w);
	  GSL_ERROR_NULL
	    ("failed to allocate space for temporary spline vector",
	     GSL_ENOMEM);
	}

      return w;
    }
}				/* gsl_bspline_alloc() */

/*
gsl_bspline_deriv_alloc()
  Allocate space for a bspline derivative workspace. The size of the
workspace is O(2k^2)

Inputs: k      - spline order (cubic = 4)

Return: pointer to workspace
*/

gsl_bspline_deriv_workspace *
gsl_bspline_deriv_alloc (const size_t k)
{
  if (k == 0)
    {
      GSL_ERROR_NULL ("k must be at least 1", GSL_EINVAL);
    }
  else
    {
      gsl_bspline_deriv_workspace *dw;

      dw =
	(gsl_bspline_deriv_workspace *)
	malloc (sizeof (gsl_bspline_deriv_workspace));

      if (dw == 0)
	{
	  GSL_ERROR_NULL ("failed to allocate space for workspace",
			  GSL_ENOMEM);
	}

      dw->A = gsl_matrix_alloc (k, k);
      if (dw->A == 0)
	{
	  free (dw);
	  GSL_ERROR_NULL
	    ("failed to allocate space for derivative work matrix",
	     GSL_ENOMEM);
	}

      dw->dB = gsl_matrix_alloc (k, k + 1);
      if (dw->dB == 0)
	{
	  gsl_matrix_free (dw->A);
	  free (dw);
	  GSL_ERROR_NULL
	    ("failed to allocate space for temporary derivative matrix",
	     GSL_ENOMEM);
	}

      dw->k = k;

      return dw;
    }
}				/* gsl_bspline_deriv_alloc() */

/* Return number of coefficients */
size_t
gsl_bspline_ncoeffs (gsl_bspline_workspace * w)
{
  return w->n;
}

/* Return order */
size_t
gsl_bspline_order (gsl_bspline_workspace * w)
{
  return w->k;
}

/* Return number of breakpoints */
size_t
gsl_bspline_nbreak (gsl_bspline_workspace * w)
{
  return w->nbreak;
}

/* Return the location of the i-th breakpoint*/
double
gsl_bspline_breakpoint (size_t i, gsl_bspline_workspace * w)
{
  size_t j = i + w->k - 1;
  return gsl_vector_get (w->knots, j);
}

/* Return the location of the i-th Greville abscissa */
double
gsl_bspline_greville_abscissa(size_t i, gsl_bspline_workspace *w)
{
#if GSL_RANGE_CHECK
  if (GSL_RANGE_COND(i >= gsl_bspline_ncoeffs(w)))
    {
      GSL_ERROR_VAL ("Greville abscissa index out of range", GSL_EINVAL, 0);
    }
#endif
  const size_t stride = w->knots->stride;
  size_t km1 = w->km1;
  double * data = w->knots->data + (i+1)*stride;

  if (km1 == 0)
    {
      /* Return interval midpoints in degenerate k = 1 case*/
      km1   = 2;
      data -= stride;
    }

  return gsl_stats_mean(data, stride, km1);
}

/*
gsl_bspline_free()
  Free a gsl_bspline_workspace.

Inputs: w - workspace to free

Return: none
*/

void
gsl_bspline_free (gsl_bspline_workspace * w)
{
  RETURN_IF_NULL (w);
  gsl_vector_free (w->knots);
  gsl_vector_free (w->deltal);
  gsl_vector_free (w->deltar);
  gsl_vector_free (w->B);
  free (w);
}				/* gsl_bspline_free() */

/*
gsl_bspline_deriv_free()
  Free a gsl_bspline_deriv_workspace.

Inputs: dw - workspace to free

Return: none
*/

void
gsl_bspline_deriv_free (gsl_bspline_deriv_workspace * dw)
{
  RETURN_IF_NULL (dw);
  gsl_matrix_free (dw->A);
  gsl_matrix_free (dw->dB);
  free (dw);
}				/* gsl_bspline_deriv_free() */

/*
gsl_bspline_knots()
  Compute the knots from the given breakpoints:

   knots(1:k) = breakpts(1)
   knots(k+1:k+l-1) = breakpts(i), i = 2 .. l
   knots(n+1:n+k) = breakpts(l + 1)

where l is the number of polynomial pieces (l = nbreak - 1) and
   n = k + l - 1
(using matlab syntax for the arrays)

The repeated knots at the beginning and end of the interval
correspond to the continuity condition there. See pg. 119
of [1].

Inputs: breakpts - breakpoints
        w        - bspline workspace

Return: success or error
*/

int
gsl_bspline_knots (const gsl_vector * breakpts, gsl_bspline_workspace * w)
{
  if (breakpts->size != w->nbreak)
    {
      GSL_ERROR ("breakpts vector has wrong size", GSL_EBADLEN);
    }
  else
    {
      size_t i;			/* looping */

      for (i = 0; i < w->k; i++)
	gsl_vector_set (w->knots, i, gsl_vector_get (breakpts, 0));

      for (i = 1; i < w->l; i++)
	{
	  gsl_vector_set (w->knots, w->k - 1 + i,
			  gsl_vector_get (breakpts, i));
	}

      for (i = w->n; i < w->n + w->k; i++)
	gsl_vector_set (w->knots, i, gsl_vector_get (breakpts, w->l));

      return GSL_SUCCESS;
    }
}				/* gsl_bspline_knots() */

/*
gsl_bspline_knots_uniform()
  Construct uniformly spaced knots on the interval [a,b] using
the previously specified number of breakpoints. 'a' is the position
of the first breakpoint and 'b' is the position of the last
breakpoint.

Inputs: a - left side of interval
        b - right side of interval
        w - bspline workspace

Return: success or error

Notes: 1) w->knots is modified to contain the uniformly spaced
          knots

       2) The knots vector is set up as follows (using octave syntax):

          knots(1:k) = a
          knots(k+1:k+l-1) = a + i*delta, i = 1 .. l - 1
          knots(n+1:n+k) = b
*/

int
gsl_bspline_knots_uniform (const double a, const double b,
			   gsl_bspline_workspace * w)
{
  size_t i;			/* looping */
  double delta;			/* interval spacing */
  double x;

  delta = (b - a) / (double) w->l;

  for (i = 0; i < w->k; i++)
    gsl_vector_set (w->knots, i, a);

  x = a + delta;
  for (i = 0; i < w->l - 1; i++)
    {
      gsl_vector_set (w->knots, w->k + i, x);
      x += delta;
    }

  for (i = w->n; i < w->n + w->k; i++)
    gsl_vector_set (w->knots, i, b);

  return GSL_SUCCESS;
}				/* gsl_bspline_knots_uniform() */

/*
gsl_bspline_eval()
  Evaluate the basis functions B_i(x) for all i. This is
a wrapper function for gsl_bspline_eval_nonzero() which
formats the output in a nice way.

Inputs: x - point for evaluation
        B - (output) where to store B_i(x) values
            the length of this vector is
            n = nbreak + k - 2 = l + k - 1 = w->n
        w - bspline workspace

Return: success or error

Notes: The w->knots vector must be initialized prior to calling
       this function (see gsl_bspline_knots())
*/

int
gsl_bspline_eval (const double x, gsl_vector * B, gsl_bspline_workspace * w)
{
  if (B->size != w->n)
    {
      GSL_ERROR ("vector B not of length n", GSL_EBADLEN);
    }
  else
    {
      size_t i;			/* looping */
      size_t istart;		/* first non-zero spline for x */
      size_t iend;		/* last non-zero spline for x, knot for x */
      int error;		/* error handling */

      /* find all non-zero B_i(x) values */
      error = gsl_bspline_eval_nonzero (x, w->B, &istart, &iend, w);
      if (error)
	{
	  return error;
	}

      /* store values in appropriate part of given vector */
      for (i = 0; i < istart; i++)
	gsl_vector_set (B, i, 0.0);

      for (i = istart; i <= iend; i++)
	gsl_vector_set (B, i, gsl_vector_get (w->B, i - istart));

      for (i = iend + 1; i < w->n; i++)
	gsl_vector_set (B, i, 0.0);

      return GSL_SUCCESS;
    }
}				/* gsl_bspline_eval() */

/*
gsl_bspline_eval_nonzero()
  Evaluate all non-zero B-spline functions at point x.
These are the B_i(x) for i in [istart, iend].
Always B_i(x) = 0 for i < istart and for i > iend.

Inputs: x      - point at which to evaluate splines
        Bk     - (output) where to store B-spline values (length k)
        istart - (output) B-spline function index of
                 first non-zero basis for given x
        iend   - (output) B-spline function index of
                 last non-zero basis for given x.
                 This is also the knot index corresponding to x.
        w      - bspline workspace

Return: success or error

Notes: 1) the w->knots vector must be initialized before calling
          this function

       2) On output, B contains

             [B_{istart,k}, B_{istart+1,k},
             ..., B_{iend-1,k}, B_{iend,k}]

          evaluated at the given x.
*/

int
gsl_bspline_eval_nonzero (const double x, gsl_vector * Bk, size_t * istart,
			  size_t * iend, gsl_bspline_workspace * w)
{
  if (Bk->size != w->k)
    {
      GSL_ERROR ("Bk vector length does not match order k", GSL_EBADLEN);
    }
  else
    {
      size_t i;			/* spline index */
      size_t j;			/* looping */
      int flag = 0;		/* interval search flag */
      int error = 0;		/* error flag */

      i = bspline_find_interval (x, &flag, w);
      error = bspline_process_interval_for_eval (x, &i, flag, w);
      if (error)
	{
	  return error;
	}

      *istart = i - w->k + 1;
      *iend = i;

      bspline_pppack_bsplvb (w->knots, w->k, 1, x, *iend, &j, w->deltal,
			     w->deltar, Bk);

      return GSL_SUCCESS;
    }
}				/* gsl_bspline_eval_nonzero() */

/*
gsl_bspline_deriv_eval()
  Evaluate d^j/dx^j B_i(x) for all i, 0 <= j <= nderiv.
This is a wrapper function for gsl_bspline_deriv_eval_nonzero()
which formats the output in a nice way.

Inputs: x      - point for evaluation
        nderiv - number of derivatives to compute, inclusive.
        dB     - (output) where to store d^j/dx^j B_i(x)
                 values. the size of this matrix is
                 (n = nbreak + k - 2 = l + k - 1 = w->n)
                 by (nderiv + 1)
        w      - bspline derivative workspace

Return: success or error

Notes: 1) The w->knots vector must be initialized prior to calling
          this function (see gsl_bspline_knots())

       2) based on PPPACK's bsplvd
*/

int
gsl_bspline_deriv_eval (const double x, const size_t nderiv, gsl_matrix * dB,
			gsl_bspline_workspace * w,
			gsl_bspline_deriv_workspace * dw)
{
  if (dB->size1 != w->n)
    {
      GSL_ERROR ("dB matrix first dimension not of length n", GSL_EBADLEN);
    }
  else if (dB->size2 < nderiv + 1)
    {
      GSL_ERROR
	("dB matrix second dimension must be at least length nderiv+1",
	 GSL_EBADLEN);
    }
  else if (dw->k < w->k) 
    {
      GSL_ERROR ("derivative workspace is too small", GSL_EBADLEN);
    }
  else
    {
      size_t i;			/* looping */
      size_t j;			/* looping */
      size_t istart;		/* first non-zero spline for x */
      size_t iend;		/* last non-zero spline for x, knot for x */
      int error;		/* error handling */

      /* find all non-zero d^j/dx^j B_i(x) values */
      error =
	gsl_bspline_deriv_eval_nonzero (x, nderiv, dw->dB, &istart, &iend, w,
					dw);
      if (error)
	{
	  return error;
	}

      /* store values in appropriate part of given matrix */
      for (j = 0; j <= nderiv; j++)
	{
	  for (i = 0; i < istart; i++)
	    gsl_matrix_set (dB, i, j, 0.0);

	  for (i = istart; i <= iend; i++)
	    gsl_matrix_set (dB, i, j, gsl_matrix_get (dw->dB, i - istart, j));

	  for (i = iend + 1; i < w->n; i++)
	    gsl_matrix_set (dB, i, j, 0.0);
	}

      return GSL_SUCCESS;
    }
}				/* gsl_bspline_deriv_eval() */

/*
gsl_bspline_deriv_eval_nonzero()
  At point x evaluate all requested, non-zero B-spline function
derivatives and store them in dB.  These are the
d^j/dx^j B_i(x) with i in [istart, iend] and j in [0, nderiv].
Always d^j/dx^j B_i(x) = 0 for i < istart and for i > iend.

Inputs: x      - point at which to evaluate splines
        nderiv - number of derivatives to request, inclusive
        dB     - (output) where to store dB-spline derivatives
                 (size k by nderiv + 1)
        istart - (output) B-spline function index of
                 first non-zero basis for given x
        iend   - (output) B-spline function index of
                 last non-zero basis for given x.
                 This is also the knot index corresponding to x.
        w      - bspline derivative workspace

Return: success or error

Notes: 1) the w->knots vector must be initialized before calling
          this function

       2) On output, dB contains

            [[B_{istart,  k}, ..., d^nderiv/dx^nderiv B_{istart  ,k}],
             [B_{istart+1,k}, ..., d^nderiv/dx^nderiv B_{istart+1,k}],
             ...
             [B_{iend-1,  k}, ..., d^nderiv/dx^nderiv B_{iend-1,  k}],
             [B_{iend,    k}, ..., d^nderiv/dx^nderiv B_{iend,    k}]]

          evaluated at x.  B_{istart, k} is stored in dB(0,0).
          Each additional column contains an additional derivative.

       3) Note that the zero-th column of the result contains the
          0th derivative, which is simply a function evaluation.

       4) based on PPPACK's bsplvd
*/

int
gsl_bspline_deriv_eval_nonzero (const double x, const size_t nderiv,
				gsl_matrix * dB, size_t * istart,
				size_t * iend, gsl_bspline_workspace * w,
				gsl_bspline_deriv_workspace * dw)
{
  if (dB->size1 != w->k)
    {
      GSL_ERROR ("dB matrix first dimension not of length k", GSL_EBADLEN);
    }
  else if (dB->size2 < nderiv + 1)
    {
      GSL_ERROR
	("dB matrix second dimension must be at least length nderiv+1",
	 GSL_EBADLEN);
    }
  else if (dw->k < w->k) 
    {
      GSL_ERROR ("derivative workspace is too small", GSL_EBADLEN);
    }
  else
    {
      size_t i;			/* spline index */
      size_t j;			/* looping */
      int flag = 0;		/* interval search flag */
      int error = 0;		/* error flag */
      size_t min_nderivk;

      i = bspline_find_interval (x, &flag, w);
      error = bspline_process_interval_for_eval (x, &i, flag, w);
      if (error)
	{
	  return error;
	}

      *istart = i - w->k + 1;
      *iend = i;

      bspline_pppack_bsplvd (w->knots, w->k, x, *iend,
			     w->deltal, w->deltar, dw->A, dB, nderiv);

      /* An order k b-spline has at most k-1 nonzero derivatives
         so we need to zero all requested higher order derivatives */
      min_nderivk = GSL_MIN_INT (nderiv, w->k - 1);
      for (j = min_nderivk + 1; j <= nderiv; j++)
	{
	  for (i = 0; i < w->k; i++)
	    {
	      gsl_matrix_set (dB, i, j, 0.0);
	    }
	}

      return GSL_SUCCESS;
    }
}				/* gsl_bspline_deriv_eval_nonzero() */

/****************************************
 *          INTERNAL ROUTINES           *
 ****************************************/

/*
bspline_find_interval()
  Find knot interval such that t_i <= x < t_{i + 1}
where the t_i are knot values.

Inputs: x    - x value
        flag - (output) error flag
        w    - bspline workspace

Return: i (index in w->knots corresponding to left limit of interval)

Notes: The error conditions are reported as follows:

       Condition                        Return value        Flag
       ---------                        ------------        ----
       x < t_0                               0               -1
       t_i <= x < t_{i+1}                    i                0
       t_i < x = t_{i+1} = t_{n+k-1}         i                0
       t_{n+k-1} < x                       l+k-1             +1
*/

static inline size_t
bspline_find_interval (const double x, int *flag, gsl_bspline_workspace * w)
{
  size_t i;

  if (x < gsl_vector_get (w->knots, 0))
    {
      *flag = -1;
      return 0;
    }

  /* find i such that t_i <= x < t_{i+1} */
  for (i = w->k - 1; i < w->k + w->l - 1; i++)
    {
      const double ti = gsl_vector_get (w->knots, i);
      const double tip1 = gsl_vector_get (w->knots, i + 1);

      if (tip1 < ti)
	{
	  GSL_ERROR ("knots vector is not increasing", GSL_EINVAL);
	}

      if (ti <= x && x < tip1)
	break;

      if (ti < x && x == tip1 && tip1 == gsl_vector_get (w->knots, w->k + w->l
							 - 1))
	break;
    }

  if (i == w->k + w->l - 1)
    *flag = 1;
  else
    *flag = 0;

  return i;
}				/* bspline_find_interval() */

/*
bspline_process_interval_for_eval()
  Consumes an x location, left knot from bspline_find_interval, flag
from bspline_find_interval, and a workspace.  Checks that x lies within
the splines' knots, enforces some endpoint continuity requirements, and
avoids divide by zero errors in the underlying bspline_pppack_* functions.
*/
static inline int
bspline_process_interval_for_eval (const double x, size_t * i, const int flag,
				   gsl_bspline_workspace * w)
{
  if (flag == -1)
    {
      GSL_ERROR ("x outside of knot interval", GSL_EINVAL);
    }
  else if (flag == 1)
    {
      if (x <= gsl_vector_get (w->knots, *i) + GSL_DBL_EPSILON)
	{
	  *i -= 1;
	}
      else
	{
	  GSL_ERROR ("x outside of knot interval", GSL_EINVAL);
	}
    }

  if (gsl_vector_get (w->knots, *i) == gsl_vector_get (w->knots, *i + 1))
    {
      GSL_ERROR ("knot(i) = knot(i+1) will result in division by zero",
		 GSL_EINVAL);
    }

  return GSL_SUCCESS;
}				/* bspline_process_interval_for_eval */

/********************************************************************
 * PPPACK ROUTINES
 *
 * The routines herein deliberately avoid using the bspline workspace,
 * choosing instead to pass all work areas explicitly.  This allows
 * others to more easily adapt these routines to low memory or
 * parallel scenarios.
 ********************************************************************/

/*
bspline_pppack_bsplvb()
  calculates the value of all possibly nonzero b-splines at x of order
jout = max( jhigh , (j+1)*(index-1) ) with knot sequence t.

Parameters:
   t      - knot sequence, of length left + jout , assumed to be
            nondecreasing.  assumption t(left).lt.t(left + 1).
            division by zero  will result if t(left) = t(left+1)
   jhigh  -
   index  - integers which determine the order jout = max(jhigh,
            (j+1)*(index-1))  of the b-splines whose values at x
            are to be returned.  index  is used to avoid
            recalculations when several columns of the triangular
            array of b-spline values are needed (e.g., in  bsplpp
            or in  bsplvd ).  precisely,

            if  index = 1 ,
               the calculation starts from scratch and the entire
               triangular array of b-spline values of orders
               1,2,...,jhigh  is generated order by order , i.e.,
               column by column .

            if  index = 2 ,
               only the b-spline values of order j+1, j+2, ..., jout
               are generated, the assumption being that biatx, j,
               deltal, deltar are, on entry, as they were on exit
               at the previous call.

            in particular, if jhigh = 0, then jout = j+1, i.e.,
            just the next column of b-spline values is generated.
   x      - the point at which the b-splines are to be evaluated.
   left   - an integer chosen (usually) so that
            t(left) .le. x .le. t(left+1).
   j      - (output) a working scalar for indexing
   deltal - (output) a working area which must be of length at least jout
   deltar - (output) a working area which must be of length at least jout
   biatx  - (output) array of length jout, with  biatx(i)
            containing the value at  x  of the polynomial of order
            jout which agrees with the b-spline b(left-jout+i,jout,t)
            on the interval (t(left), t(left+1)) .

Method:
   the recurrence relation

                      x - t(i)              t(i+j+1) - x
      b(i,j+1)(x) = -----------b(i,j)(x) + ---------------b(i+1,j)(x)
                    t(i+j)-t(i)            t(i+j+1)-t(i+1)

   is used (repeatedly) to generate the (j+1)-vector  b(left-j,j+1)(x),
   ...,b(left,j+1)(x)  from the j-vector  b(left-j+1,j)(x),...,
   b(left,j)(x), storing the new values in  biatx  over the old. the
   facts that

      b(i,1) = 1  if  t(i) .le. x .lt. t(i+1)

   and that

      b(i,j)(x) = 0  unless  t(i) .le. x .lt. t(i+j)

   are used. the particular organization of the calculations follows
   algorithm (8) in chapter x of [1].

Notes:

   (1) This is a direct translation of PPPACK's bsplvb routine with
       j, deltal, deltar rewritten as input parameters and
       utilizing zero-based indexing.

   (2) This routine contains no error checking.  Please use routines
       like gsl_bspline_eval().
*/

static void
bspline_pppack_bsplvb (const gsl_vector * t,
		       const size_t jhigh,
		       const size_t index,
		       const double x,
		       const size_t left,
		       size_t * j,
		       gsl_vector * deltal,
		       gsl_vector * deltar, gsl_vector * biatx)
{
  size_t i;			/* looping */
  double saved;
  double term;

  if (index == 1)
    {
      *j = 0;
      gsl_vector_set (biatx, 0, 1.0);
    }

  for ( /* NOP */ ; *j < jhigh - 1; *j += 1)
    {
      gsl_vector_set (deltar, *j, gsl_vector_get (t, left + *j + 1) - x);
      gsl_vector_set (deltal, *j, x - gsl_vector_get (t, left - *j));

      saved = 0.0;

      for (i = 0; i <= *j; i++)
	{
	  term = gsl_vector_get (biatx, i) / (gsl_vector_get (deltar, i)
					      + gsl_vector_get (deltal,
								*j - i));

	  gsl_vector_set (biatx, i,
			  saved + gsl_vector_get (deltar, i) * term);

	  saved = gsl_vector_get (deltal, *j - i) * term;
	}

      gsl_vector_set (biatx, *j + 1, saved);
    }

  return;
}				/* gsl_bspline_pppack_bsplvb */

/*
bspline_pppack_bsplvd()
  calculates value and derivs of all b-splines which do not vanish at x

Parameters:
   t      - the knot array, of length left+k (at least)
   k      - the order of the b-splines to be evaluated
   x      - the point at which these values are sought
   left   - an integer indicating the left endpoint of the interval
            of interest. the k b-splines whose support contains the
            interval (t(left), t(left+1)) are to be considered.
            it is assumed that t(left) .lt. t(left+1)
            division by zero will result otherwise (in  bsplvb).
            also, the output is as advertised only if
            t(left) .le. x .le. t(left+1) .
   deltal - a working area which must be of length at least k
   deltar - a working area which must be of length at least k
   a      - an array of order (k,k), to contain b-coeffs of the
            derivatives of a certain order of the k b-splines
            of interest.
   dbiatx - an array of order (k,nderiv). its entry (i,m) contains
            value of (m)th derivative of (left-k+i)-th b-spline
            of order k for knot sequence  t, i=1,...,k, m=0,...,nderiv.
   nderiv - an integer indicating that values of b-splines and
            their derivatives up to AND INCLUDING the nderiv-th
            are asked for. (nderiv is replaced internally by the
            integer mhigh in (1,k) closest to it.)

Method:
   values at x of all the relevant b-splines of order k,k-1,..., k+1-nderiv
   are generated via bsplvb and stored temporarily in dbiatx.  then, the
   b-coeffs of the required derivatives of the b-splines of interest are
   generated by differencing, each from the preceeding one of lower order,
   and combined with the values of b-splines of corresponding order in
   dbiatx  to produce the desired values .

Notes:

   (1) This is a direct translation of PPPACK's bsplvd routine with
       deltal, deltar rewritten as input parameters (to later feed them
       to bspline_pppack_bsplvb) and utilizing zero-based indexing.

   (2) This routine contains no error checking.
*/

static void
bspline_pppack_bsplvd (const gsl_vector * t,
		       const size_t k,
		       const double x,
		       const size_t left,
		       gsl_vector * deltal,
		       gsl_vector * deltar,
		       gsl_matrix * a,
		       gsl_matrix * dbiatx, const size_t nderiv)
{
  int i, ideriv, il, j, jlow, jp1mid, kmm, ldummy, m, mhigh;
  double factor, fkmm, sum;

  size_t bsplvb_j;
  gsl_vector_view dbcol = gsl_matrix_column (dbiatx, 0);

  mhigh = GSL_MIN_INT (nderiv, k - 1);
  bspline_pppack_bsplvb (t, k - mhigh, 1, x, left, &bsplvb_j, deltal, deltar,
			 &dbcol.vector);
  if (mhigh > 0)
    {
      /* the first column of dbiatx always contains the b-spline
         values for the current order. these are stored in column
         k-current order before bsplvb is called to put values
         for the next higher order on top of it.  */
      ideriv = mhigh;
      for (m = 1; m <= mhigh; m++)
	{
	  for (j = ideriv, jp1mid = 0; j < (int) k; j++, jp1mid++)
	    {
	      gsl_matrix_set (dbiatx, j, ideriv,
			      gsl_matrix_get (dbiatx, jp1mid, 0));
	    }
	  ideriv--;
	  bspline_pppack_bsplvb (t, k - ideriv, 2, x, left, &bsplvb_j, deltal,
				 deltar, &dbcol.vector);
	}

      /* at this point,  b(left-k+i, k+1-j)(x) is in  dbiatx(i,j)
         for i=j,...,k-1 and j=0,...,mhigh. in particular, the
         first column of dbiatx is already in final form. to obtain
         corresponding derivatives of b-splines in subsequent columns,
         generate their b-repr. by differencing, then evaluate at x. */
      jlow = 0;
      for (i = 0; i < (int) k; i++)
	{
	  for (j = jlow; j < (int) k; j++)
	    {
	      gsl_matrix_set (a, j, i, 0.0);
	    }
	  jlow = i;
	  gsl_matrix_set (a, i, i, 1.0);
	}

      /* at this point, a(.,j) contains the b-coeffs for the j-th of the
         k b-splines of interest here. */
      for (m = 1; m <= mhigh; m++)
	{
	  kmm = k - m;
	  fkmm = (float) kmm;
	  il = left;
	  i = k - 1;

	  /* for j=1,...,k, construct b-coeffs of (m)th  derivative
	     of b-splines from those for preceding derivative by
	     differencing and store again in  a(.,j) . the fact that
	     a(i,j) = 0  for i .lt. j  is used. */
	  for (ldummy = 0; ldummy < kmm; ldummy++)
	    {
	      factor =
		fkmm / (gsl_vector_get (t, il + kmm) -
			gsl_vector_get (t, il));
	      /* the assumption that t(left).lt.t(left+1) makes
	         denominator in factor nonzero. */
	      for (j = 0; j <= i; j++)
		{
		  gsl_matrix_set (a, i, j, factor * (gsl_matrix_get (a, i, j)
						     - gsl_matrix_get (a,
								       i - 1,
								       j)));
		}
	      il--;
	      i--;
	    }

	  /* for i=1,...,k, combine b-coeffs a(.,i) with b-spline values
	     stored in dbiatx(.,m) to get value of (m)th  derivative
	     of i-th b-spline (of interest here) at x, and store in
	     dbiatx(i,m). storage of this value over the value of a
	     b-spline of order m there is safe since the remaining
	     b-spline derivatives of the same order do not use this
	     value due to the fact that a(j,i) = 0 for j .lt. i . */
	  for (i = 0; i < (int) k; i++)
	    {
	      sum = 0;
	      jlow = GSL_MAX_INT (i, m);
	      for (j = jlow; j < (int) k; j++)
		{
		  sum +=
		    gsl_matrix_get (a, j, i) * gsl_matrix_get (dbiatx, j, m);
		}
	      gsl_matrix_set (dbiatx, i, m, sum);
	    }
	}
    }

  return;

}				/* bspline_pppack_bsplvd */
