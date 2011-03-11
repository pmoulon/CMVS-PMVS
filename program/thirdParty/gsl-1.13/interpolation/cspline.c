/* interpolation/cspline.c
 * 
 * Copyright (C) 1996, 1997, 1998, 1999, 2000, 2001, 2002, 2004 Gerard Jungman
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

/* Author:  G. Jungman
 */
#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_vector.h>
#include "integ_eval.h"
#include <gsl/gsl_interp.h>

typedef struct
{
  double * c;
  double * g;
  double * diag;
  double * offdiag;
} cspline_state_t;


/* common initialization */
static void *
cspline_alloc (size_t size)
{
  cspline_state_t * state = (cspline_state_t *) malloc (sizeof (cspline_state_t));

  if (state == NULL)
    {
      GSL_ERROR_NULL("failed to allocate space for state", GSL_ENOMEM);
    }
  
  state->c = (double *) malloc (size * sizeof (double));
  
  if (state->c == NULL)
    {
      free (state);
      GSL_ERROR_NULL("failed to allocate space for c", GSL_ENOMEM);
    }

  state->g = (double *) malloc (size * sizeof (double));
  
  if (state->g == NULL)
    {
      free (state->c);
      free (state);
      GSL_ERROR_NULL("failed to allocate space for g", GSL_ENOMEM);
    }

  state->diag = (double *) malloc (size * sizeof (double));
  
  if (state->diag == NULL)
    {
      free (state->g);
      free (state->c);
      free (state);
      GSL_ERROR_NULL("failed to allocate space for diag", GSL_ENOMEM);
    }

  state->offdiag = (double *) malloc (size * sizeof (double));
  
  if (state->offdiag == NULL)
    {
      free (state->diag);
      free (state->g);
      free (state->c);
      free (state);
      GSL_ERROR_NULL("failed to allocate space for offdiag", GSL_ENOMEM);
    }

  return state;
}


/* natural spline calculation
 * see [Engeln-Mullges + Uhlig, p. 254]
 */
static int
cspline_init (void * vstate, const double xa[], const double ya[],
              size_t size)
{
  cspline_state_t *state = (cspline_state_t *) vstate;

  size_t i;
  size_t num_points = size;
  size_t max_index = num_points - 1;  /* Engeln-Mullges + Uhlig "n" */
  size_t sys_size = max_index - 1;    /* linear system is sys_size x sys_size */

  state->c[0] = 0.0;
  state->c[max_index] = 0.0;

  for (i = 0; i < sys_size; i++)
    {
      const double h_i   = xa[i + 1] - xa[i];
      const double h_ip1 = xa[i + 2] - xa[i + 1];
      const double ydiff_i   = ya[i + 1] - ya[i];
      const double ydiff_ip1 = ya[i + 2] - ya[i + 1];
      const double g_i = (h_i != 0.0) ? 1.0 / h_i : 0.0;
      const double g_ip1 = (h_ip1 != 0.0) ? 1.0 / h_ip1 : 0.0;
      state->offdiag[i] = h_ip1;
      state->diag[i] = 2.0 * (h_ip1 + h_i);
      state->g[i] = 3.0 * (ydiff_ip1 * g_ip1 -  ydiff_i * g_i);
    }

  if (sys_size == 1)
    {
      state->c[1] = state->g[0] / state->diag[0];
      return GSL_SUCCESS;
    }
  else
    {
      gsl_vector_view g_vec = gsl_vector_view_array(state->g, sys_size);
      gsl_vector_view diag_vec = gsl_vector_view_array(state->diag, sys_size);
      gsl_vector_view offdiag_vec = gsl_vector_view_array(state->offdiag, sys_size - 1);
      gsl_vector_view solution_vec = gsl_vector_view_array ((state->c) + 1, sys_size);
      
      int status = gsl_linalg_solve_symm_tridiag(&diag_vec.vector, 
                                                 &offdiag_vec.vector, 
                                                 &g_vec.vector, 
                                                 &solution_vec.vector);
      return status;
    }
}


/* periodic spline calculation
 * see [Engeln-Mullges + Uhlig, p. 256]
 */
static int
cspline_init_periodic (void * vstate, const double xa[], const double ya[],
                       size_t size)
{
  cspline_state_t *state = (cspline_state_t *) vstate;

  size_t i;
  size_t num_points = size;
  size_t max_index = num_points - 1;  /* Engeln-Mullges + Uhlig "n" */
  size_t sys_size = max_index;    /* linear system is sys_size x sys_size */

  if (sys_size == 2) {
    /* solve 2x2 system */
    
    const double h0 = xa[1] - xa[0];
    const double h1 = xa[2] - xa[1];

    const double A = 2.0*(h0 + h1);
    const double B = h0 + h1;
    double g[2];
    double det;
    
    g[0] = 3.0 * ((ya[2] - ya[1]) / h1 - (ya[1] - ya[0]) / h0);
    g[1] = 3.0 * ((ya[1] - ya[2]) / h0 - (ya[2] - ya[1]) / h1);
    
    det = 3.0 * (h0 + h1) * (h0 + h1);
    state->c[1] = ( A * g[0] - B * g[1])/det;
    state->c[2] = (-B * g[0] + A * g[1])/det;
    state->c[0] = state->c[2];
    
    return GSL_SUCCESS;
  } else {
    
    for (i = 0; i < sys_size-1; i++) {
      const double h_i       = xa[i + 1] - xa[i];
      const double h_ip1     = xa[i + 2] - xa[i + 1];
      const double ydiff_i   = ya[i + 1] - ya[i];
      const double ydiff_ip1 = ya[i + 2] - ya[i + 1];
      const double g_i = (h_i != 0.0) ? 1.0 / h_i : 0.0;
      const double g_ip1 = (h_ip1 != 0.0) ? 1.0 / h_ip1 : 0.0;
      state->offdiag[i] = h_ip1;
      state->diag[i] = 2.0 * (h_ip1 + h_i);
      state->g[i] = 3.0 * (ydiff_ip1 * g_ip1 - ydiff_i * g_i);
    }

    i = sys_size - 1;

    {
      const double h_i       = xa[i + 1] - xa[i];
      const double h_ip1     = xa[1] - xa[0];
      const double ydiff_i   = ya[i + 1] - ya[i];
      const double ydiff_ip1 = ya[1] - ya[0];
      const double g_i = (h_i != 0.0) ? 1.0 / h_i : 0.0;
      const double g_ip1 = (h_ip1 != 0.0) ? 1.0 / h_ip1 : 0.0;
      state->offdiag[i] = h_ip1;
      state->diag[i] = 2.0 * (h_ip1 + h_i);
      state->g[i] = 3.0 * (ydiff_ip1 * g_ip1 - ydiff_i * g_i);
    }
    
    {
      gsl_vector_view g_vec = gsl_vector_view_array(state->g, sys_size);
      gsl_vector_view diag_vec = gsl_vector_view_array(state->diag, sys_size);
      gsl_vector_view offdiag_vec = gsl_vector_view_array(state->offdiag, sys_size);
      gsl_vector_view solution_vec = gsl_vector_view_array ((state->c) + 1, sys_size);
      
      int status = gsl_linalg_solve_symm_cyc_tridiag(&diag_vec.vector, 
                                                     &offdiag_vec.vector, 
                                                     &g_vec.vector, 
                                                     &solution_vec.vector);
      state->c[0] = state->c[max_index];
      
      return status;
    }
  }
}


static
void
cspline_free (void * vstate)
{
  cspline_state_t *state = (cspline_state_t *) vstate;
  
  free (state->c);
  free (state->g);
  free (state->diag);
  free (state->offdiag);
  free (state);
}

/* function for common coefficient determination
 */
static inline void
coeff_calc (const double c_array[], double dy, double dx, size_t index,  
            double * b, double * c, double * d)
{
  const double c_i = c_array[index];
  const double c_ip1 = c_array[index + 1];
  *b = (dy / dx) - dx * (c_ip1 + 2.0 * c_i) / 3.0;
  *c = c_i;
  *d = (c_ip1 - c_i) / (3.0 * dx);
}


static
int
cspline_eval (const void * vstate,
              const double x_array[], const double y_array[], size_t size,
              double x,
              gsl_interp_accel * a,
              double *y)
{
  const cspline_state_t *state = (const cspline_state_t *) vstate;

  double x_lo, x_hi;
  double dx;
  size_t index;
  
  if (a != 0)
    {
      index = gsl_interp_accel_find (a, x_array, size, x);
    }
  else
    {
      index = gsl_interp_bsearch (x_array, x, 0, size - 1);
    }
  
  /* evaluate */
  x_hi = x_array[index + 1];
  x_lo = x_array[index];
  dx = x_hi - x_lo;
  if (dx > 0.0)
    {
      const double y_lo = y_array[index];
      const double y_hi = y_array[index + 1];
      const double dy = y_hi - y_lo;
      double delx = x - x_lo;
      double b_i, c_i, d_i; 
      coeff_calc(state->c, dy, dx, index,  &b_i, &c_i, &d_i);
      *y = y_lo + delx * (b_i + delx * (c_i + delx * d_i));
      return GSL_SUCCESS;
    }
  else
    {
      *y = 0.0;
      return GSL_EINVAL;
    }
}


static
int
cspline_eval_deriv (const void * vstate,
                    const double x_array[], const double y_array[], size_t size,
                    double x,
                    gsl_interp_accel * a,
                    double *dydx)
{
  const cspline_state_t *state = (const cspline_state_t *) vstate;

  double x_lo, x_hi;
  double dx;
  size_t index;
  
  if (a != 0)
    {
      index = gsl_interp_accel_find (a, x_array, size, x);
    }
  else
    {
      index = gsl_interp_bsearch (x_array, x, 0, size - 1);
    }
  
  /* evaluate */
  x_hi = x_array[index + 1];
  x_lo = x_array[index];
  dx = x_hi - x_lo;
  if (dx > 0.0)
    {
      const double y_lo = y_array[index];
      const double y_hi = y_array[index + 1];
      const double dy = y_hi - y_lo;
      double delx = x - x_lo;
      double b_i, c_i, d_i; 
      coeff_calc(state->c, dy, dx, index,  &b_i, &c_i, &d_i);
      *dydx = b_i + delx * (2.0 * c_i + 3.0 * d_i * delx);
      return GSL_SUCCESS;
    }
  else
    {
      *dydx = 0.0;
      return GSL_FAILURE;
    }
}


static
int
cspline_eval_deriv2 (const void * vstate,
                     const double x_array[], const double y_array[], size_t size,
                     double x,
                     gsl_interp_accel * a,
                     double * y_pp)
{
  const cspline_state_t *state = (const cspline_state_t *) vstate;

  double x_lo, x_hi;
  double dx;
  size_t index;
  
  if (a != 0)
    {
      index = gsl_interp_accel_find (a, x_array, size, x);
    }
  else
    {
      index = gsl_interp_bsearch (x_array, x, 0, size - 1);
    }
  
  /* evaluate */
  x_hi = x_array[index + 1];
  x_lo = x_array[index];
  dx = x_hi - x_lo;
  if (dx > 0.0)
    {
      const double y_lo = y_array[index];
      const double y_hi = y_array[index + 1];
      const double dy = y_hi - y_lo;
      double delx = x - x_lo;
      double b_i, c_i, d_i;
      coeff_calc(state->c, dy, dx, index,  &b_i, &c_i, &d_i);
      *y_pp = 2.0 * c_i + 6.0 * d_i * delx;
      return GSL_SUCCESS;
    }
  else
    {
      *y_pp = 0.0;
      return GSL_FAILURE;
    }
}


static
int
cspline_eval_integ (const void * vstate,
                    const double x_array[], const double y_array[], size_t size,
                    gsl_interp_accel * acc,
                    double a, double b,
                    double * result)
{
  const cspline_state_t *state = (const cspline_state_t *) vstate;

  size_t i, index_a, index_b;
  
  if (acc != 0)
    {
      index_a = gsl_interp_accel_find (acc, x_array, size, a);
      index_b = gsl_interp_accel_find (acc, x_array, size, b);
    }
  else
    {
      index_a = gsl_interp_bsearch (x_array, a, 0, size - 1);
      index_b = gsl_interp_bsearch (x_array, b, 0, size - 1);
    }

  *result = 0.0;
  
  /* interior intervals */
  for(i=index_a; i<=index_b; i++) {
    const double x_hi = x_array[i + 1];
    const double x_lo = x_array[i];
    const double y_lo = y_array[i];
    const double y_hi = y_array[i + 1];
    const double dx = x_hi - x_lo;
    const double dy = y_hi - y_lo;
    if(dx != 0.0) {
      double b_i, c_i, d_i; 
      coeff_calc(state->c, dy, dx, i,  &b_i, &c_i, &d_i);
      
      if (i == index_a || i == index_b)
        {
          double x1 = (i == index_a) ? a : x_lo;
          double x2 = (i == index_b) ? b : x_hi;
          *result += integ_eval(y_lo, b_i, c_i, d_i, x_lo, x1, x2);
        }
      else
        {
          *result += dx * (y_lo + dx*(0.5*b_i + dx*(c_i/3.0 + 0.25*d_i*dx)));
        }
    }
    else {
      *result = 0.0;
      return GSL_FAILURE;
    }
  }
  
  return GSL_SUCCESS;
}

static const gsl_interp_type cspline_type = 
{
  "cspline", 
  3,
  &cspline_alloc,
  &cspline_init,
  &cspline_eval,
  &cspline_eval_deriv,
  &cspline_eval_deriv2,
  &cspline_eval_integ,
  &cspline_free
};

const gsl_interp_type * gsl_interp_cspline = &cspline_type;

static const gsl_interp_type cspline_periodic_type = 
{
  "cspline-periodic", 
  2,
  &cspline_alloc,
  &cspline_init_periodic,
  &cspline_eval,
  &cspline_eval_deriv,
  &cspline_eval_deriv2,
  &cspline_eval_integ,
  &cspline_free
};

const gsl_interp_type * gsl_interp_cspline_periodic = &cspline_periodic_type;


