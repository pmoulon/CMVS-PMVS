/* interpolation/accel.c
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

/* Author:  G. Jungman
 */
#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp.h>

gsl_interp_accel *
gsl_interp_accel_alloc (void)
{
  gsl_interp_accel *a = (gsl_interp_accel *) malloc (sizeof (gsl_interp_accel));
  if (a == 0)
    {
      GSL_ERROR_NULL("could not allocate space for gsl_interp_accel", GSL_ENOMEM);
    }

  a->cache = 0;
  a->hit_count = 0;
  a->miss_count = 0;

  return a;
}

int
gsl_interp_accel_reset (gsl_interp_accel * a)
{
  a->cache = 0;
  a->hit_count = 0;
  a->miss_count = 0;

  return GSL_SUCCESS;
}

void
gsl_interp_accel_free (gsl_interp_accel * a)
{
  RETURN_IF_NULL (a);
  free (a);
}
