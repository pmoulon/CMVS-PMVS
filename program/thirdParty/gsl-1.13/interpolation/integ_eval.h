/* interpolation/integ_eval_macro.h
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

/* function for doing the spline integral evaluation
   which is common to both the cspline and akima methods
 */

static inline double
integ_eval (double ai, double bi, double ci, double di, double xi, double a,
            double b)
{
  const double r1 = a - xi;
  const double r2 = b - xi;
  const double r12 = r1 + r2;
  const double bterm = 0.5 * bi * r12;
  const double cterm = (1.0 / 3.0) * ci * (r1 * r1 + r2 * r2 + r1 * r2);
  const double dterm = 0.25 * di * r12 * (r1 * r1 + r2 * r2);

  return (b - a) * (ai + bterm + cterm + dterm);
}
