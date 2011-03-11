/* sys/ldfrexp.c
 * 
 * Copyright (C) 2002, Gert Van den Eynde
 * Copyright (C) 2007, Brian Gough
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
#include <math.h>
#include <gsl/gsl_math.h>

double
gsl_ldexp (const double x, const int e)
{
  int ex;
  
  if (x == 0.0)
    {
      return x;
    }

  {
    double y = gsl_frexp (x, &ex);
    double e2 = e + ex, p2;
    
    if (e2 >= DBL_MAX_EXP)
      {
	y *= pow (2.0, e2 - DBL_MAX_EXP + 1);
	e2 = DBL_MAX_EXP - 1;
      }
    else if (e2 <= DBL_MIN_EXP)
      {
	y *= pow (2.0, e2 - DBL_MIN_EXP - 1);
	e2 = DBL_MIN_EXP + 1;
      }
    
    p2 = pow (2.0, e2);
    return y * p2;
  }
}

double
gsl_frexp (const double x, int *e)
{
  if (x == 0.0)
    {
      *e = 0;
      return 0.0;
    }
  else if (!finite (x))
    {
      *e = 0;
      return x;
    }
  else if (fabs (x) >= 0.5 && fabs (x) < 1)     /* Handle the common case */
    {
      *e = 0;
      return x;
    }
  else
    {
      double ex = ceil (log (fabs (x)) / M_LN2);
      int ei = (int) ex;
      double f;

      /* Prevent underflow and overflow of 2**(-ei) */
      if (ei < DBL_MIN_EXP)
        ei = DBL_MIN_EXP;

      if (ei > -DBL_MIN_EXP)
        ei = -DBL_MIN_EXP;

      f = x * pow (2.0, -ei);

      if (!finite (f))
        {
          /* This should not happen */
          *e = 0;
          return f;
        }

      while (fabs (f) >= 1.0)
        {
          ei++;
          f /= 2.0;
        }

      while (fabs (f) > 0 && fabs (f) < 0.5)
        {
          ei--;
          f *= 2.0;
        }

      *e = ei;
      return f;
    }
}
