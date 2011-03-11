/* vector/oper_source.c
 * 
 * Copyright (C) 2008 Brian Gough
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

int 
FUNCTION(gsl_vector, add) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
  const size_t N = a->size;

  if (b->size != N)
    {
      GSL_ERROR ("vectors must have same length", GSL_EBADLEN);
    }
  else 
    {
      const size_t stride_a = a->stride;
      const size_t stride_b = b->stride;

      size_t i;

      for (i = 0; i < N; i++)
        {
          a->data[2 * i * stride_a] += b->data[2 * i * stride_b];
          a->data[2 * i * stride_a + 1] += b->data[2 * i * stride_b + 1];
        }
      
      return GSL_SUCCESS;
    }
}

int 
FUNCTION(gsl_vector, sub) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
  const size_t N = a->size;

  if (b->size != N)
    {
      GSL_ERROR ("vectors must have same length", GSL_EBADLEN);
    }
  else 
    {
      const size_t stride_a = a->stride;
      const size_t stride_b = b->stride;

      size_t i;

      for (i = 0; i < N; i++)
        {
          a->data[2 * i * stride_a] -= b->data[2 * i * stride_b];
          a->data[2 * i * stride_a + 1] -= b->data[2 * i * stride_b + 1];
        }
      
      return GSL_SUCCESS;
    }
}

int 
FUNCTION(gsl_vector, mul) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
  const size_t N = a->size;

  if (b->size != N)
    {
      GSL_ERROR ("vectors must have same length", GSL_EBADLEN);
    }
  else 
    {
      const size_t stride_a = a->stride;
      const size_t stride_b = b->stride;

      size_t i;

      for (i = 0; i < N; i++)
        {
          ATOMIC ar = a->data[2 * i * stride_a];
          ATOMIC ai = a->data[2 * i * stride_a + 1];
          
          ATOMIC br = b->data[2 * i * stride_b];
          ATOMIC bi = b->data[2 * i * stride_b + 1];

          a->data[2 * i * stride_a] = ar * br - ai * bi;
          a->data[2 * i * stride_a + 1] = ar * bi + ai * br;
        }
      
      return GSL_SUCCESS;
    }
}

int 
FUNCTION(gsl_vector, div) (TYPE(gsl_vector) * a, const TYPE(gsl_vector) * b)
{
  const size_t N = a->size;

  if (b->size != N)
    {
      GSL_ERROR ("vectors must have same length", GSL_EBADLEN);
    }
  else 
    {
      const size_t stride_a = a->stride;
      const size_t stride_b = b->stride;

      size_t i;

      for (i = 0; i < N; i++)
        {
          ATOMIC ar = a->data[2 * i * stride_a];
          ATOMIC ai = a->data[2 * i * stride_a + 1];
          
          ATOMIC br = b->data[2 * i * stride_b];
          ATOMIC bi = b->data[2 * i * stride_b + 1];

          ATOMIC s = 1.0 / hypot(br, bi);
          
          ATOMIC sbr = s * br;
          ATOMIC sbi = s * bi;
          
          a->data[2 * i * stride_a] = (ar * sbr + ai * sbi) * s;
          a->data[2 * i * stride_a + 1] = (ai * sbr - ar * sbi) * s;
        }
      
      return GSL_SUCCESS;
    }
}

int 
FUNCTION(gsl_vector, scale) (TYPE(gsl_vector) * a, const BASE x)
{
  const size_t N = a->size;
  const size_t stride = a->stride;
  
  size_t i;
  
  ATOMIC xr = GSL_REAL(x);
  ATOMIC xi = GSL_IMAG(x);

  for (i = 0; i < N; i++)
    {
      ATOMIC ar = a->data[2 * i * stride];
      ATOMIC ai = a->data[2 * i * stride + 1];
          
      a->data[2 * i * stride] = ar * xr - ai * xi;
      a->data[2 * i * stride + 1] = ar * xi + ai * xr;
    }
  
  return GSL_SUCCESS;
}

int 
FUNCTION(gsl_vector, add_constant) (TYPE(gsl_vector) * a, const BASE x)
{
  const size_t N = a->size;
  const size_t stride = a->stride;
  
  size_t i;

  ATOMIC xr = GSL_REAL(x);
  ATOMIC xi = GSL_IMAG(x);
  
  for (i = 0; i < N; i++)
    {
      a->data[2 * i * stride] += xr;
      a->data[2 * i * stride + 1] += xi;
    }
  
  return GSL_SUCCESS;
}
