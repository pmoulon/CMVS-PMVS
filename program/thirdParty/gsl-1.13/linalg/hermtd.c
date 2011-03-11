/* linalg/hermtd.c
 * 
 * Copyright (C) 2001, 2007, 2009 Brian Gough
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

/* Factorise a hermitian matrix A into
 *
 * A = U T U'
 *
 * where U is unitary and T is real symmetric tridiagonal.  Only the
 * diagonal and lower triangular part of A is referenced and modified.
 *
 * On exit, T is stored in the diagonal and first subdiagonal of
 * A. Since T is symmetric the upper diagonal is not stored.
 *
 * U is stored as a packed set of Householder transformations in the
 * lower triangular part of the input matrix below the first subdiagonal.
 *
 * The full matrix for U can be obtained as the product
 *
 *       U = U_N ... U_2 U_1
 *
 * where 
 *
 *       U_i = (I - tau_i * v_i * v_i')
 *
 * and where v_i is a Householder vector
 *
 *       v_i = [0, ..., 0, 1, A(i+2,i), A(i+3,i), ... , A(N,i)]
 *
 * This storage scheme is the same as in LAPACK.  See LAPACK's
 * chetd2.f for details.
 *
 * See Golub & Van Loan, "Matrix Computations" (3rd ed), Section 8.3 */

#include <config.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_complex_math.h>

#include <gsl/gsl_linalg.h>

int 
gsl_linalg_hermtd_decomp (gsl_matrix_complex * A, gsl_vector_complex * tau)  
{
  if (A->size1 != A->size2)
    {
      GSL_ERROR ("hermitian tridiagonal decomposition requires square matrix",
                 GSL_ENOTSQR);
    }
  else if (tau->size + 1 != A->size1)
    {
      GSL_ERROR ("size of tau must be (matrix size - 1)", GSL_EBADLEN);
    }
  else
    {
      const size_t N = A->size1;
      size_t i;
  
      const gsl_complex zero = gsl_complex_rect (0.0, 0.0);
      const gsl_complex one = gsl_complex_rect (1.0, 0.0);
      const gsl_complex neg_one = gsl_complex_rect (-1.0, 0.0);

      for (i = 0 ; i < N - 1; i++)
        {
          gsl_vector_complex_view c = gsl_matrix_complex_column (A, i);
          gsl_vector_complex_view v = gsl_vector_complex_subvector (&c.vector, i + 1, N - (i + 1));
          gsl_complex tau_i = gsl_linalg_complex_householder_transform (&v.vector);
          
          /* Apply the transformation H^T A H to the remaining columns */

          if ((i + 1) < (N - 1) 
              && !(GSL_REAL(tau_i) == 0.0 && GSL_IMAG(tau_i) == 0.0)) 
            {
              gsl_matrix_complex_view m = 
                gsl_matrix_complex_submatrix (A, i + 1, i + 1, 
                                              N - (i+1), N - (i+1));
              gsl_complex ei = gsl_vector_complex_get(&v.vector, 0);
              gsl_vector_complex_view x = gsl_vector_complex_subvector (tau, i, N-(i+1));
              gsl_vector_complex_set (&v.vector, 0, one);
              
              /* x = tau * A * v */
              gsl_blas_zhemv (CblasLower, tau_i, &m.matrix, &v.vector, zero, &x.vector);

              /* w = x - (1/2) tau * (x' * v) * v  */
              {
                gsl_complex xv, txv, alpha;
                gsl_blas_zdotc(&x.vector, &v.vector, &xv);
                txv = gsl_complex_mul(tau_i, xv);
                alpha = gsl_complex_mul_real(txv, -0.5);
                gsl_blas_zaxpy(alpha, &v.vector, &x.vector);
              }
              
              /* apply the transformation A = A - v w' - w v' */
              gsl_blas_zher2(CblasLower, neg_one, &v.vector, &x.vector, &m.matrix);

              gsl_vector_complex_set (&v.vector, 0, ei);
            }
          
          gsl_vector_complex_set (tau, i, tau_i);
        }
      
      return GSL_SUCCESS;
    }
}  


/*  Form the orthogonal matrix U from the packed QR matrix */

int
gsl_linalg_hermtd_unpack (const gsl_matrix_complex * A, 
                          const gsl_vector_complex * tau,
                          gsl_matrix_complex * U, 
                          gsl_vector * diag, 
                          gsl_vector * sdiag)
{
  if (A->size1 !=  A->size2)
    {
      GSL_ERROR ("matrix A must be sqaure", GSL_ENOTSQR);
    }
  else if (tau->size + 1 != A->size1)
    {
      GSL_ERROR ("size of tau must be (matrix size - 1)", GSL_EBADLEN);
    }
  else if (U->size1 != A->size1 || U->size2 != A->size1)
    {
      GSL_ERROR ("size of U must match size of A", GSL_EBADLEN);
    }
  else if (diag->size != A->size1)
    {
      GSL_ERROR ("size of diagonal must match size of A", GSL_EBADLEN);
    }
  else if (sdiag->size + 1 != A->size1)
    {
      GSL_ERROR ("size of subdiagonal must be (matrix size - 1)", GSL_EBADLEN);
    }
  else
    {
      const size_t N = A->size1;

      size_t i;

      /* Initialize U to the identity */

      gsl_matrix_complex_set_identity (U);

      for (i = N - 1; i-- > 0;)
        {
          gsl_complex ti = gsl_vector_complex_get (tau, i);

          gsl_vector_complex_const_view c = gsl_matrix_complex_const_column (A, i);

          gsl_vector_complex_const_view h = 
            gsl_vector_complex_const_subvector (&c.vector, i + 1, N - (i+1));

          gsl_matrix_complex_view m = 
            gsl_matrix_complex_submatrix (U, i + 1, i + 1, N-(i+1), N-(i+1));

          gsl_linalg_complex_householder_hm (ti, &h.vector, &m.matrix);
        }

      /* Copy diagonal into diag */

      for (i = 0; i < N; i++)
        {
          gsl_complex Aii = gsl_matrix_complex_get (A, i, i);
          gsl_vector_set (diag, i, GSL_REAL(Aii));
        }

      /* Copy subdiagonal into sdiag */

      for (i = 0; i < N - 1; i++)
        {
          gsl_complex Aji = gsl_matrix_complex_get (A, i+1, i);
          gsl_vector_set (sdiag, i, GSL_REAL(Aji));
        }

      return GSL_SUCCESS;
    }
}

int
gsl_linalg_hermtd_unpack_T (const gsl_matrix_complex * A, 
                            gsl_vector * diag, 
                            gsl_vector * sdiag)
{
  if (A->size1 !=  A->size2)
    {
      GSL_ERROR ("matrix A must be sqaure", GSL_ENOTSQR);
    }
  else if (diag->size != A->size1)
    {
      GSL_ERROR ("size of diagonal must match size of A", GSL_EBADLEN);
    }
  else if (sdiag->size + 1 != A->size1)
    {
      GSL_ERROR ("size of subdiagonal must be (matrix size - 1)", GSL_EBADLEN);
    }
  else
    {
      const size_t N = A->size1;

      size_t i;

      /* Copy diagonal into diag */

      for (i = 0; i < N; i++)
        {
          gsl_complex Aii = gsl_matrix_complex_get (A, i, i);
          gsl_vector_set (diag, i, GSL_REAL(Aii));
        }

      /* Copy subdiagonal into sd */

      for (i = 0; i < N - 1; i++)
        {
          gsl_complex Aji = gsl_matrix_complex_get (A, i+1, i);
          gsl_vector_set (sdiag, i, GSL_REAL(Aji));
        }

      return GSL_SUCCESS;
    }
}
