#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

int
main (void)
{
  double data[] = { -1.0, 1.0, -1.0, 1.0,
                    -8.0, 4.0, -2.0, 1.0,
                    27.0, 9.0, 3.0, 1.0,
                    64.0, 16.0, 4.0, 1.0 };

  gsl_matrix_view m 
    = gsl_matrix_view_array (data, 4, 4);

  gsl_vector_complex *eval = gsl_vector_complex_alloc (4);
  gsl_matrix_complex *evec = gsl_matrix_complex_alloc (4, 4);

  gsl_eigen_nonsymmv_workspace * w = 
    gsl_eigen_nonsymmv_alloc (4);
  
  gsl_eigen_nonsymmv (&m.matrix, eval, evec, w);

  gsl_eigen_nonsymmv_free (w);

  gsl_eigen_nonsymmv_sort (eval, evec, 
                           GSL_EIGEN_SORT_ABS_DESC);
  
  {
    int i, j;

    for (i = 0; i < 4; i++)
      {
        gsl_complex eval_i 
           = gsl_vector_complex_get (eval, i);
        gsl_vector_complex_view evec_i 
           = gsl_matrix_complex_column (evec, i);

        printf ("eigenvalue = %g + %gi\n",
                GSL_REAL(eval_i), GSL_IMAG(eval_i));
        printf ("eigenvector = \n");
        for (j = 0; j < 4; ++j)
          {
            gsl_complex z = 
              gsl_vector_complex_get(&evec_i.vector, j);
            printf("%g + %gi\n", GSL_REAL(z), GSL_IMAG(z));
          }
      }
  }

  gsl_vector_complex_free(eval);
  gsl_matrix_complex_free(evec);

  return 0;
}
