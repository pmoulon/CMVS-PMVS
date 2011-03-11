#include <stdio.h>
#include <gsl/gsl_matrix.h>

int
main (void)
{
  int i, j; 
  gsl_matrix * m = gsl_matrix_alloc (10, 3);
  
  for (i = 0; i < 10; i++)
    for (j = 0; j < 3; j++)
      gsl_matrix_set (m, i, j, 0.23 + 100*i + j);
  
  for (i = 0; i < 100; i++)  /* OUT OF RANGE ERROR */
    for (j = 0; j < 3; j++)
      printf ("m(%d,%d) = %g\n", i, j, 
              gsl_matrix_get (m, i, j));

  gsl_matrix_free (m);

  return 0;
}
