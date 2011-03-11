#include <stdio.h>
#include <gsl/gsl_matrix.h>

int
main (void)
{
  int i, j, k = 0; 
  gsl_matrix * m = gsl_matrix_alloc (100, 100);
  gsl_matrix * a = gsl_matrix_alloc (100, 100);
  
  for (i = 0; i < 100; i++)
    for (j = 0; j < 100; j++)
      gsl_matrix_set (m, i, j, 0.23 + i + j);

  {  
     FILE * f = fopen ("test.dat", "wb");
     gsl_matrix_fwrite (f, m);
     fclose (f);
  }

  {  
     FILE * f = fopen ("test.dat", "rb");
     gsl_matrix_fread (f, a);
     fclose (f);
  }

  for (i = 0; i < 100; i++)
    for (j = 0; j < 100; j++)
      {
        double mij = gsl_matrix_get (m, i, j);
        double aij = gsl_matrix_get (a, i, j);
        if (mij != aij) k++;
      }

  gsl_matrix_free (m);
  gsl_matrix_free (a);

  printf ("differences = %d (should be zero)\n", k);
  return (k > 0);
}
