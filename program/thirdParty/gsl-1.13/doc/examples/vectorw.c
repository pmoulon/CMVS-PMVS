#include <stdio.h>
#include <gsl/gsl_vector.h>

int
main (void)
{
  int i; 
  gsl_vector * v = gsl_vector_alloc (100);
  
  for (i = 0; i < 100; i++)
    {
      gsl_vector_set (v, i, 1.23 + i);
    }

  {  
     FILE * f = fopen ("test.dat", "w");
     gsl_vector_fprintf (f, v, "%.5g");
     fclose (f);
  }

  gsl_vector_free (v);
  return 0;
}
