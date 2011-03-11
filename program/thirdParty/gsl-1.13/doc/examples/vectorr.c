#include <stdio.h>
#include <gsl/gsl_vector.h>

int
main (void)
{
  int i; 
  gsl_vector * v = gsl_vector_alloc (10);

  {  
     FILE * f = fopen ("test.dat", "r");
     gsl_vector_fscanf (f, v);
     fclose (f);
  }

  for (i = 0; i < 10; i++)
    {
      printf ("%g\n", gsl_vector_get(v, i));
    }

  gsl_vector_free (v);
  return 0;
}
