#include <stdio.h>
#include <gsl/gsl_multifit.h>

int
main (int argc, char **argv)
{
  int i, n;
  double xi, yi, ei, chisq;
  gsl_matrix *X, *cov;
  gsl_vector *y, *w, *c;

  if (argc != 2)
    {
      fprintf (stderr,"usage: fit n < data\n");
      exit (-1);
    }

  n = atoi (argv[1]);

  X = gsl_matrix_alloc (n, 3);
  y = gsl_vector_alloc (n);
  w = gsl_vector_alloc (n);

  c = gsl_vector_alloc (3);
  cov = gsl_matrix_alloc (3, 3);

  for (i = 0; i < n; i++)
    {
      int count = fscanf (stdin, "%lg %lg %lg",
                          &xi, &yi, &ei);

      if (count != 3)
        {
          fprintf (stderr, "error reading file\n");
          exit (-1);
        }

      printf ("%g %g +/- %g\n", xi, yi, ei);
      
      gsl_matrix_set (X, i, 0, 1.0);
      gsl_matrix_set (X, i, 1, xi);
      gsl_matrix_set (X, i, 2, xi*xi);
      
      gsl_vector_set (y, i, yi);
      gsl_vector_set (w, i, 1.0/(ei*ei));
    }

  {
    gsl_multifit_linear_workspace * work 
      = gsl_multifit_linear_alloc (n, 3);
    gsl_multifit_wlinear (X, w, y, c, cov,
                          &chisq, work);
    gsl_multifit_linear_free (work);
  }

#define C(i) (gsl_vector_get(c,(i)))
#define COV(i,j) (gsl_matrix_get(cov,(i),(j)))

  {
    printf ("# best fit: Y = %g + %g X + %g X^2\n", 
            C(0), C(1), C(2));

    printf ("# covariance matrix:\n");
    printf ("[ %+.5e, %+.5e, %+.5e  \n",
               COV(0,0), COV(0,1), COV(0,2));
    printf ("  %+.5e, %+.5e, %+.5e  \n", 
               COV(1,0), COV(1,1), COV(1,2));
    printf ("  %+.5e, %+.5e, %+.5e ]\n", 
               COV(2,0), COV(2,1), COV(2,2));
    printf ("# chisq = %g\n", chisq);
  }

  gsl_matrix_free (X);
  gsl_vector_free (y);
  gsl_vector_free (w);
  gsl_vector_free (c);
  gsl_matrix_free (cov);

  return 0;
}
