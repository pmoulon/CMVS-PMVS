#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_histogram2d.h>

int
main (void)
{
  const gsl_rng_type * T;
  gsl_rng * r;

  gsl_histogram2d * h = gsl_histogram2d_alloc (10, 10);

  gsl_histogram2d_set_ranges_uniform (h, 
                                      0.0, 1.0,
                                      0.0, 1.0);

  gsl_histogram2d_accumulate (h, 0.3, 0.3, 1);
  gsl_histogram2d_accumulate (h, 0.8, 0.1, 5);
  gsl_histogram2d_accumulate (h, 0.7, 0.9, 0.5);

  gsl_rng_env_setup ();
  
  T = gsl_rng_default;
  r = gsl_rng_alloc (T);

  {
    int i;
    gsl_histogram2d_pdf * p 
      = gsl_histogram2d_pdf_alloc (h->nx, h->ny);
    
    gsl_histogram2d_pdf_init (p, h);

    for (i = 0; i < 1000; i++) {
      double x, y;
      double u = gsl_rng_uniform (r);
      double v = gsl_rng_uniform (r);
       
      gsl_histogram2d_pdf_sample (p, u, v, &x, &y);
      
      printf ("%g %g\n", x, y);
    }

    gsl_histogram2d_pdf_free (p);
  }

  gsl_histogram2d_free (h);
  gsl_rng_free (r);

  return 0;
}
