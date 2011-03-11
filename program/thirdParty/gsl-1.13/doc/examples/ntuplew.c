#include <gsl/gsl_ntuple.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

struct data
{
  double x;
  double y;
  double z;
};

int
main (void)
{
  const gsl_rng_type * T;
  gsl_rng * r;

  struct data ntuple_row;
  int i;

  gsl_ntuple *ntuple 
    = gsl_ntuple_create ("test.dat", &ntuple_row, 
                         sizeof (ntuple_row));

  gsl_rng_env_setup ();

  T = gsl_rng_default; 
  r = gsl_rng_alloc (T);

  for (i = 0; i < 10000; i++)
    {
      ntuple_row.x = gsl_ran_ugaussian (r);
      ntuple_row.y = gsl_ran_ugaussian (r);
      ntuple_row.z = gsl_ran_ugaussian (r);
      
      gsl_ntuple_write (ntuple);
    }
  
  gsl_ntuple_close (ntuple);
  gsl_rng_free (r);

  return 0;
}
