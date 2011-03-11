#include <stdio.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_ieee_utils.h>

int
main (void)
{
  double x = 1, oldsum = 0, sum = 0; 
  int i = 0;

  gsl_ieee_env_setup (); /* read GSL_IEEE_MODE */

  do 
    {
      i++;
      
      oldsum = sum;
      sum += x;
      x = x / i;
      
      printf ("i=%2d sum=%.18f error=%g\n",
              i, sum, sum - M_E);

      if (i > 30)
         break;
    }  
  while (sum != oldsum);

  return 0;
}
