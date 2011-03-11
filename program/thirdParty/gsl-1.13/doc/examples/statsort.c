#include <stdio.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

int
main(void)
{
  double data[5] = {17.2, 18.1, 16.5, 18.3, 12.6};
  double median, upperq, lowerq;

  printf ("Original dataset:  %g, %g, %g, %g, %g\n",
         data[0], data[1], data[2], data[3], data[4]);

  gsl_sort (data, 1, 5);

  printf ("Sorted dataset: %g, %g, %g, %g, %g\n",
         data[0], data[1], data[2], data[3], data[4]);

  median 
    = gsl_stats_median_from_sorted_data (data, 
                                         1, 5);

  upperq 
    = gsl_stats_quantile_from_sorted_data (data, 
                                           1, 5,
                                           0.75);
  lowerq 
    = gsl_stats_quantile_from_sorted_data (data, 
                                           1, 5,
                                           0.25);

  printf ("The median is %g\n", median);
  printf ("The upper quartile is %g\n", upperq);
  printf ("The lower quartile is %g\n", lowerq);
  return 0;
}
