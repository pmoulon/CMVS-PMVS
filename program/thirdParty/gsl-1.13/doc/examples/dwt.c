#include <stdio.h>
#include <math.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_wavelet.h>

int
main (int argc, char **argv)
{
  int i, n = 256, nc = 20;
  double *data = malloc (n * sizeof (double));
  double *abscoeff = malloc (n * sizeof (double));
  size_t *p = malloc (n * sizeof (size_t));

  FILE * f;
  gsl_wavelet *w;
  gsl_wavelet_workspace *work;

  w = gsl_wavelet_alloc (gsl_wavelet_daubechies, 4);
  work = gsl_wavelet_workspace_alloc (n);

  f = fopen (argv[1], "r");
  for (i = 0; i < n; i++)
    {
      fscanf (f, "%lg", &data[i]);
    }
  fclose (f);

  gsl_wavelet_transform_forward (w, data, 1, n, work);

  for (i = 0; i < n; i++)
    {
      abscoeff[i] = fabs (data[i]);
    }
  
  gsl_sort_index (p, abscoeff, 1, n);
  
  for (i = 0; (i + nc) < n; i++)
    data[p[i]] = 0;
  
  gsl_wavelet_transform_inverse (w, data, 1, n, work);
  
  for (i = 0; i < n; i++)
    {
      printf ("%g\n", data[i]);
    }
  
  gsl_wavelet_free (w);
  gsl_wavelet_workspace_free (work);

  free (data);
  free (abscoeff);
  free (p);
  return 0;
}
