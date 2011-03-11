/* Paraboloid centered on (p[0],p[1]), with  
   scale factors (p[2],p[3]) and minimum p[4] */

double
my_f (const gsl_vector *v, void *params)
{
  double x, y;
  double *p = (double *)params;
  
  x = gsl_vector_get(v, 0);
  y = gsl_vector_get(v, 1);
 
  return p[2] * (x - p[0]) * (x - p[0]) +
           p[3] * (y - p[1]) * (y - p[1]) + p[4]; 
}

/* The gradient of f, df = (df/dx, df/dy). */
void 
my_df (const gsl_vector *v, void *params, 
       gsl_vector *df)
{
  double x, y;
  double *p = (double *)params;
  
  x = gsl_vector_get(v, 0);
  y = gsl_vector_get(v, 1);
 
  gsl_vector_set(df, 0, 2.0 * p[2] * (x - p[0]));
  gsl_vector_set(df, 1, 2.0 * p[3] * (y - p[1]));
}

/* Compute both f and df together. */
void 
my_fdf (const gsl_vector *x, void *params, 
        double *f, gsl_vector *df) 
{
  *f = my_f(x, params); 
  my_df(x, params, df);
}
