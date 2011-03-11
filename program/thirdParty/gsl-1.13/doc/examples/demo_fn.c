double
quadratic (double x, void *params)
{
  struct quadratic_params *p 
    = (struct quadratic_params *) params;

  double a = p->a;
  double b = p->b;
  double c = p->c;

  return (a * x + b) * x + c;
}

double
quadratic_deriv (double x, void *params)
{
  struct quadratic_params *p 
    = (struct quadratic_params *) params;

  double a = p->a;
  double b = p->b;
  double c = p->c;

  return 2.0 * a * x + b;
}

void
quadratic_fdf (double x, void *params, 
               double *y, double *dy)
{
  struct quadratic_params *p 
    = (struct quadratic_params *) params;

  double a = p->a;
  double b = p->b;
  double c = p->c;

  *y = (a * x + b) * x + c;
  *dy = 2.0 * a * x + b;
}
