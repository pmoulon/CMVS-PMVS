int
main (void)
{
  const gsl_odeiv_step_type * T 
    = gsl_odeiv_step_rk4;

  gsl_odeiv_step * s 
    = gsl_odeiv_step_alloc (T, 2);

  double mu = 10;
  gsl_odeiv_system sys = {func, jac, 2, &mu};

  double t = 0.0, t1 = 100.0;
  double h = 1e-2;
  double y[2] = { 1.0, 0.0 }, y_err[2];
  double dydt_in[2], dydt_out[2];

  /* initialise dydt_in from system parameters */
  GSL_ODEIV_FN_EVAL(&sys, t, y, dydt_in);

  while (t < t1)
    {
      int status = gsl_odeiv_step_apply (s, t, h, 
                                         y, y_err, 
                                         dydt_in, 
                                         dydt_out, 
                                         &sys);

      if (status != GSL_SUCCESS)
          break;

      dydt_in[0] = dydt_out[0];
      dydt_in[1] = dydt_out[1];

      t += h;

      printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);
    }

  gsl_odeiv_step_free (s);
  return 0;
}
