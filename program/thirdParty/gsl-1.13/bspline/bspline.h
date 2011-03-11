static inline size_t
bspline_find_interval (const double x, int *flag, gsl_bspline_workspace * w);

static inline int
bspline_process_interval_for_eval (const double x, size_t * i, const int flag,
				   gsl_bspline_workspace * w);

static void
bspline_pppack_bsplvb (const gsl_vector * t,
		       const size_t jhigh,
		       const size_t index,
		       const double x,
		       const size_t left,
		       size_t * j,
		       gsl_vector * deltal,
		       gsl_vector * deltar, gsl_vector * biatx);

static void
bspline_pppack_bsplvd (const gsl_vector * t,
		       const size_t k,
		       const double x,
		       const size_t left,
		       gsl_vector * deltal,
		       gsl_vector * deltar,
		       gsl_matrix * a,
		       gsl_matrix * dbiatx, const size_t nderiv);
