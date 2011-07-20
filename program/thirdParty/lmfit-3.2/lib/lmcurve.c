/*
 * Project:  LevenbergMarquardtLeastSquaresFitting
 *
 * File:     lmcurve.c
 *
 * Contents: Simplified wrapper for one-dimensional curve fitting,
 *           using Levenberg-Marquardt least-squares minimization.
 *
 * Usage:    see application sample demo/curve1.c
 *
 * Author:   Joachim Wuttke 2010
 * 
 * Homepage: www.messen-und-deuten.de/lmfit
 */
 

#include "lmmin.h"


typedef struct {
    const double *t;
    const double *y;
    double (*f) (double t, const double *par);
} lmcurve_data_struct;


void lmcurve_evaluate( const double *par, int m_dat, const void *data,
                       double *fvec, int *info )
{
    int i;
    for ( i = 0; i < m_dat; i++ )
	fvec[i] =
            ((lmcurve_data_struct*)data)->y[i] -
            ((lmcurve_data_struct*)data)->f(
                ((lmcurve_data_struct*)data)->t[i], par );
    // *info = *info; /* to prevent a 'unused variable' warning */
}


void lmcurve_fit( int n_par, double *par, int m_dat, 
                  const double *t, const double *y,
                  double (*f)( double t, const double *par ),
                  const lm_control_struct *control, lm_status_struct *status )
{
    lmcurve_data_struct data = { t, y, f };

    lmmin( n_par, par, m_dat, (const void*) &data,
           lmcurve_evaluate, control, status, lm_printout_std );
}

