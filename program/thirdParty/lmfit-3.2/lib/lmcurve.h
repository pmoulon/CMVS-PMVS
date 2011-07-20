/*
 * Project:  LevenbergMarquardtLeastSquaresFitting
 *
 * File:     lmcurve.h
 *
 * Contents: Simplified interface for one-dimensional curve fitting
 *
 * Author:   Joachim Wuttke 2010
 * 
 * Homepage: www.messen-und-deuten.de/lmfit
 */
 
#include<lmmin.h>

#ifndef LMCURVE_H
#define LMCURVE_H

#ifdef __cplusplus
extern "C" {
#endif

void lmcurve_fit( int n_par, double *par, int m_dat,
                  const double *t, const double *y,
                  double (*f)( double t, const double *par ),
                  const lm_control_struct *control, lm_status_struct *status );

#ifdef __cplusplus
}
#endif

#endif /* LMCURVE_H */
