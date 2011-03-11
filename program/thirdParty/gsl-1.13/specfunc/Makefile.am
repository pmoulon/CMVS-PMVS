noinst_LTLIBRARIES = libgslspecfunc.la 

pkginclude_HEADERS = gsl_sf.h gsl_sf_airy.h gsl_sf_bessel.h gsl_sf_clausen.h gsl_sf_coulomb.h gsl_sf_coupling.h gsl_sf_dawson.h gsl_sf_debye.h gsl_sf_dilog.h gsl_sf_elementary.h gsl_sf_ellint.h gsl_sf_elljac.h gsl_sf_erf.h gsl_sf_exp.h gsl_sf_expint.h gsl_sf_fermi_dirac.h gsl_sf_gamma.h gsl_sf_gegenbauer.h gsl_sf_hyperg.h gsl_sf_laguerre.h gsl_sf_lambert.h gsl_sf_legendre.h gsl_sf_log.h gsl_sf_mathieu.h gsl_sf_pow_int.h gsl_sf_psi.h gsl_sf_result.h gsl_sf_synchrotron.h gsl_sf_transport.h gsl_sf_trig.h gsl_sf_zeta.h gsl_specfunc.h

noinst_HEADERS = bessel_amp_phase.h bessel_olver.h bessel_temme.h bessel.h hyperg.h legendre.h eval.h chebyshev.h cheb_eval.c cheb_eval_mode.c check.h error.h

INCLUDES = -I$(top_srcdir)

libgslspecfunc_la_SOURCES = airy.c airy_der.c airy_zero.c atanint.c bessel.c bessel.h bessel_I0.c bessel_I1.c bessel_In.c bessel_Inu.c bessel_J0.c bessel_J1.c bessel_Jn.c bessel_Jnu.c bessel_K0.c bessel_K1.c bessel_Kn.c bessel_Knu.c bessel_Y0.c bessel_Y1.c bessel_Yn.c bessel_Ynu.c bessel_amp_phase.c bessel_amp_phase.h bessel_i.c bessel_j.c bessel_k.c bessel_olver.c bessel_temme.c bessel_y.c bessel_zero.c bessel_sequence.c beta.c beta_inc.c clausen.c coulomb.c coupling.c coulomb_bound.c dawson.c debye.c dilog.c elementary.c ellint.c elljac.c erfc.c exp.c expint.c expint3.c fermi_dirac.c gegenbauer.c gamma.c gamma_inc.c hyperg_0F1.c hyperg_2F0.c hyperg_1F1.c hyperg_2F1.c hyperg_U.c hyperg.c laguerre.c lambert.c legendre_H3d.c legendre_Qn.c legendre_con.c legendre_poly.c log.c mathieu_angfunc.c mathieu_charv.c mathieu_coeff.c mathieu_radfunc.c mathieu_workspace.c poch.c pow_int.c psi.c recurse.h result.c shint.c sinint.c synchrotron.c transport.c trig.c zeta.c

TESTS = $(check_PROGRAMS)

check_PROGRAMS = test

test_LDADD = libgslspecfunc.la ../eigen/libgsleigen.la ../linalg/libgsllinalg.la  ../sort/libgslsort.la ../matrix/libgslmatrix.la ../vector/libgslvector.la ../blas/libgslblas.la ../cblas/libgslcblas.la ../block/libgslblock.la ../complex/libgslcomplex.la ../ieee-utils/libgslieeeutils.la  ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la

test_SOURCES = test_sf.c test_sf.h test_airy.c test_bessel.c test_coulomb.c test_dilog.c test_gamma.c test_hyperg.c test_legendre.c test_mathieu.c
  

