## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libgslcdf.la

pkginclude_HEADERS= gsl_cdf.h 

INCLUDES = -I$(top_srcdir)

libgslcdf_la_SOURCES = beta.c betainv.c cauchy.c cauchyinv.c chisq.c chisqinv.c exponential.c exponentialinv.c exppow.c fdist.c fdistinv.c flat.c flatinv.c gamma.c gammainv.c gauss.c gaussinv.c gumbel1.c gumbel1inv.c gumbel2.c gumbel2inv.c laplace.c laplaceinv.c logistic.c logisticinv.c lognormal.c lognormalinv.c pareto.c paretoinv.c rayleigh.c rayleighinv.c tdist.c tdistinv.c weibull.c weibullinv.c binomial.c poisson.c geometric.c nbinomial.c pascal.c hypergeometric.c

noinst_HEADERS = beta_inc.c rat_eval.h test_auto.c error.h

TESTS = $(check_PROGRAMS)

check_PROGRAMS = test

test_SOURCES = test.c
test_LDADD = libgslcdf.la ../randist/libgslrandist.la ../rng/libgslrng.la ../specfunc/libgslspecfunc.la ../complex/libgslcomplex.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la
