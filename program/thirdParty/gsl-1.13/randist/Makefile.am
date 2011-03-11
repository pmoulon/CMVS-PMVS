noinst_LTLIBRARIES = libgslrandist.la

pkginclude_HEADERS= gsl_randist.h

INCLUDES = -I$(top_srcdir)

libgslrandist_la_SOURCES = bernoulli.c beta.c bigauss.c binomial.c cauchy.c chisq.c dirichlet.c discrete.c erlang.c exponential.c exppow.c fdist.c flat.c gamma.c gauss.c gausszig.c gausstail.c geometric.c gumbel.c hyperg.c laplace.c levy.c logarithmic.c logistic.c lognormal.c multinomial.c nbinomial.c pareto.c pascal.c poisson.c rayleigh.c shuffle.c sphere.c tdist.c weibull.c landau.c binomial_tpe.c

TESTS = $(check_PROGRAMS)

check_PROGRAMS = test

test_SOURCES = test.c
test_LDADD = libgslrandist.la ../rng/libgslrng.la ../specfunc/libgslspecfunc.la ../integration/libgslintegration.la ../complex/libgslcomplex.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la


