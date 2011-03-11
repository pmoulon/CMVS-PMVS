## Process this file with automake to produce Makefile.in

noinst_LTLIBRARIES = libgslstatistics.la

pkginclude_HEADERS = gsl_statistics.h gsl_statistics_char.h gsl_statistics_double.h gsl_statistics_float.h gsl_statistics_int.h gsl_statistics_long.h gsl_statistics_long_double.h gsl_statistics_short.h gsl_statistics_uchar.h gsl_statistics_uint.h gsl_statistics_ulong.h gsl_statistics_ushort.h

INCLUDES = -I$(top_srcdir)

libgslstatistics_la_SOURCES =  mean.c variance.c absdev.c skew.c kurtosis.c lag1.c p_variance.c minmax.c ttest.c median.c covariance.c quantiles.c wmean.c wvariance.c wabsdev.c wskew.c wkurtosis.c

noinst_HEADERS = mean_source.c variance_source.c covariance_source.c absdev_source.c skew_source.c kurtosis_source.c lag1_source.c p_variance_source.c minmax_source.c ttest_source.c median_source.c quantiles_source.c wmean_source.c wvariance_source.c wabsdev_source.c wskew_source.c wkurtosis_source.c test_float_source.c test_int_source.c

check_PROGRAMS = test
TESTS = $(check_PROGRAMS)

test_SOURCES = test.c test_nist.c
test_LDADD = libgslstatistics.la ../sort/libgslsort.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la


