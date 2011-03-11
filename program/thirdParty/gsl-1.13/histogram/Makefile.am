noinst_LTLIBRARIES = libgslhistogram.la 

pkginclude_HEADERS = gsl_histogram.h gsl_histogram2d.h

INCLUDES = -I$(top_srcdir)

libgslhistogram_la_SOURCES = add.c  get.c init.c params.c reset.c file.c pdf.c gsl_histogram.h add2d.c get2d.c init2d.c params2d.c reset2d.c file2d.c pdf2d.c gsl_histogram2d.h calloc_range.c calloc_range2d.c copy.c copy2d.c maxval.c maxval2d.c oper.c oper2d.c stat.c stat2d.c

noinst_HEADERS = urand.c find.c find2d.c

check_PROGRAMS = test
TESTS = $(check_PROGRAMS)

EXTRA_DIST = urand.c

test_SOURCES = test.c test1d.c test2d.c test1d_resample.c test2d_resample.c test1d_trap.c test2d_trap.c
test_LDADD = libgslhistogram.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la

CLEANFILES = test.txt test.dat

