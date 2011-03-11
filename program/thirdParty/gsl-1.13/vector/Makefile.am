noinst_LTLIBRARIES = libgslvector.la 

check_PROGRAMS = test test_static

pkginclude_HEADERS = gsl_vector.h gsl_vector_char.h gsl_vector_complex.h gsl_vector_complex_double.h gsl_vector_complex_float.h gsl_vector_complex_long_double.h gsl_vector_double.h gsl_vector_float.h gsl_vector_int.h gsl_vector_long.h gsl_vector_long_double.h gsl_vector_short.h gsl_vector_uchar.h gsl_vector_uint.h gsl_vector_ulong.h gsl_vector_ushort.h

INCLUDES = -I$(top_srcdir)

TESTS = $(check_PROGRAMS)

test_LDADD = libgslvector.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la
test_static_LDADD = libgslvector.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la

test_SOURCES = test.c
test_static_SOURCES = test_static.c

CLEANFILES = test.txt test.dat

noinst_HEADERS = init_source.c file_source.c copy_source.c swap_source.c prop_source.c test_complex_source.c test_source.c minmax_source.c oper_source.c oper_complex_source.c reim_source.c subvector_source.c view_source.c

libgslvector_la_SOURCES = init.c file.c vector.c copy.c swap.c prop.c minmax.c oper.c reim.c subvector.c view.c view.h
