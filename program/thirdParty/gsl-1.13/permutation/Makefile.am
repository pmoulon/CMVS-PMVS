noinst_LTLIBRARIES = libgslpermutation.la 

pkginclude_HEADERS = gsl_permutation.h gsl_permute.h gsl_permute_char.h gsl_permute_complex_double.h gsl_permute_complex_float.h gsl_permute_complex_long_double.h gsl_permute_double.h gsl_permute_float.h gsl_permute_int.h gsl_permute_long.h gsl_permute_long_double.h gsl_permute_short.h gsl_permute_uchar.h gsl_permute_uint.h gsl_permute_ulong.h gsl_permute_ushort.h gsl_permute_vector.h gsl_permute_vector_char.h gsl_permute_vector_complex_double.h gsl_permute_vector_complex_float.h gsl_permute_vector_complex_long_double.h gsl_permute_vector_double.h gsl_permute_vector_float.h gsl_permute_vector_int.h gsl_permute_vector_long.h gsl_permute_vector_long_double.h gsl_permute_vector_short.h gsl_permute_vector_uchar.h gsl_permute_vector_uint.h gsl_permute_vector_ulong.h gsl_permute_vector_ushort.h

INCLUDES = -I$(top_srcdir)

libgslpermutation_la_SOURCES = init.c file.c permutation.c permute.c canonical.c inline.c

noinst_HEADERS = permute_source.c

TESTS = $(check_PROGRAMS)

check_PROGRAMS = test

test_SOURCES = test.c

test_LDADD = libgslpermutation.la ../vector/libgslvector.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la

#CLEANFILES = test.txt test.dat


