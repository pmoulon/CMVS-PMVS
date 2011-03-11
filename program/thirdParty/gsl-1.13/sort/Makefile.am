noinst_LTLIBRARIES = libgslsort.la

pkginclude_HEADERS = gsl_heapsort.h gsl_sort.h gsl_sort_char.h gsl_sort_double.h gsl_sort_float.h gsl_sort_int.h gsl_sort_long.h gsl_sort_long_double.h gsl_sort_short.h gsl_sort_uchar.h gsl_sort_uint.h gsl_sort_ulong.h gsl_sort_ushort.h gsl_sort_vector.h gsl_sort_vector_char.h gsl_sort_vector_double.h gsl_sort_vector_float.h gsl_sort_vector_int.h gsl_sort_vector_long.h gsl_sort_vector_long_double.h gsl_sort_vector_short.h gsl_sort_vector_uchar.h gsl_sort_vector_uint.h gsl_sort_vector_ulong.h gsl_sort_vector_ushort.h

INCLUDES = -I$(top_srcdir)

libgslsort_la_SOURCES = sort.c sortind.c sortvec.c sortvecind.c subset.c subsetind.c
noinst_HEADERS = sortvec_source.c sortvecind_source.c subset_source.c subsetind_source.c test_source.c test_heapsort.c 

TESTS = $(check_PROGRAMS)

check_PROGRAMS = test

test_SOURCES = test.c
test_LDADD = libgslsort.la ../permutation/libgslpermutation.la ../vector/libgslvector.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la

