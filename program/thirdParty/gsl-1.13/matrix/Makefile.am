noinst_LTLIBRARIES = libgslmatrix.la 

check_PROGRAMS = test test_static

pkginclude_HEADERS = gsl_matrix.h gsl_matrix_char.h gsl_matrix_complex_double.h gsl_matrix_complex_float.h gsl_matrix_complex_long_double.h gsl_matrix_double.h gsl_matrix_float.h gsl_matrix_int.h gsl_matrix_long.h gsl_matrix_long_double.h gsl_matrix_short.h gsl_matrix_uchar.h gsl_matrix_uint.h gsl_matrix_ulong.h gsl_matrix_ushort.h

INCLUDES = -I$(top_srcdir)

TESTS = $(check_PROGRAMS) 

test_LDADD = libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la  ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la
test_static_LDADD = libgslmatrix.la ../vector/libgslvector.la ../block/libgslblock.la ../ieee-utils/libgslieeeutils.la  ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la ../utils/libutils.la

test_SOURCES = test.c

test_static_SOURCES = test_static.c

CLEANFILES = test.txt test.dat

noinst_HEADERS = init_source.c file_source.c rowcol_source.c swap_source.c copy_source.c test_complex_source.c test_source.c minmax_source.c prop_source.c oper_source.c getset_source.c view_source.c submatrix_source.c oper_complex_source.c

libgslmatrix_la_SOURCES = init.c matrix.c file.c rowcol.c swap.c copy.c minmax.c prop.c oper.c getset.c view.c submatrix.c view.h


