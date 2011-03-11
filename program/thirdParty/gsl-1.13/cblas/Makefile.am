lib_LTLIBRARIES = libgslcblas.la
libgslcblas_la_LDFLAGS = -version-info $(GSL_LT_CBLAS_VERSION)

MINGW32_HOST = @MINGW32_HOST@
if MINGW32_HOST
libgslcblas_la_LDFLAGS += -no-undefined
endif

pkginclude_HEADERS = gsl_cblas.h

INCLUDES = -I$(top_srcdir)

libgslcblas_la_SOURCES = sasum.c saxpy.c scasum.c scnrm2.c scopy.c sdot.c sdsdot.c sgbmv.c sgemm.c sgemv.c sger.c snrm2.c srot.c srotg.c srotm.c srotmg.c ssbmv.c sscal.c sspmv.c sspr.c sspr2.c sswap.c ssymm.c ssymv.c ssyr.c ssyr2.c ssyr2k.c ssyrk.c stbmv.c stbsv.c stpmv.c stpsv.c strmm.c strmv.c strsm.c strsv.c dasum.c daxpy.c dcopy.c ddot.c dgbmv.c dgemm.c dgemv.c dger.c dnrm2.c drot.c drotg.c drotm.c drotmg.c dsbmv.c dscal.c dsdot.c dspmv.c dspr.c dspr2.c dswap.c dsymm.c dsymv.c dsyr.c dsyr2.c dsyr2k.c dsyrk.c dtbmv.c dtbsv.c dtpmv.c dtpsv.c dtrmm.c dtrmv.c dtrsm.c dtrsv.c dzasum.c dznrm2.c caxpy.c ccopy.c cdotc_sub.c cdotu_sub.c cgbmv.c cgemm.c cgemv.c cgerc.c cgeru.c chbmv.c chemm.c chemv.c cher.c cher2.c cher2k.c cherk.c chpmv.c chpr.c chpr2.c cscal.c csscal.c cswap.c csymm.c csyr2k.c csyrk.c ctbmv.c ctbsv.c ctpmv.c ctpsv.c ctrmm.c ctrmv.c ctrsm.c ctrsv.c zaxpy.c zcopy.c zdotc_sub.c zdotu_sub.c zdscal.c zgbmv.c zgemm.c zgemv.c zgerc.c zgeru.c zhbmv.c zhemm.c zhemv.c zher.c zher2.c zher2k.c zherk.c zhpmv.c zhpr.c zhpr2.c zscal.c zswap.c zsymm.c zsyr2k.c zsyrk.c ztbmv.c ztbsv.c ztpmv.c ztpsv.c ztrmm.c ztrmv.c ztrsm.c ztrsv.c icamax.c idamax.c isamax.c izamax.c xerbla.c

noinst_HEADERS = tests.c tests.h cblas.h source_asum_c.h source_asum_r.h source_axpy_c.h source_axpy_r.h source_copy_c.h source_copy_r.h source_dot_c.h source_dot_r.h source_gbmv_c.h source_gbmv_r.h source_gemm_c.h source_gemm_r.h source_gemv_c.h source_gemv_r.h source_ger.h source_gerc.h source_geru.h source_hbmv.h source_hemm.h source_hemv.h source_her.h source_her2.h source_her2k.h source_herk.h source_hpmv.h source_hpr.h source_hpr2.h source_iamax_c.h source_iamax_r.h source_nrm2_c.h source_nrm2_r.h source_rot.h source_rotg.h source_rotm.h source_rotmg.h source_sbmv.h source_scal_c.h source_scal_c_s.h source_scal_r.h source_spmv.h source_spr.h source_spr2.h source_swap_c.h source_swap_r.h source_symm_c.h source_symm_r.h source_symv.h source_syr.h source_syr2.h source_syr2k_c.h source_syr2k_r.h source_syrk_c.h source_syrk_r.h source_tbmv_c.h source_tbmv_r.h source_tbsv_c.h source_tbsv_r.h source_tpmv_c.h source_tpmv_r.h source_tpsv_c.h source_tpsv_r.h source_trmm_c.h source_trmm_r.h source_trmv_c.h source_trmv_r.h source_trsm_c.h source_trsm_r.h source_trsv_c.h source_trsv_r.h hypot.c

check_PROGRAMS = test
TESTS = $(check_PROGRAMS)

test_LDADD = libgslcblas.la ../ieee-utils/libgslieeeutils.la ../err/libgslerr.la ../test/libgsltest.la ../sys/libgslsys.la
test_SOURCES = test.c test_amax.c test_asum.c test_axpy.c test_copy.c test_dot.c test_gbmv.c test_gemm.c test_gemv.c test_ger.c test_hbmv.c test_hemm.c test_hemv.c test_her.c test_her2.c test_her2k.c test_herk.c test_hpmv.c test_hpr.c test_hpr2.c test_nrm2.c test_rot.c test_rotg.c test_rotm.c test_rotmg.c test_sbmv.c test_scal.c test_spmv.c test_spr.c test_spr2.c test_swap.c test_symm.c test_symv.c test_syr.c test_syr2.c test_syr2k.c test_syrk.c test_tbmv.c test_tbsv.c test_tpmv.c test_tpsv.c test_trmm.c test_trmv.c test_trsm.c test_trsv.c



