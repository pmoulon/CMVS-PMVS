#######################################################################
#  This is the makefile to create a library for the BLAS.
#  The files are grouped as follows:
#
#       SBLAS1 -- Single precision real BLAS routines
#       CBLAS1 -- Single precision complex BLAS routines
#       DBLAS1 -- Double precision real BLAS routines
#       ZBLAS1 -- Double precision complex BLAS routines
#
#       CB1AUX -- Real BLAS routines called by complex routines
#       ZB1AUX -- D.P. real BLAS routines called by d.p. complex
#                 routines
#
#      ALLBLAS -- Auxiliary routines for Level 2 and 3 BLAS
#
#       SBLAS2 -- Single precision real BLAS2 routines
#       CBLAS2 -- Single precision complex BLAS2 routines
#       DBLAS2 -- Double precision real BLAS2 routines
#       ZBLAS2 -- Double precision complex BLAS2 routines
#
#       SBLAS3 -- Single precision real BLAS3 routines
#       CBLAS3 -- Single precision complex BLAS3 routines
#       DBLAS3 -- Double precision real BLAS3 routines
#       ZBLAS3 -- Double precision complex BLAS3 routines
#
#  The library can be set up to include routines for any combination
#  of the four precisions.  To create or add to the library, enter make
#  followed by one or more of the precisions desired.  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Note that these commands are not safe for parallel builds.
#
#  Alternatively, the commands
#       make all
#  or
#       make
#  without any arguments creates a library of all four precisions.
#  The name of the library is held in BLASLIB, which is set in the
#  top-level make.inc
#
#  To remove the object files after the library is created, enter
#       make clean
#  To force the source files to be recompiled, enter, for example,
#       make single FRC=FRC
#
#---------------------------------------------------------------------
#
#  Edward Anderson, University of Tennessee
#  March 26, 1990
#  Susan Ostrouchov, Last updated September 30, 1994
#  ejr, May 2006.
#
#######################################################################

#---------------------------------------------------------
#  Comment out the next 6 definitions if you already have
#  the Level 1 BLAS.
#---------------------------------------------------------
set(SBLAS1 isamax.c sasum.c saxpy.c scopy.c sdot.c snrm2.c 
	srot.c srotg.c sscal.c sswap.c sdsdot.c srotmg.c srotm.c)

set(CBLAS1 scabs1.c scasum.c scnrm2.c icamax.c caxpy.c ccopy.c 
	cdotc.c cdotu.c csscal.c crotg.c cscal.c cswap.c csrot.c)

set(DBLAS1 idamax.c dasum.c daxpy.c dcopy.c ddot.c dnrm2.c 
	drot.c drotg.c dscal.c dsdot.c dswap.c drotmg.c drotm.c)

set(ZBLAS1 dcabs1.c dzasum.c dznrm2.c izamax.c zaxpy.c zcopy.c 
	zdotc.c zdotu.c zdscal.c zrotg.c zscal.c zswap.c zdrot.c)

set(CB1AUX  isamax.c sasum.c saxpy.c scopy.c snrm2.c sscal.c)

set(ZB1AUX idamax.c dasum.c daxpy.c dcopy.c dnrm2.c dscal.c)

#---------------------------------------------------------------------
#  The following line defines auxiliary routines needed by both the
#  Level 2 and Level 3 BLAS.  Comment it out only if you already have
#  both the Level 2 and 3 BLAS.
#---------------------------------------------------------------------
set(ALLBLAS  lsame.c xerbla.c xerbla_array.c)

#---------------------------------------------------------
#  Comment out the next 4 definitions if you already have
#  the Level 2 BLAS.
#---------------------------------------------------------
set(SBLAS2 sgemv.c sgbmv.c ssymv.c ssbmv.c sspmv.c 
	strmv.c stbmv.c stpmv.c strsv.c stbsv.c stpsv.c 
	sger.c ssyr.c sspr.c ssyr2.c sspr2.c)

set(CBLAS2 cgemv.c cgbmv.c chemv.c chbmv.c chpmv.c 
	ctrmv.c ctbmv.c ctpmv.c ctrsv.c ctbsv.c ctpsv.c 
	cgerc.c cgeru.c cher.c chpr.c cher2.c chpr2.c)

set(DBLAS2 dgemv.c dgbmv.c dsymv.c dsbmv.c dspmv.c 
	dtrmv.c dtbmv.c dtpmv.c dtrsv.c dtbsv.c dtpsv.c 
	dger.c dsyr.c dspr.c dsyr2.c dspr2.c)

set(ZBLAS2 zgemv.c zgbmv.c zhemv.c zhbmv.c zhpmv.c 
	ztrmv.c ztbmv.c ztpmv.c ztrsv.c ztbsv.c ztpsv.c 
	zgerc.c zgeru.c zher.c zhpr.c zher2.c zhpr2.c)

#---------------------------------------------------------
#  Comment out the next 4 definitions if you already have
#  the Level 3 BLAS.
#---------------------------------------------------------
set(SBLAS3 sgemm.c ssymm.c ssyrk.c ssyr2k.c strmm.c strsm.c )

set(CBLAS3 cgemm.c csymm.c csyrk.c csyr2k.c ctrmm.c ctrsm.c 
	chemm.c cherk.c cher2k.c)

set(DBLAS3 dgemm.c dsymm.c dsyrk.c dsyr2k.c dtrmm.c dtrsm.c)

set(ZBLAS3 zgemm.c zsymm.c zsyrk.c zsyr2k.c ztrmm.c ztrsm.c 
	zhemm.c zherk.c zher2k.c)
# default build all of it
set(ALLOBJ ${SBLAS1} ${SBLAS2} ${SBLAS3} ${DBLAS1} ${DBLAS2} ${DBLAS3}	
	${CBLAS1} ${CBLAS2} ${CBLAS3} ${ZBLAS1} 
	${ZBLAS2} ${ZBLAS3} ${ALLBLAS})

if(BLAS_SINGLE)
  set(ALLOBJ ${SBLAS1} ${ALLBLAS} 
	${SBLAS2} ${SBLAS3})
endif()
if(BLAS_DOUBLE)
  set(ALLOBJ ${DBLAS1} ${ALLBLAS} 
	${DBLAS2} ${DBLAS3})
endif()
if(BLAS_COMPLEX)
  set(ALLOBJ  ${BLASLIB} ${CBLAS1} ${CB1AUX} 
	${ALLBLAS} ${CBLAS2})
endif()
if(BLAS_COMPLEX16)
  set(ALLOBJ  ${BLASLIB} ${ZBLAS1} ${ZB1AUX} 
	${ALLBLAS} ${ZBLAS2} ${ZBLAS3})
endif()
  
  
add_library(blas ${ALLOBJ})
if(UNIX)
  target_link_libraries(blas m)
endif()
target_link_libraries(blas f2c)
