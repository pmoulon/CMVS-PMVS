#######################################################################
#  This is the makefile to create a library for LAPACK.
#  The files are organized as follows:
#       ALLAUX -- Auxiliary routines called from all precisions
#       ALLXAUX -- Auxiliary routines called from all precisions but
#                  only from routines using extra precision.
#       SCLAUX -- Auxiliary routines called from both REAL and COMPLEX
#       DZLAUX -- Auxiliary routines called from both DOUBLE PRECISION
#                 and COMPLEX*16
#       SLASRC -- Single precision real LAPACK routines
#       SXLASRC -- Single precision real LAPACK routines using extra
#                  precision.
#       CLASRC -- Single precision complex LAPACK routines
#       CXLASRC -- Single precision complex LAPACK routines using extra
#                  precision.
#       DLASRC -- Double precision real LAPACK routines
#       DXLASRC -- Double precision real LAPACK routines using extra
#                  precision.
#       ZLASRC -- Double precision complex LAPACK routines
#       ZXLASRC -- Double precision complex LAPACK routines using extra
#                  precision.
#
#  The library can be set up to include routines for any combination
#  of the four precisions.  To create or add to the library, enter make
#  followed by one or more of the precisions desired.  Some examples:
#       make single
#       make single complex
#       make single double complex complex16
#  Alternatively, the command
#       make
#  without any arguments creates a library of all four precisions.
#  The library is called
#       lapack.a
#  and is created at the next higher directory level.
#
#  To remove the object files after the library is created, enter
#       make clean
#  On some systems, you can force the source files to be recompiled by
#  entering (for example)
#       make single FRC=FRC
#
#  ***Note***
#  The functions lsame, second, dsecnd, slamch, and dlamch may have
#  to be installed before compiling the library.  Refer to the
#  installation guide, LAPACK Working Note 41, for instructions.
#
#######################################################################

set(ALLAUX  maxloc.c ilaenv.c ieeeck.c lsamen.c  iparmq.c	
    ilaprec.c ilatrans.c ilauplo.c iladiag.c chla_transtype.c 
    ../INSTALL/ilaver.c ../INSTALL/lsame.c) # xerbla.c xerbla_array.c

set(ALLXAUX )

set(SCLAUX  
   sbdsdc.c 
   sbdsqr.c sdisna.c slabad.c slacpy.c sladiv.c slae2.c  slaebz.c 
   slaed0.c slaed1.c slaed2.c slaed3.c slaed4.c slaed5.c slaed6.c 
   slaed7.c slaed8.c slaed9.c slaeda.c slaev2.c slagtf.c 
   slagts.c slamrg.c slanst.c 
   slapy2.c slapy3.c slarnv.c 
   slarra.c slarrb.c slarrc.c slarrd.c slarre.c slarrf.c slarrj.c 
   slarrk.c slarrr.c slaneg.c 
   slartg.c slaruv.c slas2.c  slascl.c 
   slasd0.c slasd1.c slasd2.c slasd3.c slasd4.c slasd5.c slasd6.c 
   slasd7.c slasd8.c slasda.c slasdq.c slasdt.c 
   slaset.c slasq1.c slasq2.c slasq3.c slasq4.c slasq5.c slasq6.c 
   slasr.c  slasrt.c slassq.c slasv2.c spttrf.c sstebz.c sstedc.c 
   ssteqr.c ssterf.c slaisnan.c sisnan.c 
   ../INSTALL/slamch.c ${SECOND_SRC})

set(DZLAUX  
   dbdsdc.c 
   dbdsqr.c ddisna.c dlabad.c dlacpy.c dladiv.c dlae2.c  dlaebz.c 
   dlaed0.c dlaed1.c dlaed2.c dlaed3.c dlaed4.c dlaed5.c dlaed6.c 
   dlaed7.c dlaed8.c dlaed9.c dlaeda.c dlaev2.c dlagtf.c 
   dlagts.c dlamrg.c dlanst.c 
   dlapy2.c dlapy3.c dlarnv.c 
   dlarra.c dlarrb.c dlarrc.c dlarrd.c dlarre.c dlarrf.c dlarrj.c 
   dlarrk.c dlarrr.c dlaneg.c 
   dlartg.c dlaruv.c dlas2.c  dlascl.c 
   dlasd0.c dlasd1.c dlasd2.c dlasd3.c dlasd4.c dlasd5.c dlasd6.c 
   dlasd7.c dlasd8.c dlasda.c dlasdq.c dlasdt.c 
   dlaset.c dlasq1.c dlasq2.c dlasq3.c dlasq4.c dlasq5.c dlasq6.c 
   dlasr.c  dlasrt.c dlassq.c dlasv2.c dpttrf.c dstebz.c dstedc.c 
   dsteqr.c dsterf.c dlaisnan.c disnan.c 
   ../INSTALL/dlamch.c ${DSECOND_SRC})

set(SLASRC  
   sgbbrd.c sgbcon.c sgbequ.c sgbrfs.c sgbsv.c  
   sgbsvx.c sgbtf2.c sgbtrf.c sgbtrs.c sgebak.c sgebal.c sgebd2.c 
   sgebrd.c sgecon.c sgeequ.c sgees.c  sgeesx.c sgeev.c  sgeevx.c 
   sgegs.c  sgegv.c  sgehd2.c sgehrd.c sgelq2.c sgelqf.c 
   sgels.c  sgelsd.c sgelss.c sgelsx.c sgelsy.c sgeql2.c sgeqlf.c 
   sgeqp3.c sgeqpf.c sgeqr2.c sgeqrf.c sgerfs.c sgerq2.c sgerqf.c 
   sgesc2.c sgesdd.c sgesv.c  sgesvd.c sgesvx.c sgetc2.c sgetf2.c 
   sgetrf.c sgetri.c 
   sgetrs.c sggbak.c sggbal.c sgges.c  sggesx.c sggev.c  sggevx.c 
   sggglm.c sgghrd.c sgglse.c sggqrf.c 
   sggrqf.c sggsvd.c sggsvp.c sgtcon.c sgtrfs.c sgtsv.c  
   sgtsvx.c sgttrf.c sgttrs.c sgtts2.c shgeqz.c 
   shsein.c shseqr.c slabrd.c slacon.c slacn2.c 
   slaein.c slaexc.c slag2.c  slags2.c slagtm.c slagv2.c slahqr.c 
   slahrd.c slahr2.c slaic1.c slaln2.c slals0.c slalsa.c slalsd.c 
   slangb.c slange.c slangt.c slanhs.c slansb.c slansp.c 
   slansy.c slantb.c slantp.c slantr.c slanv2.c 
   slapll.c slapmt.c 
   slaqgb.c slaqge.c slaqp2.c slaqps.c slaqsb.c slaqsp.c slaqsy.c 
   slaqr0.c slaqr1.c slaqr2.c slaqr3.c slaqr4.c slaqr5.c 
   slaqtr.c slar1v.c slar2v.c ilaslr.c ilaslc.c 
   slarf.c  slarfb.c slarfg.c slarft.c slarfx.c slargv.c 
   slarrv.c slartv.c slarfp.c 
   slarz.c  slarzb.c slarzt.c slaswp.c slasy2.c slasyf.c 
   slatbs.c slatdf.c slatps.c slatrd.c slatrs.c slatrz.c slatzm.c 
   slauu2.c slauum.c sopgtr.c sopmtr.c sorg2l.c sorg2r.c 
   sorgbr.c sorghr.c sorgl2.c sorglq.c sorgql.c sorgqr.c sorgr2.c 
   sorgrq.c sorgtr.c sorm2l.c sorm2r.c 
   sormbr.c sormhr.c sorml2.c sormlq.c sormql.c sormqr.c sormr2.c 
   sormr3.c sormrq.c sormrz.c sormtr.c spbcon.c spbequ.c spbrfs.c 
   spbstf.c spbsv.c  spbsvx.c 
   spbtf2.c spbtrf.c spbtrs.c spocon.c spoequ.c sporfs.c sposv.c  
   sposvx.c spotf2.c spotrf.c spotri.c spotrs.c spstrf.c spstf2.c 
   sppcon.c sppequ.c 
   spprfs.c sppsv.c  sppsvx.c spptrf.c spptri.c spptrs.c sptcon.c 
   spteqr.c sptrfs.c sptsv.c  sptsvx.c spttrs.c sptts2.c srscl.c  
   ssbev.c  ssbevd.c ssbevx.c ssbgst.c ssbgv.c  ssbgvd.c ssbgvx.c 
   ssbtrd.c sspcon.c sspev.c  sspevd.c sspevx.c sspgst.c 
   sspgv.c  sspgvd.c sspgvx.c ssprfs.c sspsv.c  sspsvx.c ssptrd.c 
   ssptrf.c ssptri.c ssptrs.c sstegr.c sstein.c sstev.c  sstevd.c sstevr.c 
   sstevx.c ssycon.c ssyev.c  ssyevd.c ssyevr.c ssyevx.c ssygs2.c 
   ssygst.c ssygv.c  ssygvd.c ssygvx.c ssyrfs.c ssysv.c  ssysvx.c 
   ssytd2.c ssytf2.c ssytrd.c ssytrf.c ssytri.c ssytrs.c stbcon.c 
   stbrfs.c stbtrs.c stgevc.c stgex2.c stgexc.c stgsen.c 
   stgsja.c stgsna.c stgsy2.c stgsyl.c stpcon.c stprfs.c stptri.c 
   stptrs.c 
   strcon.c strevc.c strexc.c strrfs.c strsen.c strsna.c strsyl.c 
   strti2.c strtri.c strtrs.c stzrqf.c stzrzf.c sstemr.c 
   slansf.c spftrf.c spftri.c spftrs.c ssfrk.c stfsm.c stftri.c stfttp.c 
   stfttr.c stpttf.c stpttr.c strttf.c strttp.c 
   sgejsv.c  sgesvj.c  sgsvj0.c  sgsvj1.c 
   sgeequb.c ssyequb.c spoequb.c sgbequb.c)

set(SXLASRC  sgesvxx.c sgerfsx.c sla_gerfsx_extended.c sla_geamv.c		
   sla_gercond.c sla_rpvgrw.c ssysvxx.c ssyrfsx.c			
   sla_syrfsx_extended.c sla_syamv.c sla_syrcond.c sla_syrpvgrw.c	
   sposvxx.c sporfsx.c sla_porfsx_extended.c sla_porcond.c		
   sla_porpvgrw.c sgbsvxx.c sgbrfsx.c sla_gbrfsx_extended.c		
   sla_gbamv.c sla_gbrcond.c sla_gbrpvgrw.c sla_lin_berr.c slarscl2.c	
   slascl2.c sla_wwaddw.c)

set(CLASRC 
   cbdsqr.c cgbbrd.c cgbcon.c cgbequ.c cgbrfs.c cgbsv.c  cgbsvx.c 
   cgbtf2.c cgbtrf.c cgbtrs.c cgebak.c cgebal.c cgebd2.c cgebrd.c 
   cgecon.c cgeequ.c cgees.c  cgeesx.c cgeev.c  cgeevx.c 
   cgegs.c  cgegv.c  cgehd2.c cgehrd.c cgelq2.c cgelqf.c 
   cgels.c  cgelsd.c cgelss.c cgelsx.c cgelsy.c cgeql2.c cgeqlf.c cgeqp3.c 
   cgeqpf.c cgeqr2.c cgeqrf.c cgerfs.c cgerq2.c cgerqf.c 
   cgesc2.c cgesdd.c cgesv.c  cgesvd.c cgesvx.c cgetc2.c cgetf2.c cgetrf.c 
   cgetri.c cgetrs.c 
   cggbak.c cggbal.c cgges.c  cggesx.c cggev.c  cggevx.c cggglm.c 
   cgghrd.c cgglse.c cggqrf.c cggrqf.c 
   cggsvd.c cggsvp.c 
   cgtcon.c cgtrfs.c cgtsv.c  cgtsvx.c cgttrf.c cgttrs.c cgtts2.c chbev.c  
   chbevd.c chbevx.c chbgst.c chbgv.c  chbgvd.c chbgvx.c chbtrd.c 
   checon.c cheev.c  cheevd.c cheevr.c cheevx.c chegs2.c chegst.c 
   chegv.c  chegvd.c chegvx.c cherfs.c chesv.c  chesvx.c chetd2.c 
   chetf2.c chetrd.c 
   chetrf.c chetri.c chetrs.c chgeqz.c chpcon.c chpev.c  chpevd.c 
   chpevx.c chpgst.c chpgv.c  chpgvd.c chpgvx.c chprfs.c chpsv.c  
   chpsvx.c 
   chptrd.c chptrf.c chptri.c chptrs.c chsein.c chseqr.c clabrd.c 
   clacgv.c clacon.c clacn2.c clacp2.c clacpy.c clacrm.c clacrt.c cladiv.c 
   claed0.c claed7.c claed8.c 
   claein.c claesy.c claev2.c clags2.c clagtm.c 
   clahef.c clahqr.c 
   clahrd.c clahr2.c claic1.c clals0.c clalsa.c clalsd.c clangb.c clange.c clangt.c 
   clanhb.c clanhe.c 
   clanhp.c clanhs.c clanht.c clansb.c clansp.c clansy.c clantb.c 
   clantp.c clantr.c clapll.c clapmt.c clarcm.c claqgb.c claqge.c 
   claqhb.c claqhe.c claqhp.c claqp2.c claqps.c claqsb.c 
   claqr0.c claqr1.c claqr2.c claqr3.c claqr4.c claqr5.c 
   claqsp.c claqsy.c clar1v.c clar2v.c ilaclr.c ilaclc.c 
   clarf.c  clarfb.c clarfg.c clarft.c clarfp.c 
   clarfx.c clargv.c clarnv.c clarrv.c clartg.c clartv.c 
   clarz.c  clarzb.c clarzt.c clascl.c claset.c clasr.c  classq.c 
   claswp.c clasyf.c clatbs.c clatdf.c clatps.c clatrd.c clatrs.c clatrz.c 
   clatzm.c clauu2.c clauum.c cpbcon.c cpbequ.c cpbrfs.c cpbstf.c cpbsv.c  
   cpbsvx.c cpbtf2.c cpbtrf.c cpbtrs.c cpocon.c cpoequ.c cporfs.c 
   cposv.c  cposvx.c cpotf2.c cpotrf.c cpotri.c cpotrs.c cpstrf.c cpstf2.c 
   cppcon.c cppequ.c cpprfs.c cppsv.c  cppsvx.c cpptrf.c cpptri.c cpptrs.c 
   cptcon.c cpteqr.c cptrfs.c cptsv.c  cptsvx.c cpttrf.c cpttrs.c cptts2.c 
   crot.c   cspcon.c cspmv.c  cspr.c   csprfs.c cspsv.c  
   cspsvx.c csptrf.c csptri.c csptrs.c csrscl.c cstedc.c 
   cstegr.c cstein.c csteqr.c csycon.c csymv.c  
   csyr.c   csyrfs.c csysv.c  csysvx.c csytf2.c csytrf.c csytri.c 
   csytrs.c ctbcon.c ctbrfs.c ctbtrs.c ctgevc.c ctgex2.c 
   ctgexc.c ctgsen.c ctgsja.c ctgsna.c ctgsy2.c ctgsyl.c ctpcon.c 
   ctprfs.c ctptri.c 
   ctptrs.c ctrcon.c ctrevc.c ctrexc.c ctrrfs.c ctrsen.c ctrsna.c 
   ctrsyl.c ctrti2.c ctrtri.c ctrtrs.c ctzrqf.c ctzrzf.c cung2l.c cung2r.c 
   cungbr.c cunghr.c cungl2.c cunglq.c cungql.c cungqr.c cungr2.c 
   cungrq.c cungtr.c cunm2l.c cunm2r.c cunmbr.c cunmhr.c cunml2.c 
   cunmlq.c cunmql.c cunmqr.c cunmr2.c cunmr3.c cunmrq.c cunmrz.c 
   cunmtr.c cupgtr.c cupmtr.c icmax1.c scsum1.c cstemr.c 
   chfrk.c ctfttp.c clanhf.c cpftrf.c cpftri.c cpftrs.c ctfsm.c ctftri.c 
   ctfttr.c ctpttf.c ctpttr.c ctrttf.c ctrttp.c 
   cgeequb.c cgbequb.c csyequb.c cpoequb.c cheequb.c)

set(CXLASRC     cgesvxx.c cgerfsx.c cla_gerfsx_extended.c cla_geamv.c 
   cla_gercond_c.c cla_gercond_x.c cla_rpvgrw.c 
   csysvxx.c csyrfsx.c cla_syrfsx_extended.c cla_syamv.c 
   cla_syrcond_c.c cla_syrcond_x.c cla_syrpvgrw.c 
   cposvxx.c cporfsx.c cla_porfsx_extended.c 
   cla_porcond_c.c cla_porcond_x.c cla_porpvgrw.c 
   cgbsvxx.c cgbrfsx.c cla_gbrfsx_extended.c cla_gbamv.c 
   cla_gbrcond_c.c cla_gbrcond_x.c cla_gbrpvgrw.c 
   chesvxx.c cherfsx.c cla_herfsx_extended.c cla_heamv.c 
   cla_hercond_c.c cla_hercond_x.c cla_herpvgrw.c 
   cla_lin_berr.c clarscl2.c clascl2.c cla_wwaddw.c)

set(DLASRC 
   dgbbrd.c dgbcon.c dgbequ.c dgbrfs.c dgbsv.c  
   dgbsvx.c dgbtf2.c dgbtrf.c dgbtrs.c dgebak.c dgebal.c dgebd2.c 
   dgebrd.c dgecon.c dgeequ.c dgees.c  dgeesx.c dgeev.c  dgeevx.c 
   dgegs.c  dgegv.c  dgehd2.c dgehrd.c dgelq2.c dgelqf.c 
   dgels.c  dgelsd.c dgelss.c dgelsx.c dgelsy.c dgeql2.c dgeqlf.c 
   dgeqp3.c dgeqpf.c dgeqr2.c dgeqrf.c dgerfs.c dgerq2.c dgerqf.c 
   dgesc2.c dgesdd.c dgesv.c  dgesvd.c dgesvx.c dgetc2.c dgetf2.c 
   dgetrf.c dgetri.c 
   dgetrs.c dggbak.c dggbal.c dgges.c  dggesx.c dggev.c  dggevx.c 
   dggglm.c dgghrd.c dgglse.c dggqrf.c 
   dggrqf.c dggsvd.c dggsvp.c dgtcon.c dgtrfs.c dgtsv.c  
   dgtsvx.c dgttrf.c dgttrs.c dgtts2.c dhgeqz.c 
   dhsein.c dhseqr.c dlabrd.c dlacon.c dlacn2.c 
   dlaein.c dlaexc.c dlag2.c  dlags2.c dlagtm.c dlagv2.c dlahqr.c 
   dlahrd.c dlahr2.c dlaic1.c dlaln2.c dlals0.c dlalsa.c dlalsd.c 
   dlangb.c dlange.c dlangt.c dlanhs.c dlansb.c dlansp.c 
   dlansy.c dlantb.c dlantp.c dlantr.c dlanv2.c 
   dlapll.c dlapmt.c 
   dlaqgb.c dlaqge.c dlaqp2.c dlaqps.c dlaqsb.c dlaqsp.c dlaqsy.c 
   dlaqr0.c dlaqr1.c dlaqr2.c dlaqr3.c dlaqr4.c dlaqr5.c 
   dlaqtr.c dlar1v.c dlar2v.c iladlr.c iladlc.c 
   dlarf.c  dlarfb.c dlarfg.c dlarft.c dlarfx.c dlargv.c 
   dlarrv.c dlartv.c dlarfp.c 
   dlarz.c  dlarzb.c dlarzt.c dlaswp.c dlasy2.c dlasyf.c 
   dlatbs.c dlatdf.c dlatps.c dlatrd.c dlatrs.c dlatrz.c dlatzm.c dlauu2.c 
   dlauum.c dopgtr.c dopmtr.c dorg2l.c dorg2r.c 
   dorgbr.c dorghr.c dorgl2.c dorglq.c dorgql.c dorgqr.c dorgr2.c 
   dorgrq.c dorgtr.c dorm2l.c dorm2r.c 
   dormbr.c dormhr.c dorml2.c dormlq.c dormql.c dormqr.c dormr2.c 
   dormr3.c dormrq.c dormrz.c dormtr.c dpbcon.c dpbequ.c dpbrfs.c 
   dpbstf.c dpbsv.c  dpbsvx.c 
   dpbtf2.c dpbtrf.c dpbtrs.c dpocon.c dpoequ.c dporfs.c dposv.c  
   dposvx.c dpotf2.c dpotrf.c dpotri.c dpotrs.c dpstrf.c dpstf2.c 
   dppcon.c dppequ.c 
   dpprfs.c dppsv.c  dppsvx.c dpptrf.c dpptri.c dpptrs.c dptcon.c 
   dpteqr.c dptrfs.c dptsv.c  dptsvx.c dpttrs.c dptts2.c drscl.c  
   dsbev.c  dsbevd.c dsbevx.c dsbgst.c dsbgv.c  dsbgvd.c dsbgvx.c 
   dsbtrd.c  dspcon.c dspev.c  dspevd.c dspevx.c dspgst.c 
   dspgv.c  dspgvd.c dspgvx.c dsprfs.c dspsv.c  dspsvx.c dsptrd.c 
   dsptrf.c dsptri.c dsptrs.c dstegr.c dstein.c dstev.c  dstevd.c dstevr.c 
   dstevx.c dsycon.c dsyev.c  dsyevd.c dsyevr.c 
   dsyevx.c dsygs2.c dsygst.c dsygv.c  dsygvd.c dsygvx.c dsyrfs.c 
   dsysv.c  dsysvx.c 
   dsytd2.c dsytf2.c dsytrd.c dsytrf.c dsytri.c dsytrs.c dtbcon.c 
   dtbrfs.c dtbtrs.c dtgevc.c dtgex2.c dtgexc.c dtgsen.c 
   dtgsja.c dtgsna.c dtgsy2.c dtgsyl.c dtpcon.c dtprfs.c dtptri.c 
   dtptrs.c 
   dtrcon.c dtrevc.c dtrexc.c dtrrfs.c dtrsen.c dtrsna.c dtrsyl.c 
   dtrti2.c dtrtri.c dtrtrs.c dtzrqf.c dtzrzf.c dstemr.c 
   dsgesv.c dsposv.c dlag2s.c slag2d.c dlat2s.c 
   dlansf.c dpftrf.c dpftri.c dpftrs.c dsfrk.c dtfsm.c dtftri.c dtfttp.c 
   dtfttr.c dtpttf.c dtpttr.c dtrttf.c dtrttp.c 
   dgejsv.c  dgesvj.c  dgsvj0.c  dgsvj1.c 
   dgeequb.c dsyequb.c dpoequb.c dgbequb.c)

set(DXLASRC dgesvxx.c dgerfsx.c dla_gerfsx_extended.c dla_geamv.c		
   dla_gercond.c dla_rpvgrw.c dsysvxx.c dsyrfsx.c			
   dla_syrfsx_extended.c dla_syamv.c dla_syrcond.c dla_syrpvgrw.c	
   dposvxx.c dporfsx.c dla_porfsx_extended.c dla_porcond.c		
   dla_porpvgrw.c dgbsvxx.c dgbrfsx.c dla_gbrfsx_extended.c		
   dla_gbamv.c dla_gbrcond.c dla_gbrpvgrw.c dla_lin_berr.c dlarscl2.c	
   dlascl2.c dla_wwaddw.c)

set(ZLASRC 
   zbdsqr.c zgbbrd.c zgbcon.c zgbequ.c zgbrfs.c zgbsv.c  zgbsvx.c 
   zgbtf2.c zgbtrf.c zgbtrs.c zgebak.c zgebal.c zgebd2.c zgebrd.c 
   zgecon.c zgeequ.c zgees.c  zgeesx.c zgeev.c  zgeevx.c 
   zgegs.c  zgegv.c  zgehd2.c zgehrd.c zgelq2.c zgelqf.c 
   zgels.c  zgelsd.c zgelss.c zgelsx.c zgelsy.c zgeql2.c zgeqlf.c zgeqp3.c 
   zgeqpf.c zgeqr2.c zgeqrf.c zgerfs.c zgerq2.c zgerqf.c 
   zgesc2.c zgesdd.c zgesv.c  zgesvd.c zgesvx.c zgetc2.c zgetf2.c zgetrf.c 
   zgetri.c zgetrs.c 
   zggbak.c zggbal.c zgges.c  zggesx.c zggev.c  zggevx.c zggglm.c 
   zgghrd.c zgglse.c zggqrf.c zggrqf.c 
   zggsvd.c zggsvp.c 
   zgtcon.c zgtrfs.c zgtsv.c  zgtsvx.c zgttrf.c zgttrs.c zgtts2.c zhbev.c  
   zhbevd.c zhbevx.c zhbgst.c zhbgv.c  zhbgvd.c zhbgvx.c zhbtrd.c 
   zhecon.c zheev.c  zheevd.c zheevr.c zheevx.c zhegs2.c zhegst.c 
   zhegv.c  zhegvd.c zhegvx.c zherfs.c zhesv.c  zhesvx.c zhetd2.c 
   zhetf2.c zhetrd.c 
   zhetrf.c zhetri.c zhetrs.c zhgeqz.c zhpcon.c zhpev.c  zhpevd.c 
   zhpevx.c zhpgst.c zhpgv.c  zhpgvd.c zhpgvx.c zhprfs.c zhpsv.c  
   zhpsvx.c 
   zhptrd.c zhptrf.c zhptri.c zhptrs.c zhsein.c zhseqr.c zlabrd.c 
   zlacgv.c zlacon.c zlacn2.c zlacp2.c zlacpy.c zlacrm.c zlacrt.c zladiv.c 
   zlaed0.c zlaed7.c zlaed8.c 
   zlaein.c zlaesy.c zlaev2.c zlags2.c zlagtm.c 
   zlahef.c zlahqr.c 
   zlahrd.c zlahr2.c zlaic1.c zlals0.c zlalsa.c zlalsd.c zlangb.c zlange.c 
   zlangt.c zlanhb.c 
   zlanhe.c 
   zlanhp.c zlanhs.c zlanht.c zlansb.c zlansp.c zlansy.c zlantb.c 
   zlantp.c zlantr.c zlapll.c zlapmt.c zlaqgb.c zlaqge.c 
   zlaqhb.c zlaqhe.c zlaqhp.c zlaqp2.c zlaqps.c zlaqsb.c 
   zlaqr0.c zlaqr1.c zlaqr2.c zlaqr3.c zlaqr4.c zlaqr5.c 
   zlaqsp.c zlaqsy.c zlar1v.c zlar2v.c ilazlr.c ilazlc.c 
   zlarcm.c zlarf.c  zlarfb.c 
   zlarfg.c zlarft.c zlarfp.c 
   zlarfx.c zlargv.c zlarnv.c zlarrv.c zlartg.c zlartv.c 
   zlarz.c  zlarzb.c zlarzt.c zlascl.c zlaset.c zlasr.c  
   zlassq.c zlaswp.c zlasyf.c 
   zlatbs.c zlatdf.c zlatps.c zlatrd.c zlatrs.c zlatrz.c zlatzm.c zlauu2.c 
   zlauum.c zpbcon.c zpbequ.c zpbrfs.c zpbstf.c zpbsv.c  
   zpbsvx.c zpbtf2.c zpbtrf.c zpbtrs.c zpocon.c zpoequ.c zporfs.c 
   zposv.c  zposvx.c zpotf2.c zpotrf.c zpotri.c zpotrs.c zpstrf.c zpstf2.c 
   zppcon.c zppequ.c zpprfs.c zppsv.c  zppsvx.c zpptrf.c zpptri.c zpptrs.c 
   zptcon.c zpteqr.c zptrfs.c zptsv.c  zptsvx.c zpttrf.c zpttrs.c zptts2.c 
   zrot.c   zspcon.c zspmv.c  zspr.c   zsprfs.c zspsv.c  
   zspsvx.c zsptrf.c zsptri.c zsptrs.c zdrscl.c zstedc.c 
   zstegr.c zstein.c zsteqr.c zsycon.c zsymv.c  
   zsyr.c   zsyrfs.c zsysv.c  zsysvx.c zsytf2.c zsytrf.c zsytri.c 
   zsytrs.c ztbcon.c ztbrfs.c ztbtrs.c ztgevc.c ztgex2.c 
   ztgexc.c ztgsen.c ztgsja.c ztgsna.c ztgsy2.c ztgsyl.c ztpcon.c 
   ztprfs.c ztptri.c 
   ztptrs.c ztrcon.c ztrevc.c ztrexc.c ztrrfs.c ztrsen.c ztrsna.c 
   ztrsyl.c ztrti2.c ztrtri.c ztrtrs.c ztzrqf.c ztzrzf.c zung2l.c 
   zung2r.c zungbr.c zunghr.c zungl2.c zunglq.c zungql.c zungqr.c zungr2.c 
   zungrq.c zungtr.c zunm2l.c zunm2r.c zunmbr.c zunmhr.c zunml2.c 
   zunmlq.c zunmql.c zunmqr.c zunmr2.c zunmr3.c zunmrq.c zunmrz.c 
   zunmtr.c zupgtr.c 
   zupmtr.c izmax1.c dzsum1.c zstemr.c 
   zcgesv.c zcposv.c zlag2c.c clag2z.c zlat2c.c 
   zhfrk.c ztfttp.c zlanhf.c zpftrf.c zpftri.c zpftrs.c ztfsm.c ztftri.c 
   ztfttr.c ztpttf.c ztpttr.c ztrttf.c ztrttp.c 
   zgeequb.c zgbequb.c zsyequb.c zpoequb.c zheequb.c)

set(ZXLASRC  zgesvxx.c zgerfsx.c zla_gerfsx_extended.c zla_geamv.c		
   zla_gercond_c.c zla_gercond_x.c zla_rpvgrw.c zsysvxx.c zsyrfsx.c	
   zla_syrfsx_extended.c zla_syamv.c zla_syrcond_c.c zla_syrcond_x.c	
   zla_syrpvgrw.c zposvxx.c zporfsx.c zla_porfsx_extended.c		
   zla_porcond_c.c zla_porcond_x.c zla_porpvgrw.c zgbsvxx.c zgbrfsx.c	
   zla_gbrfsx_extended.c zla_gbamv.c zla_gbrcond_c.c zla_gbrcond_x.c	
   zla_gbrpvgrw.c zhesvxx.c zherfsx.c zla_herfsx_extended.c		
   zla_heamv.c zla_hercond_c.c zla_hercond_x.c zla_herpvgrw.c		
   zla_lin_berr.c zlarscl2.c zlascl2.c zla_wwaddw.c)


if( USEXBLAS)
  set(ALLXOBJ ${SXLASRC} ${DXLASRC} ${CXLASRC} ${ZXLASRC} ${ALLXAUX})
endif()

set(ALLOBJ ${SLASRC} ${DLASRC} ${CLASRC} ${ZLASRC} ${SCLAUX} ${DZLAUX}	
	${ALLAUX})
if(BUILD_SINGLE)
set(ALLOBJ  ${SLASRC} ${ALLAUX} 
	${SCLAUX})
endif()
if(BUILD_DOUBLE)
  set(ALLOBJ ${DLASRC} ${ALLAUX} ${DZLAUX})
endif()
if(BUILD_COMPLEX)
  set(ALLOBJ ${CLASRC} ${ALLAUX} ${SCLAUX})
endif()
if(BUILD_COMPLEX16)
  set(ALLOBJ  ${ZLASRC} ${ALLAUX} ${DZLAUX})
endif()
add_library(lapack ${ALLOBJ} ${ALLXOBJ})
target_link_libraries(lapack blas)

