include ../make.inc
F2CLIB       = ../F2CLIBS/libf2c.a

.SUFFIXES : .o .c 
all:  testlsame testslamch testdlamch testsecond testdsecnd testieee testversion 

testlsame:  lsame.o lsametst.o
	$(CC) $(LOADOPTS)   -o testlsame lsame.o lsametst.o $(F2CLIB) -lm

testslamch: slamch.o lsame.o slamchtst.o
	$(CC) $(LOADOPTS)   -o testslamch slamch.o lsame.o slamchtst.o $(F2CLIB) -lm

testdlamch: dlamch.o lsame.o dlamchtst.o
	$(CC) $(LOADOPTS)   -o testdlamch dlamch.o lsame.o dlamchtst.o $(F2CLIB) -lm

testsecond: second.o secondtst.o
	$(CC) $(LOADOPTS)   -o testsecond second.o secondtst.o $(F2CLIB) -lm

testdsecnd: dsecnd.o dsecndtst.o
	$(CC) $(LOADOPTS)   -o testdsecnd dsecnd.o dsecndtst.o $(F2CLIB) -lm

testieee: tstiee.o
	$(CC) $(LOADOPTS)   -o testieee tstiee.o $(F2CLIB) -lm

testversion: ilaver.o LAPACK_version.o
	    $(CC) $(LOADOPTS)   -o testversion ilaver.o LAPACK_version.o $(F2CLIB) -lm

clean:
	rm -f *.o

slamch.o: slamch.c ; $(CC) $(NOOPT) -I../INCLUDE -c $< -o $@
dlamch.o: dlamch.c ; $(CC) $(NOOPT) -I../INCLUDE -c $< -o $@

.c.o: ; $(CC) $(CFLAGS) -I../INCLUDE -c $< -o $@
