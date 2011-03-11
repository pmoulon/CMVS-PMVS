include ../Makefile.in

INCLUDES = -I../metisLib -I/p/include

CFLAGS = $(COPTIONS) $(OPTFLAGS) $(INCLUDES)
LD = $(CC) -L. 


OBJS = wkkm.o mlkkm.o 

.c.o:
	$(CC) $(CFLAGS) -c $*.c

../libmultilevel.a: $(OBJS)
	$(AR) $@ $(OBJS)
	$(RANLIB) $@

clean:
	rm -f *.o

realclean:
	rm -f *.o ; rm -f ../libmetis.a
