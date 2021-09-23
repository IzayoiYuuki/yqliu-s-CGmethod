include ../make.inc

.PHONY: clean

OBJS = cgmethod.o
LIB  = ../libcg.a

lib: $(OBJS)
	$(AR) $(ARFLAGS) $(LIB) $(OBJS)
	$(RANLIB) $(LIB)

clean:
	rm -f *.o *.exe

.SUFFIXES: .o .c
.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<
