include ../make.inc

.PHONY: clean

LIB  = ../libcg.a

main.exe: main.o
	$(LOADER) $(LDFLAGS) -o $@ $^ $(LIB) $(LIBS)

tryssymv.exe: tryssymv.o
	$(LOADER) $(LDFLAGS) -o $@ $^ $(LIB) $(LIBS)

clean:
	rm -f *.o *.exe

.SUFFIXES: .o .c
.c.o:
	$(CC) $(CFLAGS) -o $@ -c $<
