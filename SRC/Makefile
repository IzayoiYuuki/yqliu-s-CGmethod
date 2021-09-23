CGtest:main.o cgmethod.o
	gfortran -o CGtest main.o cgmethod.o -lblas

main.o:main.c
	gcc -c main.c

cgmethod.o:cgmethod.c
	gcc -c cgmethod.c 

.PHONY:tryssymv

tryssymv:tryssymv.o
	gfortran -o tryssymv tryssymv.o -lblas

tryssymv.o:tryssymv.c
	gcc -c tryssymv.c	

.PHONY:clean

clean:
	rm -f tryssymv tryssymv.o main.o cgmethod.o CGtest
