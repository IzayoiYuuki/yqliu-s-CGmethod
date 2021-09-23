LIB = ./libcg.a

all:
	cd SRC; make; cd ..
	cd TEST; make; cd ..

clean:
	cd SRC; make clean; cd ..
	cd TEST; make clean; cd ..
	$(RM) $(LIB)
