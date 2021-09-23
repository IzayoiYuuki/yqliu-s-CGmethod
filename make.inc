DEFS      =
BLASLIB   = -lblas
LAPACKLIB =
LIBS      = $(BLASLIB)

CC        = gcc
CFLAGS    = -g --std=c89

FC        = gfortran
FFLAGS    = -g

LOADER    = $(FC)
LDFLAGS   =

ARCH      = ar
ARCHFLAGS = cr
RANLIB    = ranlib

RM        = rm -f
