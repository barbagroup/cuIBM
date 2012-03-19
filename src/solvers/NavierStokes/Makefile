include ../../../Makefile.inc

.PHONY: all clean

SRCS=$(wildcard *.cu)
OBJS=$(SRCS:.cu=.o)
LIBNAME=libNavierStokesSolvers.a

all: $(OBJS)
	$(AR) $(ARFLAGS) ruv ../$(LIBNAME) $?
	$(RANLIB) ../$(LIBNAME)

clean:
	@rm -f *.o
