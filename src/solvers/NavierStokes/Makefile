include ../../../Makefile.inc

.PHONY: all clean kernels

SRCS=$(wildcard *.cu)
OBJS=$(SRCS:.cu=.o)
LIBNAME=libNavierStokesSolvers.a

all: objects kernels force_look

kernels:
	cd kernels; $(MAKE) $(MFLAGS)

objects: $(OBJS)
	$(AR) $(ARFLAGS) ruv ../$(LIBNAME) $?
	$(RANLIB) ../$(LIBNAME)

force_look:
	true

clean:
	cd kernels; $(MAKE) $(MFLAGS) clean
	@rm -f *.o
