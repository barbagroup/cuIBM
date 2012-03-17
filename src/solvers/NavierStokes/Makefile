include ../../../Makefile.inc

.PHONY: clean NavierStokes local

SRCS=$(wildcard *.cu)
OBJS=$(SRCS:.cu=.o)
LIBNAME=libNavierStokesSolvers.a

all: NavierStokes local

local: $(OBJS)
	$(AR) $(ARFLAGS) ruv ../$(LIBNAME) $?
	$(RANLIB) ../$(LIBNAME)

NavierStokes: 
	cd NavierStokes; $(MAKE) $(MFLAGS)

FadlunEtAll: 
	cd FadlunEtAl; $(MAKE) $(MFLAGS)

clean:
	@rm -f *.o
	cd NavierStokes; $(MAKE) $(MFLAGS) clean
