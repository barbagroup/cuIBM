include ../../../../Makefile.inc

.PHONY: all clean

SRCS=$(wildcard *.cu)
OBJS=$(SRCS:.cu=.o)
LIBNAME=libNavierStokesSolvers.a

all: objects force_look

objects: $(OBJS)
	$(AR) $(ARFLAGS) ruv ../../$(LIBNAME) $?
	$(RANLIB) ../../$(LIBNAME)

force_look:
	true

clean:
	@rm -f *.o
