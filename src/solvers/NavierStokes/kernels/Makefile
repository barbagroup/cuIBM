include ../../../../Makefile.inc

.PHONY: all clean

SRCS=$(wildcard *.cu)
OBJS=$(SRCS:.cu=.o)
LIBNAME=libNavierStokesSolvers.a
LIBTESTS=libTests.a

all: objects force_look

objects: $(OBJS)
	$(AR) $(ARFLAGS) ruv ../../$(LIBNAME) $?
	$(RANLIB) ../../$(LIBNAME)

force_look:
	true

clean:
	@rm -f *.o

tests: tests_objects force_look

tests_objects: generateA.o generateBC1.o generateBC2.o generateM.o generateQT.o generateRN.o
	$(AR) $(ARFLAGS) ruv ../../${LIBTESTS} $?
	$(RANLIB) ../../$(LIBTESTS)