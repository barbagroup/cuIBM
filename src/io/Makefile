include ../../Makefile.inc

.PHONY: all objects install clean

SRC=$(wildcard *.cu)
OBJ=$(SRC:.cu=.o)

LIBNAME=libIO.a

all: objects install force_look

objects: $(OBJ)
	$(AR) $(ARFLAGS) ruv $(LIBNAME) $?
	$(RANLIB) $(LIBNAME)

install: objects
	@mv $(LIBNAME) ../../lib/$(LIBNAME)

force_look:
	true

clean:
	@rm *.o
