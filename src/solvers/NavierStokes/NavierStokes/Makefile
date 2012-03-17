include ../../../../Makefile.inc

.PHONY: clean

SRCS=$(wildcard *.cu)
OBJ =$(SRCS:.cu=.o)

all: $(OBJ)
	$(AR) $(ARFLAGS) ruv ../../libNavierStokesSolvers.a $?
	$(RANLIB) ../../libNavierStokesSolvers.a

clean:
	@rm -f *.o

