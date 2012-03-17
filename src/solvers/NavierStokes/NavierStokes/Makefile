include ../../../../Makefile.inc

SRCS=$(wildcard *.cu)
OBJ =$(SRCS:.cu=.o)

all: $(OBJ)
	$(AR) $(ARFLAGS) ruv ../../libsolvers.a $?
	$(RANLIB) ../../libsolvers.a

clean:
	-rm -f $(OBJS)

