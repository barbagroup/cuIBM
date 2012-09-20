include ../../../Makefile.inc

.PHONY: all clean kernels

# list of all files with extension .cu
SRCS=$(wildcard *.cu)

# takes the list of files in SRCS, replaces .cu with .o
OBJS=$(SRCS:.cu=.o)
LIBNAME=libNavierStokesSolvers.a

# all is the default target that is always run
all: objects kernels force_look

kernels:
	cd kernels; $(MAKE) $(MFLAGS)

# all the prerequisites are files with extension .o
# the implicit rule .cu.o defined in Makefile.inc is used to compile these
# then ar is used to create a library from the object files
# $^ lists all dependencies ($(OBJS) in this case)
# $? lists only those dependencies that have been updated
objects: $(OBJS)
	$(AR) $(ARFLAGS) ruv ../$(LIBNAME) $?
	$(RANLIB) ../$(LIBNAME)

force_look:
	true

clean:
	cd kernels; $(MAKE) $(MFLAGS) clean
	@rm -f *.o
