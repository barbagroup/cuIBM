.PHONY: lib
.SUFFIXES: .cpp

CC = nvcc ${OSXOPTS}
OPTS = -arch=compute_20 -O3 -I../include

SRC=$(wildcard *.cpp)
OBJ=$(SRC:.cpp=.o)

LIBNAME=libyaml-cpp.a

lib: $(OBJ)
	$(AR) ruv  ../$(LIBNAME) $?
	ranlib ../$(LIBNAME)

clean:
	@rm -f *.o

.cpp.o:
	$(CC) $(OPTS) -c $< -o $@


