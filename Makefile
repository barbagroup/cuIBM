include Makefile.inc

.PHONY: cuibm src/solvers/libsolvers.a clean

main: cuibm

LIBS = lib/libNavierStokesSolvers.a lib/libCuspWrapped.a lib/libIO.a
EXT_LIBS = external/lib/libyaml-cpp.a
FINAL_LIB = lib/libcuIBM.a

cuibm: src/parameterDB.o $(LIBS) $(EXT_LIBS)  src/cuIBM.o src/bodies.o 
	nvcc $? -o bin/cuIBM

#lib/libcuIBM.a: $(LIBS) $(EXT_LIBS)
#	cd lib; libtool -static -o libcuIBM.a $?

lib/libNavierStokesSolvers.a: force_look
	cd src/solvers; $(MAKE) $(MFLAGS)

lib/libCuspWrapped.a: force_look
	cd src/cusp/; $(MAKE) $(MFLAGS)
  
lib/libIO.a: force_look
	cd src/io/; $(MAKE) $(MFLAGS)

external/lib/libyaml-cpp.a:
	cd external; $(MAKE) $(MFLAGS) all

clean:
	@rm -f lib/*.a bin/cuIBM src/*.o
	cd src/solvers; $(MAKE) $(MFLAGS) clean
	cd src/cusp; $(MAKE) $(MFLAGS) clean
	cd src/io; $(MAKE) $(MFLAGS) clean
	cd external; $(MAKE) $(MFLAGS) clean

force_look:
	true

.PHONY: cuibm
