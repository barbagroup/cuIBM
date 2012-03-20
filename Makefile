include Makefile.inc

.PHONY: cuibm src/solvers/libsolvers.a clean

main: cuibm

cuibm: lib/libNavierStokesSolvers.a lib/libCuspWrapped.a src/cuIBM.o src/io.o 
	nvcc $? -o bin/cuIBM

#src/solvers/libsolvers.a:  force_look
lib/libNavierStokesSolvers.a: force_look
	cd src/solvers; $(MAKE) $(MFLAGS)

lib/libCuspWrapped.a: force_look
	cd src/cusp/; $(MAKE) $(MFLAGS)

clean:
	@rm -f lib/*.a bin/cuIBM src/*.o
	cd src/solvers; $(MAKE) $(MFLAGS) clean
	cd src/cusp; $(MAKE) $(MFLAGS) clean

force_look:
	true

.PHONY: cuibm
