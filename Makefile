include Makefile.inc

.PHONY: cuibm src/solvers/libsolvers.a clean

main: cuibm

cuibm: lib/libNavierStokesSolvers.a lib/libCuspWrapped.a src/cuIBM.o src/io.o src/bodies.o
	nvcc $? -o bin/cuIBM

#src/solvers/libsolvers.a:  force_look
lib/libNavierStokesSolvers.a: force_look
	cd src/solvers; $(MAKE) $(MFLAGS)

lib/libCuspWrapped.a: force_look
	cd src/cusp/; $(MAKE) $(MFLAGS)

external:
	cd external; $(MAKE) $(MFLAGS) all

clean:
	@rm -f lib/*.a bin/cuIBM src/*.o
	cd src/solvers; $(MAKE) $(MFLAGS) clean
	cd src/cusp; $(MAKE) $(MFLAGS) clean
	cd external; $(MAKE) $(MFLAGS) clean

force_look:
	true

.PHONY: cuibm
