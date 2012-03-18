include Makefile.inc

.PHONY: cuibm src/solvers/libsolvers.a clean

main: cuibm

cuibm: lib/libNavierStokesSolvers.a src/cuIBM.o src/io.o
	nvcc $? -o bin/cuIBM

#src/solvers/libsolvers.a:  force_look
lib/libNavierStokesSolvers.a:
	cd src/solvers; $(MAKE) $(MFLAGS)

clean:
	@rm -f lib/*.a bin/cuIBM src/*.o
	cd src/solvers; $(MAKE) $(MFLAGS) clean

force_look:
	true

.PHONY: cuibm
