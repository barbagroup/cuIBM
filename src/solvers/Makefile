include ../../Makefile.inc

.PHONY: clean NavierStokes

all: NavierStokes

NavierStokes: force_look 
	cd NavierStokes; $(MAKE) $(MFLAGS)
	mv libNavierStokesSolvers.a ../../lib

tests: force_look
	cd NavierStokes; $(MAKE) $(MFLAGS) tests
	mv libTests.a ../../lib

clean:
	cd NavierStokes; $(MAKE) $(MFLAGS) clean

force_look:
	true
