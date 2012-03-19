include ../../Makefile.inc

.PHONY: clean NavierStokes

all: NavierStokes

#NavierStokes: force_look
NavierStokes: force_look 
	cd NavierStokes; $(MAKE) $(MFLAGS)
	@mv libNavierStokesSolvers.a ../../lib

clean:
	cd NavierStokes; $(MAKE) $(MFLAGS) clean

force_look:
	true
