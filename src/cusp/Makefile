include ../../Makefile.inc

.PHONY: wrappers clean

all: wrappers

wrappers: force_look
	cd wrapped; $(MAKE) $(MFLAGS)
	@mv libCuspWrapped.a ../../lib

clean:
	cd wrapped; $(MAKE) $(MFLAGS) clean

force_look:
	true
