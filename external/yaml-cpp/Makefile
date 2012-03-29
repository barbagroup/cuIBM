.PHONY: lib install clean

lib: force_look
	cd src; $(MAKE) $(MFLAGS)

install: force_look
	@mv libyaml-cpp.a ../lib

clean: force_look
	cd src; $(MAKE) $(MFLAGS) clean
	rm ../lib/libyaml-cpp.a

force_look:
	true
