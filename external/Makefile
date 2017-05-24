YAMLCPP_DIR="yaml-cpp-0.5.1"

.PHONY: libyaml-cpp

all: libyaml-cpp

libyaml-cpp:
	cd $(YAMLCPP_DIR); $(MAKE) $(MFLAGS); $(MAKE) $(MFLAGS) install

clean:
	cd $(YAMLCPP_DIR); $(MAKE) $(MFLAGS) clean

force_look:
	true
