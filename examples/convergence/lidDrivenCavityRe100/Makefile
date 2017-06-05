CURRENT_DIR=$(shell pwd)

.PHONY: all navierstokes directforcing fadlunetal feamodified dfmodified dfimproved diffusion

all: navierstokes directforcing fadlunetal feamodified dfmodified dfimproved diffusion

navierstokes:
	cd $(CURRENT_DIR)/NavierStokes; make all

directforcing:
	cd $(CURRENT_DIR)/DirectForcing; make all

fadlunetal:
	cd $(CURRENT_DIR)/FadlunEtAl; make all

feamodified:
	cd $(CURRENT_DIR)/FEAModified; make all

dfmodified:
	cd $(CURRENT_DIR)/DFModified; make all

dfimproved:
	cd $(CURRENT_DIR)/DFImproved; make all

diffusion:
	cd $(CURRENT_DIR)/Diffusion; make all
