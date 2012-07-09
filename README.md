To compile, set the environment variable CUSP_DIR to point to the directory 
with the Cusp library. In a bash shell, add the following line to  
~/.bashrc or ~/.bash_profile:

	export CUSP_DIR=/path/to/cusp_folder

And reload the file:

	source ~/.bashrc
	
Compile all files:
	
	make

Run the test:
	
	bin/cuIBM

Known issues: (July 9th, 2012)

ATTENTION: The documentation and the code base will be made more complete in a couple of week's time. Some issues still need to be sorted out. Everything that is available now will solve flows for arbitrary bodies, but options for moving bodies, time stepping and preconditioners do not work currently.

* Linker step fails on Mac OS X
* Only the first order explicit Euler method for convection and the first order implicit Euler method for diffusion have been tested
* Forces are not calculated for the direct forcing method (FadlunEtAlSolver)