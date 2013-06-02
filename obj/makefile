#
# include teh definition of object and include files
#

include ../makeinclude

#
# the bin directory
#

BIN=$(HOME)/bin

#
# define linker
#

LD = f90

#
# -Wl flags for linker
#
#     -M MAP,e,a     - make a load map sort by address and cross ref
#                      see man mppldrLDFLAGS = -Wl'-M MAP,e,a -X1'
#
#     -Xn            - number of processors to compile for
#                      if not set binary will be relinked at run time.
#                      X1 is special as it will not complain if run
#                      on a different number (for SPMD message passing),
#                      but the values of N$PES which is used in this 
#                      version is the value set at compile time, so it
#                      will not work here: each PE will see N$PES = 1.
#
#     -f zeros       - initialises all variables with 0
#        nan               "                "         Not a Number
#
#     -lmpi          - link to the MPI library (was behind -lapp)
#     -lapp          - link to apprentice libraries

LDFLAGS = -Wl'-X$(NP)' 

#
# Libraries
#

LIBS    = -lmpi -lapp

#
#  names of the executables
#

HYBEXE=$(HOME)/bin/higgs_$(KIND)_$(NP)

#
# how to create the executable
#

$(HYBEXE): $(OBJ) makefile
	@echo ''                                
	@echo '********************************'
	@echo '**linking higgs_$(KIND)_$(NP)'**  
	@echo '********************************'
	@echo ''                                
	@$(LD) $(LDFLAGS) -o $@ $(OBJ) $(LIBS)
	@echo ''                                
	@echo '********************************'
	@echo '************ READY *************'
	@echo '********************************'



