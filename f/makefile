#
# include the definition of object and include files
#

include ../makeinclude

#
# das parent directory
#

PARENT=..

#
# define fortran compiler
#

FC = f90

#
# flags for the compiler
#
# -c     as usual
#
# -g     debugging (cripples perfomance)
#
# -dp    disable DOUBLE PRECISION map to REAL
#
# -O 0,1,2,3 optimisation level
#    aggress 
#    unroll2
#    split2
#
# -e A   use apprentice (perfomance monitoring) needs -lapp
#        some impact on performance
#

offstr = 
onstr  =

#FFLAGS =  -c -d$(offstr) -e$(onstr) -I INCLUDE -N 132

FFLAGS =  -c  -I INCLUDE  -N 132 -R abcnsCE

all: $(OBJ) 

$(PARENT)/obj/$(KIND)/%.o:%.f $(INC_FILES) makefile $(PARENT)/makeinclude
	@echo ''                                
	@echo '********************************'
	@echo '********** compiling ***********'
	@echo '********************************'
	@echo ''                                
	$(FC) $(FFLAGS) $(OPTFLAGS) $<
	@mv $*.o `dirname $@`

