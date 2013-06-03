FFLAGS =  -ffixed-line-length-256 
LFLAGS = 

OBJDIR = .
SRCDIR = ../f
EXEDIR = ../runs

#
# directory where include files may be found
#

INC     = $(SRCDIR)/INCLUDE

# works on plm using mpich
#
#   F90 = /usr/lib64/mpich2/bin/mpif90
#   INCPATH = -I $(INC)

# works on plm using openmpi
#   F90 = mpif90
#   INCPATH = -I $(INC)

F90 = openmpif90
INCPATH = -I $(INC)

#endif

SRC1 = \
$(SRCDIR)/adresses.f \
$(SRCDIR)/a_force.f \
$(SRCDIR)/agauge.f \
$(SRCDIR)/ahiggs2.f \
$(SRCDIR)/ahiggs3.f \
$(SRCDIR)/confget.f \
$(SRCDIR)/confsave.f \
$(SRCDIR)/correlator.f \
$(SRCDIR)/finalize.f \
$(SRCDIR)/gauss.f \
$(SRCDIR)/heat_kpi.f \
$(SRCDIR)/initialize.f \
$(SRCDIR)/it_gauge_fixing.f \
$(SRCDIR)/landau_gauge.f \
$(SRCDIR)/luxury.f \
$(SRCDIR)/matmul.f \
$(SRCDIR)/meas1.f \
$(SRCDIR)/meas2a.f \
$(SRCDIR)/measinit.f \
$(SRCDIR)/or_gauge_fixing1.f \
$(SRCDIR)/paraminit.f \
$(SRCDIR)/paratypeconstruct.f \
$(SRCDIR)/polyakov.f \
$(SRCDIR)/readdata.f \
$(SRCDIR)/renorm.f \
$(SRCDIR)/scalarcorr.f \
$(SRCDIR)/staples.f \
$(SRCDIR)/testmax.f \
$(SRCDIR)/topology.f \
$(SRCDIR)/typeconstruct.f \
$(SRCDIR)/u_f_renorm.f \



#
#Names of the INCLUDE FILES
#

INC_FILES=$(INC)/const.inc                \
          $(INC)/filenames.inc            \
          $(INC)/gaugefield.inc           \
          $(INC)/input_parameter.inc      \
          $(INC)/measurement.inc          \
          $(INC)/parallel_mpi_types.inc   \
          $(INC)/parallel_parameter.inc   \
          $(INC)/parameter.inc            \
          $(INC)/pointer.inc              \
          $(INC)/times.inc                \


MAIN1 =  \
$(SRCDIR)/main.f

SRC1O     = $(patsubst $(SRCDIR)%.f, $(OBJDIR)%.o, $(SRC1))


all: $(EXEDIR)/main.exe

$(EXEDIR)/main.exe: $(OBJDIR)/main.o $(SRC1O) 
	$(F90) $(LFLAGS) $(FFLAGS) $(OBJDIR)/main.o $(SRC1O) -o $@ $(LIBS) \


$(OBJDIR)/%.o : $(SRCDIR)/%.f $(INC_FILES)
	$(F90) $(INCPATH) -c $(FFLAGS) $<  -o $@

