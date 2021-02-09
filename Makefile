# Specify a particular compiler with "make COMPILER=pgi", etc.
# Specify debugging flags with "make MODE=debug"


# Compiler specific flags

# Use intel compatible flags by default
FFLAGS = -O3
#FFLAGS = -O3 -heap-arrays 64 -ipo -xHost # Optimised (B)
#FFLAGS = -O3 -heap-arrays 64 -ipo -xAVX  # Optimised (W)
#FFLAGS = -O0 -fpe0 -nothreads -traceback -fltconsistency \
#         -C -g -heap-arrays 64 -warn -save-temps -fpic -Wl,-no_pie # Debug
MODULEFLAG = -module $(OBJDIR)

ifeq ($(strip $(MODE)),debug)
  # Autodetect OSX
  SYSTEM := $(shell uname -s)
  FFLAGS = -O0 -fpe0 -nothreads -traceback -fltconsistency \
           -C -g -heap-arrays 64 -warn -fpic
  ifeq ($(strip $(SYSTEM)),Darwin)
    FFLAGS += -Wl,-no_pie
  endif
endif

MPIF90 ?= mpif90
D = -D

# PGI
# ===
ifeq ($(strip $(COMPILER)),pgi)
  FFLAGS = -r8 -fast -fastsse -O3 -Mipa=fast,inline -Minfo # Optimised
  ifeq ($(strip $(MODE)),debug)
    FFLAGS = -Mbounds -g                                     # Debug
  endif
  MODULEFLAG = -module $(OBJDIR)
endif

# Intel
# =====
ifeq ($(strip $(COMPILER)),intel)
  FFLAGS = -O3
  #FFLAGS = -O3 -opt-report 3
  #FFLAGS = -O3 -heap-arrays 64 -ipo -xHost # Optimised (B)
  #FFLAGS = -O3 -heap-arrays 64 -ipo -xAVX  # Optimised (W)
  ifeq ($(strip $(MODE)),debug)
    FFLAGS = -O0 -fpe0 -nothreads -traceback -fltconsistency \
             -C -g -heap-arrays 64 -warn -fpic
    ifeq ($(strip $(SYSTEM)),Darwin)
      FFLAGS += -Wl,-no_pie
    endif
  endif
  MODULEFLAG = -module $(OBJDIR)
  MPIF90 = ftn
endif

# gfortran
# ========
ifeq ($(strip $(COMPILER)),gfortran)
  FFLAGS = -O3
  ifeq ($(strip $(MODE)),debug)
    FFLAGS = -O0 -g -Wall -Wextra -pedantic -fbounds-check \
             -ffpe-trap=invalid,zero,overflow -Wno-unused-parameter \
             -ffpe-trap=underflow,denormal

    GNUVER := $(shell gfortran -dumpversion | head -1 \
        | sed 's/[^0-9\.]*\([0-9\.]\+\).*/\1/')
    GNUMAJOR := $(shell echo $(GNUVER) | cut -f1 -d\.)
    GNUMINOR := $(shell echo $(GNUVER) | cut -f2 -d\.)

    # gfortran-4.3
    GNUGE43 := $(shell expr $(GNUMAJOR) \>= 4 \& $(GNUMINOR) \>= 3)
    ifeq "$(GNUGE43)" "1"
      FFLAGS += -fbacktrace -fdump-core

      # gfortran-4.6
      GNUGE46 := $(shell expr $(GNUMINOR) \>= 6)
      ifeq "$(GNUGE46)" "1"
        FFLAGS += -Wno-unused-dummy-argument

        # gfortran-4.8
        GNUGE48 := $(shell expr $(GNUMINOR) \>= 8)
        ifeq "$(GNUGE48)" "1"
          FFLAGS += -Wno-target-lifetime
        endif
      endif
    endif
  endif
  MODULEFLAG = -I$(OBJDIR) -J$(OBJDIR)
  INFO_FLAGS = -Wno-conversion -fno-range-check
  MPIF90 = mpif90
endif

# g95
# ========
ifeq ($(strip $(COMPILER)),g95)
  FFLAGS = -O3
  ifeq ($(strip $(MODE)),debug)
    FFLAGS = -O0 -g                                        # Debug
  endif
  MODULEFLAG = -fmod=$(OBJDIR)
  MPIF90 = ftn
endif

# IBM Bluegene
# ============
ifeq ($(strip $(COMPILER)),ibm)
  FFLAGS = -O5 -qhot -qipa # Optimised
  ifeq ($(strip $(MODE)),debug)
    FFLAGS = -O0 -C -g -qfullpath -qinfo #-qkeepparm -qflttrap \
          -qnosmp -qxflag=dvz -Q! -qnounwind -qnounroll # Debug
    #FFLAGS = -O0 -qarch=qp -qtune=qp
    #FFLAGS = -qthreaded -qsmp=noauto -qsmp=omp # Hybrid stuff
  endif
  MODULEFLAG = -I$(OBJDIR) -qmoddir=$(OBJDIR)
  MPIF90 = mpixlf90_r

  # IBM compiler needs a -WF to recognise preprocessor directives
  D = -WF,-D
endif

# ARCHER
# ========
ifeq ($(strip $(COMPILER)),archer)
  FFLAGS = -O3
  ifeq ($(strip $(MODE)),debug)
    FFLAGS = -O0 -g -ea -ec -eC -eD -eI -en -hfp_trap -Ktrap=fp -m0 -M1438,7413
  endif
  MODULEFLAG = -em -I/usr/include -I$(OBJDIR) -J$(OBJDIR)
  MPIF90 = ftn
endif

FFLAGS += $(MODULEFLAG)

# Set some of the build parameters
TARGET = epoch3d

# Set pre-processor defines
DEFINES := $(DEF)

# The following are a list of pre-processor defines which can be added to
# the above line modifying the code behaviour at compile time.

# delta-f method
DEFINES += $(D)DELTAF_METHOD
#DEFINES += $(D)DELTAF_DEBUG

# minEPOCH only works for BSPLINE3 but adding this include allow 
# reuse of main EPOCH code without modification
DEFINES += $(D)PARTICLE_SHAPE_BSPLINE3

# Keep global particle counts up to date.
#DEFINES += $(D)PARTICLE_COUNT_UPDATE

# Turn on debugging.
#DEFINES += $(D)PARSER_DEBUG $(D)MPI_DEBUG $(D)SIMPLIFY_DEBUG

# Don't generate any output at all. Useful for benchmarking.
#DEFINES += $(D)NO_IO

# Perform checks on evaluated deck expressions.
# This slows down the code but may be required if floating point exceptions
# are enabled.
#DEFINES += $(D)PARSER_CHECKING

# Enable linking with Trilinos
#DEFINES += $(D)TRILINOS

ifneq (,$(findstring TRILINOS,$(DEFINES)))
  # Up to the user to create, copy, or link this file
  include Makefile.export.Trilinos
  # Trilinos libraries and header files
  INCLUDE = $(Trilinos_INCLUDE_DIRS) $(Trilinos_TPL_INCLUDE_DIRS)
  LIBDIR = $(Trilinos_LIBRARY_DIRS) $(Trilinos_TPL_LIBRARY_DIRS)
  LIB = $(Trilinos_LIBRARIES) $(Trilinos_TPL_LIBRARIES)
  LDFLAGS = $(LIBDIR) $(LIB) $(Trilinos_EXTRA_LD_FLAGS)
  CXX = $(Trilinos_CXX_COMPILER)
  # Check for clang. If so use CXX as linker
  CLANG := $(shell $(CXX) --version | grep 'clang')
  ifneq ($(CLANG),)
    LD = $(CXX)
  else
    LDFLAGS+= -lstdc++
  endif
endif

# Use FC for LD unless specified otherwise
ifeq ($(origin LD),default)
  LD = $(FC)
  LDFLAGS += $(FFLAGS)
endif

# --------------------------------------------------
# Shouldn't need to touch below here
# --------------------------------------------------

all: main

SRCDIR = src
OBJDIR = obj
BINDIR = bin
INCDIR = $(SRCDIR)/include
FC = $(MPIF90)
DATE := $(shell date +%s)
MACHINE := $(shell uname -n)
PREPROFLAGS = $(DEFINES) $(D)_COMMIT='"$(COMMIT)"' $(D)_DATE=$(DATE) \
  $(D)_MACHINE='"$(MACHINE)"' -I$(INCDIR)

FC_INFO := $(shell ${FC} --version 2>/dev/null \
    || ${FC} -V 2>&1 | grep '[a-zA-Z]' | head -n 1)

SRCFILES = balance.F90 boundary.f90 calc_df.F90 current_deposition.F90 \
  custom_laser.f90 deck.f90 diagnostics.F90 epoch3d.F90 fields.f90 finish.f90 \
  helper.F90 ic_module.f90 laser.f90 mpi_routines.F90 mpi_subtype_control.f90 \
  particle_init.F90 particle_temperature.F90 particles.F90 partlist.F90 \
  problem_setup.f90 random_generator.f90 redblack_module.f90 setup.F90 \
  shared_data.F90 strings.f90 timer.f90 utilities.F90 version_data.F90 \
  welcome.F90 deltaf_loader.F90

OBJFILES := $(SRCFILES:.f90=.o)
OBJFILES := $(OBJFILES:.F90=.o)

# C++/Trilinos Interface Files
ifneq (,$(findstring TRILINOS,$(DEFINES)))
  CPPDIR = $(SRCDIR)/cpp

  SRCFILESCPP = JFNKInterface.o JFNKSolver.cpp TrilinosInterface.cpp

  OBJFILES += $(SRCFILESCPP:.cpp=.o)

  # Dependencies
  TRILINOSDEPS = JFNKInterface.o JFNKSolver.o TrilinosInterface.o
endif

INCLUDES = $(INCDIR)/bspline3/b_part.inc $(INCDIR)/bspline3/e_part.inc \
  $(INCDIR)/bspline3/gx.inc $(INCDIR)/bspline3/gxfac.inc \
  $(INCDIR)/bspline3/hx_dcell.inc

OBJFILES := $(OBJFILES)

FULLTARGET = $(BINDIR)/$(TARGET)

VPATH = $(SRCDIR):$(SRCDIR)/deck:$(SRCDIR)/housekeeping:$(SRCDIR)/io:\
  $(SRCDIR)/parser:$(SRCDIR)/user_interaction:$(CPPDIR):$(OBJDIR)

$(SRCDIR)/COMMIT: FORCE
	@sh $(SRCDIR)/gen_commit_string.sh || $(MAKE) $(MAKECMDGOALS)
-include $(SRCDIR)/COMMIT

# Rule to build the fortran files

%.o: %.f90
	$(FC) -c $(FFLAGS) -o $(OBJDIR)/$@ $<

%.o: %.F90
	$(FC) -c $(FFLAGS) -o $(OBJDIR)/$@ $(PREPROFLAGS) $<

%.o: %.cpp
	$(CXX) -c $(CFLAGS) $(Trilinos_INCLUDE_DIRS) $(Trilinos_TPL_INCLUDE_DIRS) -o $(OBJDIR)/$@ $<

main: $(FULLTARGET)
$(FULLTARGET): $(OBJFILES)
	@mkdir -p $(BINDIR)
	$(LD) -o $@ $(addprefix $(OBJDIR)/,$(OBJFILES)) $(LDFLAGS)

clean:
	@rm -rf $(BINDIR) $(OBJDIR)

cleanall: tidy

tidy:
	@rm -rf $(OBJDIR) *~ *.pbs.* *.sh.* $(SRCDIR)/*~ *.log

datatidy:
	@rm -rf Data/*

tarball:
	@cd ..; sh make_tarball.sh

docs:
	@cd Docs; latexmk -pdf report_implemented.tex

cleandocs:
	@cd Docs; rm -f *.aux *.bbl *.blg *.fdb_latexmk *.fls *.log *.out *.pdf

$(OBJFILES): | $(OBJDIR)

$(OBJDIR):
	@mkdir -p $(OBJDIR)

.PHONY: clean cleanall tidy datatidy visit visitclean main FORCE docs cleandocs

# All the dependencies

balance.o: balance.F90 boundary.o redblack_module.o timer.o
boundary.o: boundary.f90 laser.o mpi_subtype_control.o particle_temperature.o \
  partlist.o
calc_df.o: calc_df.F90 shared_data.o
current_deposition.o: current_deposition.F90 shared_data.o utilities.o
custom_laser.o: custom_laser.f90 shared_data.o
deck.o: deck.f90 shared_data.o timer.o fields.o
diagnostics.o: diagnostics.F90 calc_df.o shared_data.o strings.o timer.o
epoch3d.o: epoch3d.F90 balance.o deck.o diagnostics.o fields.o finish.o \
  helper.o ic_module.o mpi_routines.o particles.o problem_setup.o setup.o \
  shared_data.o $(TRILINOSDEPS) welcome.o pat_mpi_lib_interface.o
fields.o: fields.f90 boundary.o
finish.o: finish.f90 laser.o partlist.o
helper.o: helper.F90 boundary.o deltaf_loader.o particle_init.o partlist.o \
  strings.o utilities.o
ic_module.o: ic_module.f90 helper.o setup.o shared_data.o fields.o
laser.o: laser.f90 custom_laser.o
mpi_routines.o: mpi_routines.F90 helper.o
mpi_subtype_control.o: mpi_subtype_control.f90 shared_data.o
particle_init.o: particle_init.F90 shared_data.o random_generator.o
particle_temperature.o: particle_temperature.F90 particle_init.o \
  shared_data.o
particles.o: particles.F90 boundary.o current_deposition.o \
  deltaf_loader.o partlist.o utilities.o
deltaf_loader.o: deltaf_loader.F90 shared_data.o
partlist.o: partlist.F90 shared_data.o
problem_setup.o: problem_setup.f90 fields.o ic_module.o shared_data.o \
  strings.o timer.o
random_generator.o: random_generator.f90
redblack_module.o: redblack_module.f90 partlist.o
setup.o: setup.F90 fields.o laser.o strings.o timer.o version_data.o
shared_data.o: shared_data.F90
strings.o: strings.f90 shared_data.o
timer.o: timer.f90 shared_data.o
utilities.o: utilities.F90 shared_data.o
version_data.o: version_data.F90 $(SRCDIR)/COMMIT
welcome.o: welcome.F90 shared_data.o version_data.o
