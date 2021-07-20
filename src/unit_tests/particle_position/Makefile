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
LDFLAGS = $(FFLAGS)

# Set some of the build parameters
TARGET = particle_load_test

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


# --------------------------------------------------
# Shouldn't need to touch below here
# --------------------------------------------------

all: main

SRCDIR = src
MAINSRCDIR = ../..
OBJDIR = obj
BINDIR = bin
INCDIR = $(MAINSRCDIR)/include
FC = $(MPIF90)
DATE := $(shell date +%s)
MACHINE := $(shell uname -n)
PREPROFLAGS = $(DEFINES) -I$(INCDIR)

FC_INFO := $(shell ${FC} --version 2>/dev/null \
    || ${FC} -V 2>&1 | grep '[a-zA-Z]' | head -n 1)

SRCFILES = particle_load_test.F90 shared_data.F90 random_generator.f90 \
  particle_init.F90

OBJFILES := $(SRCFILES:.f90=.o)
OBJFILES := $(OBJFILES:.F90=.o)

INCLUDES = $(INCDIR)/bspline3/b_part.inc $(INCDIR)/bspline3/e_part.inc \
  $(INCDIR)/bspline3/gx.inc $(INCDIR)/bspline3/gxfac.inc \
  $(INCDIR)/bspline3/hx_dcell.inc

OBJFILES := $(OBJFILES)

FULLTARGET = $(BINDIR)/$(TARGET)

VPATH = $(SRCDIR):$(MAINSRCDIR):$(MAINSRCDIR)/housekeeping:$(MAINSRCDIR)/user_interaction:$(OBJDIR)

# Rule to build the fortran files

%.o: %.f90
	$(FC) -c $(FFLAGS) -o $(OBJDIR)/$@ $<

%.o: %.F90
	$(FC) -c $(FFLAGS) -o $(OBJDIR)/$@ $(PREPROFLAGS) $<

main: $(FULLTARGET)
$(FULLTARGET): $(OBJFILES)
	@mkdir -p $(BINDIR)
	$(FC) -o $@ $(addprefix $(OBJDIR)/,$(OBJFILES)) $(LDFLAGS)

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

particle_load_test.o: particle_load_test.F90 random_generator.o shared_data.o \
  particle_init.o
particle_lading.o: particle_init.F90 random_generator.o shared_data.o
random_generator.o: random_generator.f90
shared_data.o: shared_data.F90
