#
# Define relative path...
#
SRCDIR  = sources/main/
SRCDIR2 = sources/variabledef/
SRCDIR3 = sources/wave/
SRCDIR4 = sources/fourier/
SRCDIR5 = sources/timeintegration/
SRCDIR6 = sources/IO/
OBJDIR  = obj/
LIBDIR  = /usr/local/lib/
BINDIR  = bin/
LINKLIB = $(LIBDIR)libfftw3.a $(LIBDIR)liblapack.a $(LIBDIR)librefblas.a
#
# ifort compiler
#
#FC=ifort
#
#OPTFLAGS= -O3 -xHOST -ipo -ip -module $(OBJDIR) # -fp-model precise this option to keep exactly same fp description
#DBFLAGS= -O0 -traceback -check all -warn all -module $(OBJDIR)
#
# gfortran compiler (optimized for Mac Pro with intel corei7: avx not working on Mac)
FC=gfortran
#
OPTFLAGS= -O3 -march=corei7 -msse2 -funroll-loops -fno-protect-parens -ffast-math
DBFLAGS = -O0 -Wline-truncation -fbounds-check -fpack-derived -fbacktrace -ffpe-summary=all -fimplicit-none -fcheck=all -Wall -Wtabs -Wextra -Wunderflow  -Wno-zerotrip -finit-integer=inf -finit-real=nan -ffpe-trap=zero,overflow
FLAGMOD1= -J $(OBJDIR) #Flag for writing modules in $(OBJ)
FLAGMOD2= -I $(OBJDIR) #Flag for reading modules in $(OBJ)
#
nomcode=HOS-NWT
#
#define all modules and subroutines
#
module=\
	wavemaking\
	filtering\
	dealiasing\
	resol_HOS\
	resol_wmkr\
	initial_condition

module2=\
	variables\
	common_vars\
	type
	
module3=\
	definition2\
	config_def\
	read_ocean_txt\
	wave_def\
	linear_wave

module4=\
	fourier_FFTW

module5=\
	runge_kutta\
	variablechange

module6=\
	input\
	vol_energy\
	velocities\
	output	

subroutine=\

# define the list of source files.
SRCS1 = $(addprefix $(SRCDIR), $(addsuffix .f90, $(module))) $(addprefix $(SRCDIR), $(addsuffix .f90, $(subroutine))) 
SRCS2 = $(addprefix $(SRCDIR2), $(addsuffix .f90, $(module2))) 
SRCS3 = $(addprefix $(SRCDIR3), $(addsuffix .f90, $(module3))) 
SRCS4 = $(addprefix $(SRCDIR4), $(addsuffix .f90, $(module4)))
SRCS5 = $(addprefix $(SRCDIR5), $(addsuffix .f90, $(module5)))
SRCS6 = $(addprefix $(SRCDIR6), $(addsuffix .f90, $(module6)))
#
SRCS0 = $(SRCS1) $(SRCS2) $(SRCS3) $(SRCS4) $(SRCS5) $(SRCS6)
SRCS  = $(addprefix $(SRCDIR), $(addsuffix .f90, $(nomcode))) $(SRCS0)

build:Release

# top-level rule, to compile everything.
all:clean depend Release

# Targets for compilation

Release: FFLAGS1 = $(OPTFLAGS) $(FLAGMOD1)
Release: FFLAGS2 = $(OPTFLAGS) $(FLAGMOD2)
Release: $(BINDIR)$(nomcode) $(SRCS)

Debug: FFLAGS1 = $(DBFLAGS) $(FLAGMOD1)
Debug: FFLAGS2 = $(DBFLAGS) $(FLAGMOD2)
Debug: $(BINDIR)$(nomcode) $(SRCS)

# now add a line to include the dependency list.
include $(OBJDIR).depend

Objlink = $(addprefix $(OBJDIR), $(addsuffix .o, $(module))) $(addprefix $(OBJDIR), $(addsuffix .o, $(subroutine))) $(addprefix $(OBJDIR), $(addsuffix .o, $(module2))) $(addprefix $(OBJDIR), $(addsuffix .o, $(module3))) $(addprefix $(OBJDIR), $(addsuffix .o, $(module4))) $(addprefix $(OBJDIR), $(addsuffix .o, $(module5))) $(addprefix $(OBJDIR), $(addsuffix .o, $(module6)))

# rule to link the program
$(BINDIR)$(nomcode): $(SRCDIR)$(nomcode:%=%.f90) $(Objlink) $(OBJDIR).depend
	$(FC) $(FFLAGS2) -o $(BINDIR)$(nomcode) $(SRCDIR)$(nomcode:%=%.f90) $(Objlink) $(LINKLIB)

# now comes a meta-rule for compiling any ".f90" source file.
$(OBJDIR)%.o: $(SRCDIR)%.f90 $(OBJDIR).depend
	$(FC) $(FFLAGS1) -c $< -o $@

$(OBJDIR)%.o: $(SRCDIR2)%.f90 $(OBJDIR).depend
	$(FC) $(FFLAGS1) -c $< -o $@

$(OBJDIR)%.o: $(SRCDIR3)%.f90 $(OBJDIR).depend
	$(FC) $(FFLAGS1) -c $< -o $@

$(OBJDIR)%.o: $(SRCDIR4)%.f90 $(OBJDIR).depend
	$(FC) $(FFLAGS1) -c $< -o $@

$(OBJDIR)%.o: $(SRCDIR5)%.f90 $(OBJDIR).depend
	$(FC) $(FFLAGS1) -c $< -o $@

$(OBJDIR)%.o: $(SRCDIR6)%.f90 $(OBJDIR).depend
	$(FC) $(FFLAGS1) -c $< -o $@

clean:
	rm -f $(OBJDIR)*
	rm -f $(OBJDIR).depend
	rm -f $(BINDIR)$(nomcode)
	@echo "all cleaned up!"

# rule for building dependency lists, and writing them to a file
# named ".depend".

depend $(OBJDIR).depend:
	@if test -d bin; then echo "exists"; else mkdir bin; fi
	@if test -d obj; then echo "exists"; else mkdir Obj; fi
	@if test -d bin/Results; then echo "exists"; else mkdir bin/Results; fi
	rm -f $(OBJDIR).depend
	makedepf90 -b $(OBJDIR) $(SRCS0) > $(OBJDIR).depend

.PHONY: clean depend
