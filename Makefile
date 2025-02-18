F90 := gfortran
FFLAGS := -fdefault-real-8 -fdefault-double-8 -fimplicit-none -fPIC -pipe -std=f2018 -J./build
LDFLAGS := -lm

# Og allows some O1 optimizations while making debugging cleaner/clearer than O0
# also, plain `-g` is equivalent to `-g2`
DBG1 := -Og -g1 -ffpe-trap=invalid,zero,overflow
DBG2 := -Og -g2 -ffpe-trap=invalid,zero,overflow,underflow\
        -Wall -Wno-unused-variable -Wno-maybe-uninitialized -Wno-unused-dummy-argument\
		-fcheck=all,no-array-temps -fmax-errors=5
DBG3 := -O0 -g2 -ffpe-trap=invalid,zero,overflow,underflow -Wall -Wextra -Werror -fcheck=all -fmax-errors=5

SYNTAX := -fsyntax-only -Wall -Wextra -Werror -fmax-errors=5

optW := -Warray-temporaries -Wconversion-extra -Wimplicit-interface -Wimplicit-procedure\
	    -Wrealloc-lhs-all -Wfrontend-loop-interchange
OPT1 := -O1 $(optW)
OPT2 := -O2 $(optW)
OPT3 := -O3
OPT4 := -O3 -ffast-math -fno-protect-parens

# ----------------------------------------------------------------------------------------
SRC_DIR := ./src
vpath %.f90 $(SRC_DIR)
target := ./test/miniapp.exe
MOD_SRCS := src/pprk4.f90 src/miniapp_rkc.f90 src/ncarles_rkc.f90 src/integrators.f90 src/chemistry.f90
MODS := $(MOD_SRCS:.f90=.mod)

# the first recipe in this list is the default when running `make` without specifying a recipe.
debug: $(MODS)
	$(F90) $(FFLAGS) $(DBG1) $(MOD_SRCS) src/chem_ode_miniapp.f90 -o $(target) $(LDFLAGS)

fast: $(MODS)
	$(F90) $(FFLAGS) $(OPT3) $(MOD_SRCS) src/chem_ode_miniapp.f90 -o $(target) $(LDFLAGS)

syntax: $(MODS)
	$(F90) $(FFLAGS) $(SYNTAX) $(MOD_SRCS) src/chem_ode_miniapp.f90

# make step for module files
%.mod : %.f90
	$(MPI) $(FFLAGS) -c $<

.PHONY: format
format:
	fprettify --indent 4 -w 4 -l 300 --whitespace-intrinsics true\
	 --enable-decl --enable-replacements --c-relations *.f90

# This removes everything in the .build/ directory
.PHONY: clean
clean:
	rm -rf ./build/*
	rm -rf $(target) $(target).dSYM

# this removes everything in the .build/ directory AND any accidental outputs made in the project directory
.PHONY: realclean
realclean: clean
	rm -rf *.o *.mod *.dSYM
