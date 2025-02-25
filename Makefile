.SUFFIXES:

F90 := mpif90
FFLAGS := -fdefault-real-8 -fdefault-double-8 -fimplicit-none -fPIC -pipe -std=f2018 -J./.build
LDFLAGS := -limf

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
OPT5 := -r8 -132 -O2

target := miniapp
sources := pprk4.f90 miniapp_rkc.f90 ncarles_rkc.f90 integrators.f90 chemistry.f90 chem_ode_miniapp.f90

# the first recipe in this list is the default when running `make` without specifying a recipe.
fast:
	$(F90) $(OPT5) -o $(target) $(sources) $(LDFLAGS)

debug:
	$(F90) $(FFLAGS) $(DBG1) $(sources) -o $(target) $(LDFLAGS)

syntax:
	$(F90) $(FFLAGS) $(SYNTAX) $(sources)

.PHONY: debug1
debug1: debug

debug2:
	$(F90) $(FFLAGS) $(DBG2) $(sources) -o $(target) $(LDFLAGS)

debug3:
	$(F90) $(FFLAGS) $(DBG3) $(sources) -o $(target) $(LDFLAGS)

opt1:
	$(F90) $(FFLAGS) $(OPT1) $(sources) -o $(target) $(LDFLAGS)

opt2:
	$(F90) $(FFLAGS) $(OPT2) $(sources) -o $(target) $(LDFLAGS)

.PHONY: opt3
opt3: fast

opt4:
	$(F90) $(FFLAGS) $(OPT4) $(SOURCES) -o $(target) $(LDFLAGS)
	$(F90) $(FFLAGS) $(OPT4) $(SOURCES) -o $(target) $(LDFLAGS)

.PHONY: format
format:
	fprettify --indent 4 -w 4 -l 300 --whitespace-intrinsics true\
	 --enable-decl --enable-replacements --c-relations *.f90

# This removes everything in the .build/ directory
.PHONY: clean
clean:
	rm -rf .build/*

# this removes everything in the .build/ directory AND any accidental outputs made in the project directory
.PHONY: realclean
realclean: clean
	rm -rf $(exec) *.o *.mod *.dSYM
