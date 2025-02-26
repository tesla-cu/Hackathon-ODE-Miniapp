COMPILER := nvidia
FC := nvfortran # OR mpif90, etc.
SRCDIR := src
BUILDDIR := build

SOURCES := src/miniapp_rkc.f90 src/chemistry.f90 src/chem_ode_miniapp.f90

# Generate corresponding object file paths in the build directory
OBJECTS := $(patsubst $(SRCDIR)/%.f90, $(BUILDDIR)/%.o, $(SOURCES))

# Define the executable name
EXECUTABLE := test/miniapp.exe

# ----------------------------------------------------------------------------------------
ifeq ($(COMPILER),nvidia)
FFLAGS := -r8 -module ./build -acc=multicore -Minfo=ftn,all
OPT2 := -g -O2

else ifeq ($(COMPILER),cray)

FFLAGS := -s default64 -f PIC -ef -J ./build
LDFLAGS := # I don't think anything is necessary here

# USE THESE FOR ACTUAL DEBUGGING
SYNTAX := -g # NO SYNTAX-ONLY OPTION!?
DBG1 := -O2 -G1 -m3  # O2 placed after g1 and debug, otherwise they override O2 to O0
DBG2 := -O1 -G1 -Ktrap=fp -eD -h fp1,scalar1,vector1,add_paren
DBG3 := -O0 -G0 -Ktrap=fp -eD -m1 -h add_paren # -O0 implies fp0, scalar0, vector0, etc.

# USE THESE FOR CODE PROFILING
# -- Very unoptimized floating-point operations.
#    Every way in which you can force Cray to do math slower is turned on
OPT1 := -O0 -G0 -h add_paren # -O0 implies fp0, scalar0, vector0, etc.
# -- "Normal" floating-point operations
OPT2 := -O2 -G2
# -- Very optimized FLOPs, any way in which you can trade accuracy for speed is turned on.
OPT3 := -O2 -G2 -h scalar3,vector3,fp4 # fma on at fp1 or higher

# USE THESE FOR MAXIMUM COMPILER OPTIMIZATION
# -- Once O3 and ipo are turned on, I don't think you can use -g anymore.
OPT4 := -O3 -eo
OPT5 := -O3 -h aggress,scalar3,vector3,cache3,fp4

else ifeq ($(COMPILER),intel)

FFLAGS := -real-size 64 -extend-source 132 -fpic -stand f18 -module ./build\
 -I"${MKLROOT}/include"
# -march=core-avx2 is enabled on Derecho for AMD Epyc Milan CPUs
LDFLAGS := -L${MKLROOT}/lib/intel64 -lmkl_rt -lpthread -lm -ldl

# USE THESE FOR ACTUAL DEBUGGING
SYNTAX := -syntax-only -warn all,errors -diag-error-limit=5
DBG1 := -g2 -debug all -traceback -O2  # O2 placed after g1 and debug, otherwise they override O2 to O0
DBG2 := -g2 -debug all -traceback -fpe-all=0 -warn nounused -check all,noarg_temp_created -fp-model=strict  # -g2 == -g
DBG3 := -g3 -debug all -traceback -prec-div -fprotect-parens -fpe-all=0 -warn nounused -check all -fp-model=strict,source

# USE THESE FOR CODE PROFILING
# -- Very unoptimized floating-point operations.
#    Every way in which you can force intel to do math slower is turned on
OPT1 := -g -O2 -fp-model=strict,source -fprotect-parens -prec-div
# -- "Normal" floating-point operations
OPT2 := -g -O2
# -- Very optimized FLOPs, any way in which you can trade accuracy for speed is turned on.
OPT3 := -g -O2 -fp-model=fast -fast-transcendentals -fma -no-prec-div -nostandard-realloc-lhs

# USE THESE FOR MAXIMUM COMPILER OPTIMIZATION
# -- Once O3 and ipo are turned on, I don't think you can use -g anymore.
OPT4 := -O3 -ipo # -qopt-zmm-usage=high may help or hurt if added here
OPT5 := -O3 -ipo -fast-transcendentals -no-prec-div -nostandard-realloc-lhs # -fimf-precision=high or link to MKL!
# Other stuff: it's possible settings like -qopt-zmm-usage=high, and
# -mcmodel=medium could help performance and/or avoid runtime memory errors.

else ifeq ($(COMPILER),gnu)

FFLAGS := -fdefault-real-8 -fdefault-double-8 -fimplicit-none -fPIC -pipe -std=f2018 -J./build
LDFLAGS := -lm

SYNTAX := -fsyntax-only -Wall -Wextra -Werror -fmax-errors=5
# Og allows some O1 optimizations while making debugging cleaner/clearer than O0
# also, plain `-g` is equivalent to `-g2`
DBG1 := -Og -g1 -ffpe-trap=invalid,zero,overflow
DBG2 := -Og -g2 -ffpe-trap=invalid,zero,overflow,underflow\
		-Wall -Wno-unused-variable -Wno-maybe-uninitialized -Wno-unused-dummy-argument\
		-fcheck=all,no-array-temps -fmax-errors=5
DBG3 := -O0 -g2 -ffpe-trap=invalid,zero,overflow,underflow -Wall -Wextra -Werror -fcheck=all -fmax-errors=5

# CODE PROFILING
warn := -Warray-temporaries -Wconversion-extra -Wimplicit-interface -Wimplicit-procedure\
		-Wrealloc-lhs-all -Wfrontend-loop-interchange
WRN1 := -g2 -O1 $(dbgW)
WRN2 := -g2 -O2 $(dbgW)

OPT1 := -g2 -O2
OPT2 := -g2 -O2 -ffast-math -fno-protect-parens

# HARDCORE OPTIMIZATION
OPT3 := -O3
OPT4 := -O3 -ffast-math -fno-protect-parens

endif

# the first recipe in this list is the default when running `make` without specifying a recipe.
.PHONY: debug
debug: OPTIONS=$(DBG1)
debug: $(EXECUTABLE)

.PHONY: profile
profile: OPTIONS=$(OPT2)
profile: $(EXECUTABLE)

.PHONY: fast
fast: OPTIONS=$(OPT4)
fast: $(EXECUTABLE)

.PHONY: syntax
syntax: OPTIONS=$(SYNTAX)
syntax: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(FC) $(FFLAGS) $(OPTIONS) $(OBJECTS) -o $@ $(LDFLAGS)

# Rule to compile source files into object files
$(BUILDDIR)/%.o: $(SRCDIR)/%.f90
	@mkdir -p $(BUILDDIR)
	$(FC) $(FFLAGS) $(OPTIONS) -c $< -o $@

.PHONY: format
format:
	fprettify --indent 4 -w 4 -l 300 --whitespace-intrinsics true\
	 --enable-decl --enable-replacements --c-relations src/*.f90

# This removes everything in the .build/ directory
.PHONY: clean
clean:
	rm -rf $(BUILDDIR) $(EXECUTABLE) $(EXECUTABLE).dSYM

# this removes everything in the .build/ directory AND any accidental outputs made in the project directory
.PHONY: realclean
realclean: clean
	rm -rf *.o *.mod *.dSYM
