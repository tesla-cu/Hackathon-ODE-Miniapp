FC := ifort
SRCDIR := src
BUILDDIR := build

# Find all C/C++ source files in the source directory
SOURCES := $(wildcard $(SRCDIR)/*.f90)

# Generate corresponding object file paths in the build directory
OBJECTS := $(patsubst $(SRCDIR)/%.f90, $(BUILDDIR)/%.o, $(SOURCES))

# Define the executable name
EXECUTABLE := test/miniapp.exe

# ----------------------------------------------------------------------------------------
ifeq ($(FC),gfortran)

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

else ifeq ($(FC),ifort)

FFLAGS := -real-size 64 -extend-source 132 -fpic -stand f18 -module ./build # -xCORE-AVX512 AVX512 setting highly recommended for Frontera and Stampede3
LDFLAGS := -limf -lm

DBG1 := -g1 -debug all -traceback -O2  # O2 placed after g1 and debug, otherwise they override O2 to O0
DBG2 := -O0 -g2 -debug all -traceback -prec-div -fpe-all=0 -warn nounused -check all,noarg_temp_created -fp-model=strict  # -g2 == -g
DBG3 := -O0 -g3 -debug all -traceback -prec-div -fpe-all=0 -warn nounused -check all -fp-model=strict,source

SYNTAX := -syntax-only -warn all,errors -diag-error-limit=5

# Intel fortran option -O1 optimizes for small executable size, not faster performance, unlike gfortran
OPT1 := -O2 -fp-model=precise -fprotect-parens
OPT2 := -O2 -fprotect-parens
OPT3 := -O3 -ipo # -qopt-zmm-usage=high may help or hurt if added here
OPT4 := -O3 -no-prec-div -ipo # -qopt-zmm-usage=high may help or hurt if added here
# Other stuff: it's possible settings like -heap-arrays, and
# -mcmodel=medium could help performance and/or avoid runtime memory errors.

endif

# the first recipe in this list is the default when running `make` without specifying a recipe.
.PHONY: debug
debug: OPTIONS=$(DBG2)
debug: $(EXECUTABLE)

.PHONY: fast
fast: OPTIONS=$(OPT2)
fast: $(EXECUTABLE)

.PHONY: syntax
syntax: OPTIONS=$(SYNTAX)
syntax: $(EXECUTABLE)

$(EXECUTABLE): $(OBJECTS)
	$(FC) $(FFLAGS) $(OPTIONS) $(OBJECTS) -o $@ $(LDFLAGS)

# Rule to compile source files into object files
$(BUILDDIR)/%.o: $(SRCDIR)/%.f90
	@mkdir -p $(BUILDDIR)
	$(FC) $(FFLAGS) -c $< -o $@

.PHONY: format
format:
	fprettify --indent 4 -w 4 -l 300 --whitespace-intrinsics true\
	 --enable-decl --enable-replacements --c-relations *.f90

# This removes everything in the .build/ directory
.PHONY: clean
clean:
	rm -rf $(BUILDDIR) $(EXECUTABLE) $(EXECUTABLE).dSYM

# this removes everything in the .build/ directory AND any accidental outputs made in the project directory
.PHONY: realclean
realclean: clean
	rm -rf *.o *.mod *.dSYM
