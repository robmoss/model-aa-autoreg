#
# Define the project layout, source code modules, and target binary.
#

SRC_DIR = ./src
OBJ_DIR = ./obj
BIN_DIR = ./bin

SRC = $(wildcard $(SRC_DIR)/*.f90) $(wildcard $(SRC_DIR)/*.f)
NAMES = $(basename $(notdir $(SRC)))
OBJ = $(NAMES:%=$(OBJ_DIR)/%.o)
MOD = $(NAMES:%=$(OBJ_DIR)/%.mod)

TARGET = $(BIN_DIR)/run_model

#
# Compilation settings, including optimisations and debugging checks.
#

# Use the GCC Fortran compiler.
FC = gfortran

# Default to compiling with optimisations and without debug information.
# Debugging can be enabled from the command-line:
#
# 1. Remove all object files:
#        make distclean
#
# 2. Recompile the model with debugging enabled:
#        DEBUG=true make
#
DEBUG ?= false

# Check whether the user has overridden the DEBUG setting.
ifeq ($(strip $(DEBUG)), true)
    OPT = -g -fbacktrace
else
    OPT = -O3
endif

# Check the gfortran version.
GCC_48 := $(shell expr `gcc -dumpversion | sed -e 's/\.\([0-9][0-9]\)/\1/g' -e 's/\.\([0-9]\)/0\1/g' -e 's/^[0-9]\{3,4\}$$/&00/'` \>= 40800)

# Search the compilation directory for object and module files (*.o & *.mod).
FFLAGS := -J $(OBJ_DIR)
# Enable common warnings, including extensions to Fortran 95, and disable
# implicit typing (except where overridden by explicit IMPLICIT statements).
FFLAGS += -Wall $(OPT) -pedantic -Wno-unused-variable -fimplicit-none
# Warn when uninitialised variables are used.
FFLAGS += -Wuninitialized
ifeq "$(GCC_48)" "1"
    # Initialise REAL variables to a signalling NaN, so that the use of an
    # uninitialised REAL variable will trigger a floating point exception.
    FFLAGS += -finit-real=snan
    # Enable all run-time checks (e.g., bounds checking for array subscripts).
    FFLAGS += -fcheck=all
endif
# Abort the program (SIGFPE) when a floating point exception occurs:
#   - An invalid operation (e.g., SQRT(-1.0));
#   - Division by zero;
#   - Floating point overflow.
FFLAGS += -ffpe-trap=invalid,zero,overflow

#
# The default target is to compile the binary.
#

.PHONY: default
default: $(TARGET)

#
# Program compilation.
#

$(TARGET): $(OBJ)
	@$(FC) -fopenmp -o $@ $(OBJ)

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f90 | $(OBJ_DIR)
	@$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR)/%.o: $(SRC_DIR)/%.f | $(OBJ_DIR)
	@$(FC) $(FFLAGS) -c $< -o $@

$(OBJ_DIR):
	@mkdir $(OBJ_DIR)


#
# Module compilation order.
#

# The object files are intermediate files; if the target binary is
# newer than the source files, the object files are not needed.
.SECONDARY: $(OBJ)

# Dependencies that control the compilation order of the modules.
$(OBJ_DIR)/main.o: $(OBJ_DIR)/tests.o $(OBJ_DIR)/sngfr.o $(OBJ_DIR)/cmdline.o
$(OBJ_DIR)/main.o: $(OBJ_DIR)/cmdline.o
$(OBJ_DIR)/tests.o: $(OBJ_DIR)/model_base.o $(OBJ_DIR)/model_glom.o

#
# Removal of temporary and generated files.
#

# Remove all temporary files.
.PHONY: clean
clean:
	@rm -f $(OBJ) $(MOD)

# Remove all generated files, including the target binary.
.PHONY: distclean
distclean: clean
	@rm -f $(TARGET) $(wildcard *.ssv) $(wildcard *.pdf)
