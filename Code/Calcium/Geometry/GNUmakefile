# -*- Mode: Makefile -*- 

## Define the variables needed by Make.example

## Note: requires USE_EB=TRUE

# trace the chain of included makefiles
makefiles += ChomboVCell_Code_Calcium_Geometry

# the base name(s) of the application(s) in this directory
ebase = vcellOne

# the location of Chombo lib dir
CHOMBO_HOME = ../../../../Chombo/lib
#CHOMBO_HOME = $(HOME)/_cvschombo/lib
# names of Chombo libraries needed by this program, in order of search.
LibNames = MFTools Workshop EBAMRTools EBTools AMRTools BoxTools BaseTools

# relative paths to source code directories
base_dir = .
src_dirs = ../../src

# input file for 'run' target
## NOTE: this relies on the 'foreach' loop used in the 'run-only' target in "Make.rules"
INPUT = test.inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

# application-specific targets

