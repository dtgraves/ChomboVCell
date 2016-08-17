# -*- Mode: Makefile -*- 

## Define the variables needed by Make.example

USE_EB=TRUE
USE_MF=TRUE
# trace the chain of included makefiles
makefiles += ChomboVCell_Code_Mitochondria

# the base name(s) of the application(s) in this directory
ebase = mitochondria

# the location of Chombo lib dir
CHOMBO_HOME = ../../chombo_lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = Workshop EBAMRElliptic EBAMRTools EBTools AMRElliptic AMRTools BoxTools

# relative paths to source code directories
base_dir = .
src_dirs = ../src  ../MFTools

# input file for 'run' target
## NOTE: this relies on the 'foreach' loop used in the 'run-only' target in "Make.rules"
INPUT = test.inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

# application-specific targets

