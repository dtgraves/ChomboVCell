# -*- Mode: Makefile -*- 

## Define the variables needed by Make.example

## Note: requires USE_EB=TRUE

# trace the chain of included makefiles
makefiles += example_EBAMRPoisson_tests_amr_solutionError

# the base name(s) of the application(s) in this directory
ebase = backwardEulerHeat

# the location of Chombo lib dir
CHOMBO_HOME = ../../Chombo/lib
# CHOMBO_HOME = $(HOME)/_svnchombo/lib

# names of Chombo libraries needed by this program, in order of search.
LibNames = EBAMRElliptic EBAMRTools Workshop EBTools AMRElliptic AMRTools BoxTools

# relative paths to source code directories
base_dir = .
src_dirs = .

# input file for 'run' target
INPUT = regression.inputs

# shared code for building example programs
include $(CHOMBO_HOME)/mk/Make.example

# application-specific variables

# application-specific targets
