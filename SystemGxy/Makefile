###############################################################################
### Specific makefile for Boltzmann with Hermite-discretizations
###############################################################################

### include generic definitions
### choose which pardiso version to use

#include ../MakefileMKL.def
include ../MakefilePARDISO.def

### $PROBLEM must correspond to the system directory
PROBLEM = SystemGxy

SRCS    = $(SRC_ALL) \
          Equations.cpp \
          SystemGxy.cpp 

include ../Makefile.all
