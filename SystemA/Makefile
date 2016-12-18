###############################################################################
### Specific makefile for system (A)
###############################################################################

### include generic definitions
### choose which pardiso version to use

#include ../MakefileMKL.def
include ../MakefilePARDISO.def

### $PROBLEM must correspond to the system directory
PROBLEM = SystemA

SRCS    = $(SRC_ALL) \
          EquationsA.cpp \
          SystemA.cpp 

include ../Makefile.all
