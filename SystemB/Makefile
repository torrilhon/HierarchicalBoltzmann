###############################################################################
### Specific makefile for system (B)
###############################################################################

### include generic definitions
### choose which pardiso version to use

#include ../MakefileMKL.def
include ../MakefilePARDISO.def

### $PROBLEM must correspond to the system directory
PROBLEM = SystemB

SRCS    = $(SRC_ALL) \
          EquationsB.cpp \
          SystemB.cpp 

include ../Makefile.all
