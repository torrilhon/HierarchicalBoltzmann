###############################################################################
### Specific makefile for R13 system 
###############################################################################

### include generic definitions
### choose which pardiso version to use

#include ../MakefileMKL.def
include ../MakefilePARDISO.def

### $PROBLEM must correspond to the system directory
PROBLEM = SystemR13

SRCS    = $(SRC_ALL) \
          EquationsR13.cpp \
          SystemR13.cpp 

include ../Makefile.all
