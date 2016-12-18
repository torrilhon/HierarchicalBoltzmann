###############################################################################
### Specific makefile for NSF system 
###############################################################################

### include generic definitions
### choose which pardiso version to use

#include ../MakefileMKL.def
include ../MakefilePARDISO.def

### $PROBLEM must correspond to the system directory
PROBLEM = SystemNSF

SRCS    = $(SRC_ALL) \
          EquationsNSF.cpp \
          SystemNSF.cpp 

include ../Makefile.all
