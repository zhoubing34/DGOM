# $Id: macros.make.in 2264 2015-12-22 15:42:59Z wkliao $

# The purpose of this file is to contain common make(1) macros.
# It should be processed by every execution of that utility.



# POSIX shell.  Shouldn't be necessary -- but is under IRIX 5.3.
SHELL		= /bin/sh
RM		= rm
LN_S		= ln -s

# Installation Directories:
# SRCDIR	= @SRCDIR@
prefix		= /Users/mac/Documents/Model/DGOM/thirdLibrary/install
INCDIR		= $(prefix)/include
LIBDIR		= $(prefix)/lib
BINDIR		= $(prefix)/bin
MANDIR		= $(prefix)/man
BUILDDIR	= /Users/mac/Documents/Model/DGOM/thirdLibrary/parallel-netcdf-1.7.0
LIBRARY		= /Users/mac/Documents/Model/DGOM/thirdLibrary/parallel-netcdf-1.7.0/src/lib/libpnetcdf.a


# Useful tools
M4		= m4
M4FLAGS		= 
EGREP		= /usr/bin/grep -E

# AC_PROG_SED and AC_PROG_GREP are only available on autoconf 2.60 and later
# SED		= @SED@
# GREP		= /usr/bin/grep
SED		= sed
GREP		= grep

# Preprocessing:
DEFS		= -DHAVE_CONFIG_H
FC_DEFINE	= 
CPP		= /usr/local/Cellar/mpich2/1.4.1/bin/mpicc -E
CPPFLAGS	= $(INCLUDES) $(DEFS) 
CXXCPPFLAGS     = $(INCLUDES) $(DEFS) 
FPP		= 
FPPFLAGS	= $(INCLUDES)  


# Compilation:
MPICC		= /usr/local/Cellar/mpich2/1.4.1/bin/mpicc
MPICXX		= 
MPIF77		= 
MPIF90		= 

SEQ_CC          = /usr/bin/gcc

# debugging and optimization options for compiling and linking
CFLAGS		= -O3
CXXFLAGS	= 
F77FLAGS	=  
F90FLAGS	=  

# compiler options for different file extensions: .f .F .f90 .F90
F77FLAGS_f	= 
F77FLAGS_F	= 
F90FLAGS_f90	= 
F90FLAGS_F90	= 

# preprocessor options for different file extensions: .f .F .f90 .F90
F77PPFLAGS_f	= 
F77PPFLAGS_F	= 
F90PPFLAGS_f90	= 
F90PPFLAGS_F90	= 

# NETCDF.MOD	= @NETCDF_MOD@
CC_MAKEDEPEND	= false

COMPILE.c	= $(MPICC)  $(CFLAGS)       $(CPPFLAGS) -c
COMPILE.cxx	= $(MPICXX) $(CXXFLAGS)     $(CXXCPPFLAGS) -c
COMPILE.f	= $(MPIF77) $(F77FLAGS_f)   $(FPPFLAGS) $(F77FLAGS) -c
COMPILE.f90	= $(MPIF90) $(F90FLAGS_f90) $(FPPFLAGS) $(F90FLAGS) -c
COMPILE.F	= $(MPIF77) $(F77FLAGS_F)   $(FPPFLAGS) $(F77FLAGS) $(F77PPFLAGS_F) -c
COMPILE.F90	= $(MPIF90) $(F90FLAGS_F90) $(FPPFLAGS) $(F90FLAGS) $(F90PPFLAGS_F90) -c
# In PnetCDF, we follow the file extension convention that .F and .F90 files
# require preprocessing, while .f and .f90 do not.


# Linking:
FLIBS		= 
FCLIBS		= 
F90LIBS		= 
FLDFLAGS	= 
F90LDFLAGS	= 
LDFLAGS		=  
LIBS		= 

LINK.c		= $(MPICC)  $(CFLAGS)   -o $@
LINK.cxx	= $(MPICXX) $(CXXFLAGS) -o $@
LINK.F77	= $(MPIF77) $(F77FLAGS) -o $@
LINK.F90	= $(MPIF90) $(F90FLAGS) -o $@

TEST_MPIRUN	= mpiexec -n NP
TEST_OUTDIR	= .
TEST_SEQRUN	= 

# Manual pages:
WHATIS		= whatis
# The following macro should be empty on systems that don't
# allow users to create their own manual-page indexes.
MAKEWHATIS_CMD	= 


# Misc. Utilities:
AR		= ar
ARFLAGS		= cru
AWK		= @AWK@
RANLIB		= ranlib
INSTALL 	= /usr/bin/install -c
INSTALL_DATA	= ${INSTALL} -m 644
TARFLAGS	= -chf


# Dummy macros: used only as placeholders to silence GNU make.  They are
# redefined, as necessary, in subdirectory makefiles.
HEADER		= dummy_header
HEADER1		= dummy_header1
HEADER2		= dummy_header2
HEADER3		= dummy_header3
MANUAL		= dummy_manual
PROGRAM		= dummy_program


# Distribution macros:
FTPDIR		= /home/ftp/pub/$(PACKAGE)
FTPBINDIR	= @FTPBINDIR@

PNETCDF_VERSION_MAJOR = 1
PNETCDF_VERSION_MINOR = 7
PNETCDF_VERSION_SUB   = 0
PNETCDF_VERSION_PRE   = 
PNETCDF_VERSION       = 1.7.0

