####################################################################
#  CLAPACK make include file.                                      #
#  CLAPACK, Version 3.1.1.1                                        #
#  October 2008                                                    #
####################################################################
#
#  The machine (platform) identifier to append to the library names
#
PLAT = _LINUX
PREFIX =/Users/mac/Documents/Model/DGOM/externlib/install
#  
#  Modify the CC and CFLAGS definitions to refer to the
#  compiler and desired compiler options for your machine.  NOOPT
#  refers to the compiler options desired when NO OPTIMIZATION is
#  selected.  Define LOADER and LOADOPTS to refer to the loader and
#  desired load options for your machine.
#
#######################################################
# This is used to compile C libary
CC        = /usr/local/Cellar/mpich2/1.4.1/bin/mpicc
# if no wrapping of the blas library is needed, uncomment next line
#CC        = gcc -DNO_BLAS_WRAP
CFLAGS    = -O3
LOADER    = /usr/local/Cellar/mpich2/1.4.1/bin/mpicc
LOADOPTS  = 
NOOPT     = -O0 
DRVCFLAGS = $(CFLAGS)
F2CCFLAGS = $(CFLAGS)
#######################################################################
#
#  The archiver and the flag(s) to use when building archive (library)
#  If you system has no ranlib, set RANLIB = echo.
#
ARCH     = ar
ARCHFLAGS= cr
RANLIB   = ranlib
#
#  The location of the libraries to which you will link.  (The 
#  machine-specific, optimized BLAS library should be used whenever
#  possible.)
#

# Using reference BLAS
BLASLIB      = ../../libblas.a
# Using Vectib
#BLASLIB = -framework vecLib
# Using Fortran BLAS interface of ATLAS without wrapper (CLAPACK has to be compiled with -DNO_BLAS_WRAP)
#BLASLIB = -lf77blas -latlas
# Using Fortran BLAS interface of ATLAS with wrapper
#BLASLIB = ../../libfblaswr.a -lf77blas -latlas
# Using CBLAS interface of ATLAS with wrapper
#BLASLIB = ../../libcblaswr.a -lcblas -latlas
# Using Goto without wrapper (CLAPACK has to be compiled with -DNO_BLAS_WRAP)
#BLASLIB = -lgoto -lpthread
# Using Goto with wrapper
#BLASLIB = ../../libfblaswr.a -lgoto -lpthread 
LAPACKLIB    = liblapack.a
F2CLIB       = ../../F2CLIBS/libf2c.a
TMGLIB       = libtmglib.a
EIGSRCLIB    = libeigsrc.a
LINSRCLIB    = liblinsrc.a

