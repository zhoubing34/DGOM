#
# Copyright (C) 2014, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id: Makefile.in 2264 2015-12-22 15:42:59Z wkliao $
#
# src/utils/pnetcdf_version/Makefile.  Generated from Makefile.in by configure.

srcdir = .


include ../../../macros.make

DEFS  = -DMPICC="\"$(MPICC)\""   -DCFLAGS="\"$(CFLAGS)\""
DEFS += -DCONFIGURE_ARGS_CLEAN="\"--with-mpi=/usr/local/Cellar/mpich2/1.4.1/ --prefix=/Users/mac/Documents/Model/DGOM/externlib/install --disable-cxx --disable-fortran MPICC=/usr/local/Cellar/mpich2/1.4.1/bin/mpicc CC=/usr/local/Cellar/mpich2/1.4.1/bin/mpicc CFLAGS=-O3\""
DEFS += -DPNETCDF_VERSION="\"1.7.0\""
DEFS += -DPNETCDF_RELEASE_DATE="\"03 Mar 2016\""

ifeq (no, yes)
DEFS += -DMPICXX="\"$(MPICXX)\"" -DCXXFLAGS="\"$(CXXFLAGS)\""
endif

ifeq (no, yes)
DEFS += -DMPIF77="\"$(MPIF77)\"" -DF77FLAGS="\"$(F77FLAGS)\""
DEFS += -DMPIF90="\"$(MPIF90)\"" -DF90FLAGS="\"$(F90FLAGS)\""
endif

C_SOURCES	= pnetcdf_version.c

PROGRAM		= pnetcdf_version
MANUAL		= pnetcdf_version.1

PACKING_LIST    = $(C_SOURCES) $(MANUAL) \
                  Makefile.in depend

GARBAGE		= $(PROGRAM)

all: $(PROGRAM)

pnetcdf_version: pnetcdf_version.c
	$(SEQ_CC) $(DEFS) -o $@ $<

install: $(PROGRAM) $(MANUAL)
	$(INSTALL) -d -m 755 $(MANDIR)/man1
	$(INSTALL_DATA) $(srcdir)/$(MANUAL) $(MANDIR)/man1/$(MANUAL)
	$(INSTALL) -d $(BINDIR)
	$(INSTALL) -m 755 $(PROGRAM) $(BINDIR)/$(PROGRAM)

uninstall:
	$(RM) -f $(BINDIR)/$(PROGRAM)
	$(RM) -f $(MANDIR)/man1/$(MANUAL)

include $(srcdir)/../../../rules.make
include $(srcdir)/depend

.PHONY: $(LIBRARY)
