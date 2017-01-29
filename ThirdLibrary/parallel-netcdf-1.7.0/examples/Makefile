#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id: Makefile.in 2318 2016-02-04 00:18:26Z wkliao $
#
# examples/Makefile.  Generated from Makefile.in by configure.

srcdir  = .


include ../macros.make

GARBAGE =

SUBDIRS = C tutorial

ifeq (no, yes)
SUBDIRS += CXX
endif

ifeq (no, yes)
SUBDIRS += F77
ifeq (, yes)
SUBDIRS += F90
endif
endif

PACKING_LIST = Makefile.in README

PACKING_SUBDIRS = C tutorial F77 F90 CXX

all: $(SUBDIRS)
$(SUBDIRS):
	$(MAKE) $(MFLAGS) -C $@

ptest: all
ifeq (no, yes)
	echo "Parallel test is disabled because coverage analysis was enabled"
else
	(echo "=============================================")  && \
	(echo "    Parallel testing on 4 MPI processes")  && \
	(echo "=============================================")  && \
	( for d in $(SUBDIRS) ; do \
		$(MAKE) $(MFLAGS) -C $$d ptest4 ; \
	done ) ;
endif

ptests: all
ifeq (no, yes)
	echo "Parallel test is disabled because coverage analysis was enabled"
else
	for i in 3 4 8 ; do \
	(echo "=============================================")  && \
	(echo "    Parallel testing on $$i MPI processes")  && \
	(echo "=============================================")  && \
	( for d in $(SUBDIRS) ; do \
	$(MAKE) $(MFLAGS) -C $$d ptest$$i ; \
	done ) ; \
	done ;
endif

INSTALLDIRS = $(SUBDIRS:%=install-%)
install: all $(INSTALLDIRS)
$(INSTALLDIRS): 
	$(MAKE) $(MFLAGS) -C $(@:install-%=%) install

UNINSTALLDIRS = $(SUBDIRS:%=uninstall-%)
uninstall: $(UNINSTALLDIRS)
$(UNINSTALLDIRS): 
	$(MAKE) $(MFLAGS) -C $(@:uninstall-%=%) uninstall

include $(srcdir)/../rules.make

