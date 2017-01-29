#
# Copyright (C) 2012, Northwestern University and Argonne National Laboratory
# See COPYRIGHT notice in top-level directory.
#
# $Id: Makefile.in 1468 2013-10-26 16:53:18Z wkliao $
#
# doc/Makefile.  Generated from Makefile.in by configure.

srcdir = .


include ../macros.make

LATEX  = /usr/texbin/latex
DVIPDF = /usr/local/bin/dvipdf
ifeq (yes, yes)
    API_GUIDE = pnetcdf-api.pdf
else
    API_GUIDE =
endif

all: $(API_GUIDE)

pnetcdf-api.pdf: $(srcdir)/pnetcdf-api.tex
	$(LATEX) $<
	$(LATEX) $<
	$(DVIPDF) pnetcdf-api

PACKING_LIST = Makefile.in \
               c_api.tex \
               data_mode_api.tex \
               flexible_api.tex \
               latex8.sty \
               Makefile.in \
               multiple_io.tex \
               nonblocking.tex \
               pnetcdf-api.bbl \
               pnetcdf-api.tex \
               porting_notes.txt \
               symbol_renaming.txt

GARBAGE = *.dvi *.log *.aux *.ps $(API_GUIDE)

include $(srcdir)/../rules.make

