## Process this file with automake to produce Makefile.in
#
bin_PROGRAMS    =
noinst_PROGRAMS =
if ENABLE_CDI_LIB
  bin_PROGRAMS    += cdi
else
  noinst_PROGRAMS += cdi
endif
#
noinst_PROGRAMS    += createtable
#
cdi_SOURCES         = cdi.c printinfo.h
cdi_LDADD           = $(top_builddir)/src/libcdi.la @LIBS@
if ENABLE_ALL_STATIC
cdi_LDFLAGS         = -all-static
endif
#
createtable_SOURCES = createtable.c
createtable_LDADD   = $(top_builddir)/src/libcdi.la
if ENABLE_ALL_STATIC
createtable_LDFLAGS = -all-static
endif
#
AM_CPPFLAGS         = -I$(top_srcdir)/src
#
CLEANFILES          = `ls *~`
