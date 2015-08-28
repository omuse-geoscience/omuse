## Process this file with automake to produce Makefile.in
#
EXTRA_DIST = cdi_read_f2003.f90 cdi_write_f2003.f90
#
noinst_PROGRAMS = cdi_write cdi_write_ens cdi_write_hybrid cdi_read cdi_copy
if CREATE_ISOC
  noinst_PROGRAMS += cdi_read_f2003 cdi_write_f2003
endif
#
AM_CPPFLAGS              = -I$(top_srcdir)/src
AM_FCFLAGS               = $(FPP_INCOPT)$(top_srcdir)/src $(FC_MOD_FLAG)../src
#
cdi_write_SOURCES        = cdi_write.c
cdi_write_LDADD          = $(top_builddir)/src/libcdi.la
#
cdi_write_ens_SOURCES    = cdi_write_ens.c
cdi_write_ens_LDADD      = $(top_builddir)/src/libcdi.la
#
cdi_write_hybrid_SOURCES = cdi_write_hybrid.c
cdi_write_hybrid_LDADD   = $(top_builddir)/src/libcdi.la
#
cdi_read_SOURCES         = cdi_read.c
cdi_read_LDADD           = $(top_builddir)/src/libcdi.la
#
cdi_copy_SOURCES         = cdi_copy.c
cdi_copy_LDADD           = $(top_builddir)/src/libcdi.la
#
cdi_read_f2003_SOURCES   = cdi_read_f2003.f90
cdi_read_f2003_LDADD     = $(top_builddir)/src/libcdi.la $(top_builddir)/src/mo_cdi.o
cdi_read_f2003.$(OBJEXT): $(top_builddir)/src/mo_cdi.$(FCMODEXT)
#
cdi_write_f2003_SOURCES  = cdi_write_f2003.f90
cdi_write_f2003_LDADD    = $(top_builddir)/src/libcdi.la $(top_builddir)/src/mo_cdi.$(OBJEXT)
cdi_write_f2003.$(OBJEXT): $(top_builddir)/src/mo_cdi.$(FCMODEXT)

#
if ENABLE_ALL_STATIC
cdi_write_LDFLAGS        = -all-static
cdi_read_LDFLAGS         = -all-static
cdi_copy_LDFLAGS         = -all-static
cdi_read_f2003_LDFLAGS   = -all-static
cdi_write_f2003_LDFLAGS  = -all-static
endif

