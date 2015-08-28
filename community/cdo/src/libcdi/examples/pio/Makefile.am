noinst_PROGRAMS=collectData collectDataNStreams

AM_CPPFLAGS = -I$(top_srcdir)/src $(YAXT_CFLAGS)
AM_FCFLAGS  = $(FPP_INCOPT)$(top_srcdir)/src $(MPI_FC_INCLUDE)
AM_CFLAGS = $(MPI_C_INCLUDE)

if USE_MPI
noinst_PROGRAMS+=compareResourcesArray
AM_FCFLAGS += $(FPP_DEFOPT)USE_MPI
endif

if USE_FC
 noinst_PROGRAMS+=collectData2003
endif

collectData_SOURCES=collectData.c

if USE_MPI
LDADD=$(top_builddir)/src/libcdipio.la $(MPI_C_LIB)
collectData2003_LDADD=$(top_builddir)/src/libcdipio.la $(MPI_FC_LIB)
else
LDADD=$(top_builddir)/src/libcdi.la
collectData2003_LDADD=$(top_builddir)/src/libcdi.la
endif

collectData2003_SOURCES=collectData2003.F90

collectDataNStreams_SOURCES = collectDataNStreams.c

compareResourcesArray_SOURCES = compareResourcesArray.c
