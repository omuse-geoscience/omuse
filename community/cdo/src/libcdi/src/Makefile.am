## Process this file with automake to produce Makefile.in
#
EXTRA_DIST = cdilib.c mo_cdi.f90

noinst_LTLIBRARIES = libcdiresunpack.la

if ENABLE_CDI_LIB
  lib_LTLIBRARIES = libcdi.la
  include_HEADERS = cdi.h cdi.inc
else
  noinst_LTLIBRARIES += libcdi.la
endif

AM_CFLAGS = $(PPM_CORE_CFLAGS) $(YAXT_CFLAGS) $(MPI_C_INCLUDE)

libcdi_la_SOURCES = 	 \
	basetime.c     	 \
	basetime.h	 \
	binary.c	 \
	binary.h	 \
	calendar.c 	 \
	calendar.h	 \
	cdf.c            \
	cdf.h	 	 \
	cdf_int.c	 \
	cdf_int.h	 \
	cdi.h	 	 \
	cdi_error.c      \
	cdi_limits.h	 \
	cdi_util.c       \
	cgribex.h	 \
	cgribexlib.c  	 \
	datetime.h	 \
	dmemory.c      	 \
	dmemory.h	 \
	cksum.c		\
	cksum.h		\
	cdi_cksum.c	\
	cdi_cksum.h	\
	create_uuid.h	 \
	dtypes.h	 \
	error.c        	 \
	error.h	 	 \
	extra.h	 	 \
	extralib.c       \
	file.c         	 \
	file.h	 	 \
	gaussgrid.c      \
	gaussgrid.h	 \
	gribapi.c  	 \
	gribapi.h	 \
	gribapi_utilities.c \
	gribapi_utilities.h \
	grid.c           \
	grid.h	 	 \
	ieg.h	 	 \
	ieglib.c         \
	input_file.c	\
	input_file.h	\
	institution.c  	 \
	institution.h  	 \
	model.c        	 \
	model.h        	 \
	namespace.c      \
	namespace.h      \
	serialize.h	\
	serialize.c	\
	referenceCounting.c \
	referenceCounting.h \
	resource_handle.c\
	resource_handle.h\
	service.h	 \
	servicelib.c     \
	stream_cdf.c     \
	stream_cdf.h	 \
	stream_cgribex.c \
	stream_cgribex.h \
	stream_ext.c     \
	stream_ext.h	 \
	stream_grb.c     \
	stream_grb.h     \
	stream_gribapi.c \
	stream_gribapi.h \
	stream_history.c \
	stream_ieg.c     \
	stream_ieg.h	 \
	stream_fcommon.c \
	stream_fcommon.h \
	cdi_int.c        \
	cdi_int.h	 \
	stream_record.c  \
	stream_srv.c     \
	stream_srv.h	 \
	stream_var.c     \
        subtype.c        \
        subtype.h        \
	swap.h	 	 \
	table.c        	 \
	table.h	 	 \
	tablepar.h	 \
	taxis.c          \
	taxis.h	         \
	timebase.c 	 \
	timebase.h	 \
	tsteps.c         \
	util.c         	 \
	varscan.c      	 \
	varscan.h        \
	version.c      	 \
	vlist.c 	 \
	vlist.h	         \
	vlist_att.c 	 \
	vlist_att.h 	 \
	vlist_var.c 	 \
	vlist_var.h	 \
	zaxis.c		 \
	zaxis.h		 \
	stream.c         \
	swap.c          \
	iterator.c          \
	iterator.h          \
	iterator_fallback.c \
	iterator_fallback.h \
	iterator_grib.c     \
	iterator_grib.h

# cfortran.h is an optional part of libcdi
libcdi_la_USE_FC_extra_sources = \
	cfortran.h       \
	cdiFortran.c

if USE_FC
libcdi_la_SOURCES += 	 \
	$(libcdi_la_USE_FC_extra_sources)
endif
EXTRA_DIST += \
	$(libcdi_la_USE_FC_extra_sources)



libcdiresunpack_la_SOURCES = \
	resource_unpack.c

if USE_MPI
if ENABLE_CDI_LIB
lib_LTLIBRARIES += libcdipio.la
include_HEADERS += cdipio.h cdipio.inc
else
noinst_LTLIBRARIES += libcdipio.la
endif
endif

libcdipio_la_SOURCES =							\
	cdipio.h							\
	cdipioFortran.c							\
	cfortran.h							\
	pio.c								\
	pio.h								\
	pio_comm.c							\
	pio_comm.h							\
	pio_dbuffer.c							\
	pio_impl.h							\
	pio_interface.c							\
	pio_interface.h							\
	pio_mpinonb.c							\
	pio_record_send.c						\
	pio_posixasynch.c						\
	pio_posixfpguardsendrecv.c					\
	pio_posixnonb.c							\
	pio_list_set.c							\
	resource_unpack.h						\
	resource_unpack.c						\
	pio_rpc.h							\
	pio_client.c							\
	pio_client.h							\
	pio_server.c							\
	pio_server.h							\
	pio_serialize.h							\
	pio_serialize.c							\
	pio_util.c							\
	pio_util.h							\
	pio_cdf_int.c							\
	pio_cdf_int.h

libcdipio_la_LIBADD    = libcdi.la $(PPM_CORE_LIBS) $(YAXT_LIBS) $(MPI_C_LIB)

libcdi_la_LIBADD    = $(UUID_C_LIB)
#
#cdilib.c:
#	$(top_srcdir)/src/make_cdilib $(top_srcdir)/src
#
#cdilib.o: cdilib.c
#	$(COMPILE) -c $<

LOCALTARGETS  =
if ENABLE_CDI_LIB
LOCALTARGETS += pkgconfig/cdi.pc
if USE_MPI
LOCALTARGETS += pkgconfig/cdipio.pc
endif
endif

if CREATE_ISOC
LOCALTARGETS += mo_cdi.$(OBJEXT)
endif
#
mo_cdi.f90: $(top_srcdir)/src/cdi.h $(top_srcdir)/interfaces/f2003/bindGen.rb
	ruby $(top_srcdir)/interfaces/f2003/bindGen.rb $(top_srcdir)/src/cdi.h $@
#
mo_cdi.$(OBJEXT): mo_cdi.f90
	$(FC) $(FCFLAGS) -c $(FCFLAGS_f90) $<
#
pkgconfig/cdi.pc: pkgconfig/cdi.pc.in ../config.status
	(cd .. ; ./config.status src/pkgconfig/cdi.pc)

pkgconfig/cdipio.pc: pkgconfig/cdipio.pc.in ../config.status
	(cd .. ; ./config.status src/pkgconfig/cdipio.pc)
#
all-local: $(LOCALTARGETS)
#
CLEANFILES  = `ls *~`
#CLEANFILES += cdilib.c
if CREATE_ISOC
CLEANFILES += mo_cdi.$(FCMODEXT) mo_cdi.$(OBJEXT)
endif

if ENABLE_CDI_LIB
PKGCONFIG_FILES = pkgconfig/cdi.pc
CLEANFILES += pkgconfig/cdi.pc
if USE_MPI
CLEANFILES += pkgconfig/cdipio.pc
PKGCONFIG_FILES += pkgconfig/cdipio.pc
endif

install-exec-local: $(PKGCONFIG_FILES)
	$(mkinstalldirs) "$(DESTDIR)$(libdir)/pkgconfig"
	@for pkgconfigfile in $(PKGCONFIG_FILES) ; do \
		echo $(install_sh_DATA) $$pkgconfigfile \
		"$(DESTDIR)$(libdir)/$$pkgconfigfile" ; \
		$(install_sh_DATA) $$pkgconfigfile \
		"$(DESTDIR)$(libdir)/$$pkgconfigfile" ; \
	done

uninstall-local:
	@for pkgconfigfile in $(PKGCONFIG_FILES) ; do \
		echo rm -f "$(DESTDIR)$(libdir)/$$pkgconfigfile" ; \
		rm -f "$(DESTDIR)$(libdir)/$$pkgconfigfile" ; \
	done
	-@rmdir "$(DESTDIR)$(libdir)/pkgconfig"

endif

install-exec-hook:
	-@rmdir "$(DESTDIR)$(libdir)"
install-data-hook:
	-@rmdir "$(DESTDIR)$(includedir)"

