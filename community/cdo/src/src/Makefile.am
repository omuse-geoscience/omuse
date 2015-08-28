## Process this file with automake to produce Makefile.in
noinst_LTLIBRARIES = libcdo.la
libcdo_la_SOURCES =            \
               cdo_int.h       \
               compare.h       \
               cdo_pthread.c   \
               cdo_vlist.c     \
               cdo_getopt.c    \
               cdo_getopt.h    \
               cdo_history.c   \
               after_sptrans.c \
               after_fctrans.c \
               after_dvtrans.c \
               after_vertint.c \
               after_vertint.h \
               after_namelist.c\
               afterburnerlib.c\
               afterburner.h   \
               vct_l191.h      \
               constants.h     \
               constants.c     \
               color.c         \
               color.h         \
               commandline.c   \
               const.h         \
               counter.h       \
               datetime.c      \
               datetime.h      \
               dmemory.h       \
               ecacore.c       \
               ecacore.h       \
               ecautil.c       \
               ecautil.h       \
               error.h         \
               etopo.h         \
               temp.h          \
               mask.h          \
               exception.c     \
               expr.c          \
               expr.h          \
               expr_lex.c      \
               expr_yacc.c     \
               expr_yacc.h     \
               features.c      \
               field.c         \
               field.h         \
               field2.c        \
               fieldc.c        \
               fieldmem.c      \
               fieldmer.c      \
               fieldzon.c      \
               functs.h        \
               gradsdeslib.c   \
               gradsdeslib.h   \
               grid.c          \
               grid.h          \
               grid_area.c     \
               grid_gme.c      \
               grid_lcc.c      \
               grid_rot.c      \
               gridreference.c \
               griddes.c       \
               griddes.h       \
               griddes_h5.c    \
               griddes_nc.c    \
               hetaeta.c       \
               hetaeta.h       \
               institution.c   \
               interpol.c      \
               interpol.h      \
               job.c           \
               juldate.c       \
               kdtree.c        \
               kdtree.h        \
               kdtree_fast.c   \
               kdtree_fast.h   \
               grid_search.c   \
               grid_search.h   \
               kvlist.c        \
               kvlist.h        \
               list.c          \
               list.h          \
	       merge_sort2.c   \
	       merge_sort2.h   \
               modules.c       \
               modules.h       \
               namelist.c      \
               namelist.h      \
               normal.c        \
               nth_element.c   \
               nth_element.h   \
               operator_help.h \
               par_io.c        \
               par_io.h        \
               percentiles.c   \
               percentiles.h   \
               pipe.c          \
               pipe.h          \
               pragma_omp_atomic_update.h \
               printinfo.h     \
               process.c       \
               process.h       \
               pstream.c       \
               pstream.h       \
               pstream_write.h \
               pstream_int.h   \
               pthread_debug.c \
               pthread_debug.h \
               readline.c      \
               realtime.c      \
               remap.h         \
               remaplib.c      \
               remapsort.c     \
               remap_scrip_io.c \
               remap_search_reg2d.c \
               remap_search_latbins.c \
               remap_store_link.c \
               remap_store_link.h \
               remap_store_link_cnsrv.c \
               remap_store_link_cnsrv.h \
               remap_conserv.c \
               remap_conserv_scrip.c \
               remap_distwgt_scrip.c \
               remap_bicubic_scrip.c \
               remap_bilinear_scrip.c \
               stdnametable.c  \
               stdnametable.h  \
               specspace.c     \
               specspace.h     \
               statistic.c     \
               statistic.h     \
               table.c         \
               text.c          \
               text.h          \
               timebase.h      \
               timer.c         \
               userlog.c       \
               util.c          \
               util.h          \
               zaxis.c
libcdo_la_SOURCES +=           \
               clipping/clipping.c           \
               clipping/clipping.h           \
               clipping/area.c               \
               clipping/area.h               \
               clipping/ensure_array_size.c  \
               clipping/ensure_array_size.h  \
               clipping/geometry.h           \
               clipping/grid.h               \
               clipping/points.h             \
               clipping/dep_list.h           \
               clipping/grid_cell.c          \
               clipping/grid_cell.h          \
               clipping/intersection.c       \
               clipping/utils.c              \
               clipping/utils.h 
#
bin_PROGRAMS = cdo
#
cdo_SOURCES  = cdo.c
cdo_SOURCES += Adisit.c        \
               Afterburner.c   \
               Arith.c         \
               Arithc.c        \
               Arithdays.c     \
               Arithlat.c      \
               CDItest.c       \
               CDIread.c       \
               CDIwrite.c      \
               Cat.c           \
               Change.c        \
               Change_e5slm.c  \
               Cloudlayer.c    \
               Collgrid.c      \
               Command.c       \
               Comp.c          \
               Compc.c         \
               Complextorect.c \
               Cond.c          \
               Cond2.c         \
               Condc.c         \
               Consecstat.c    \
               Copy.c          \
               Deltime.c       \
               Derivepar.c     \
               Detrend.c       \
               Diff.c          \
               Distgrid.c      \
               Duplicate.c     \
               EOFs.c          \
	       Eof3d.c         \
               EcaIndices.c    \
               Echam5ini.c     \
               Enlarge.c       \
               Enlargegrid.c   \
               Ensstat.c       \
               Ensstat3.c      \
               Ensval.c        \
               Eofcoeff.c      \
	       Eofcoeff3d.c    \
               Exprf.c         \
               FC.c            \
               Filedes.c       \
               Fillmiss.c      \
               Filter.c        \
               Fldrms.c        \
               Fldstat.c       \
               Fldstat2.c      \
               Fourier.c       \
               Gengrid.c       \
               Gradsdes.c      \
               Gridboxstat.c   \
               Gridcell.c      \
               Gridsearch.c    \
               Harmonic.c      \
               Hi.c            \
               Histogram.c     \
               Importamsr.c    \
               Importbinary.c  \
               Importcmsaf.c   \
               Importobs.c     \
               Info.c          \
               Input.c         \
               Intgrid.c       \
               Intgridtraj.c   \
               Intlevel.c      \
               Intlevel3d.c    \
               Intntime.c      \
               Inttime.c       \
               Intyear.c       \
               Invert.c        \
               Invertlev.c     \
               Isosurface.c    \
               Kvl.c           \
               Log.c           \
               Maskbox.c       \
               Mastrfu.c       \
               Math.c          \
               Merge.c         \
               Mergegrid.c     \
               Mergetime.c     \
               Merstat.c       \
               Monarith.c      \
               Mrotuv.c        \
               Mrotuvb.c       \
               Ninfo.c         \
               Nmltest.c       \
               Output.c        \
               Outputgmt.c     \
               Pack.c          \
               Pinfo.c         \
               Pressure.c      \
               Regres.c        \
               Remap.c         \
               Remapeta.c      \
               Replace.c       \
               Replacevalues.c \
               Rhopot.c        \
               Rotuv.c         \
               Runpctl.c       \
               Runstat.c       \
               Seascount.c     \
               Seaspctl.c      \
               Seasstat.c      \
               Selbox.c        \
               Select.c        \
               Seloperator.c   \
               Selrec.c        \
               Seltime.c       \
               Selvar.c        \
               Set.c           \
               Setbox.c        \
               Setgatt.c       \
               Setgrid.c       \
               Sethalo.c       \
               Setmiss.c       \
               Setpartab.c     \
               Setrcaname.c    \
               Settime.c       \
               Setzaxis.c      \
               Showinfo.c      \
               Sinfo.c         \
               Smooth9.c       \
               Sort.c          \
               Sorttimestamp.c \
               Specinfo.c      \
               Spectral.c      \
               Spectrum.c      \
               Split.c         \
               Splitrec.c      \
               Splitsel.c      \
               Splittime.c     \
               Splityear.c     \
               SSOpar.c        \
               Subtrend.c      \
               Tee.c           \
               Templates.c     \
               Test.c          \
               Tests.c         \
               Timcount.c      \
               Timpctl.c       \
               Timselpctl.c    \
               Timselstat.c    \
               Timsort.c       \
               Timstat.c       \
               Timstat2.c      \
               Timstat3.c      \
               Tinfo.c         \
               Tocomplex.c     \
               Transpose.c     \
               Trend.c         \
               Trms.c          \
               Tstepcount.c    \
               Vardup.c        \
               Vargen.c        \
               Varrms.c        \
               Vertintml.c     \
               Vertintap.c     \
               Vertstat.c      \
               Vertwind.c      \
               Wct.c           \
               Wind.c          \
               Writegrid.c     \
               Writerandom.c   \
               YAR.c           \
               Yearmonstat.c   \
               Ydayarith.c     \
               Ydaypctl.c      \
               Ydaystat.c      \
               Ydrunpctl.c     \
               Ydrunstat.c     \
               Yhourarith.c    \
               Yhourstat.c     \
               Ymonarith.c     \
               Ymonpctl.c      \
               Ymonstat.c      \
               Yseaspctl.c     \
               Yseasstat.c     \
               Zonstat.c       \
               cdo.h


if ENABLE_MAGICS
cdo_SOURCES += Magplot.c       \
               Magvector.c     \
               Maggraph.c      \
               template_parser.h   \
               template_parser.c   \
               results_template_parser.h   \
               results_template_parser.c   \
               magics_template_parser.h   \
               magics_template_parser.c   \
               StringUtilities.h   \
               StringUtilities.c   \
               CdoMagicsMapper.h   \
               CdoMagicsMapper.c
endif

cdo_CPPFLAGS = -I$(top_srcdir)/libcdi/src
cdo_LDADD    = libcdo.la $(top_builddir)/libcdi/src/libcdi.la
cdo_LDFLAGS  =

if ENABLE_ALL_STATIC
cdo_LDFLAGS  += -all-static
endif

libcdo_la_CPPFLAGS = $(cdo_CPPFLAGS)
noinst_PROGRAMS  = cdotest
cdotest_SOURCES  = cdo_int.h	\
	           cdotest.c
cdotest_LDADD    = $(cdo_LDADD)
cdotest_CPPFLAGS = $(cdo_CPPFLAGS)
cdotest_LDFLAGS  = $(cdo_LDFLAGS)
if ENABLE_ALL_STATIC
cdotest_LDFLAGS  += -all-static
endif

#cdo-userlog.o: userlog.c config.h
#	$(COMPILE) -DLOGPATH=${exec_prefix}/log -c -o cdo-userlog.o `test -f 'userlog.c' || echo '$(srcdir)/'`userlog.c
#cdo_static-userlog.o: userlog.c config.h
#	$(COMPILE) -DLOGPATH=${exec_prefix}/log -c -o cdo_static-userlog.o `test -f 'userlog.c' || echo '$(srcdir)/'`userlog.c
#
CLEANFILES  = `ls *~ 2> /dev/null`
