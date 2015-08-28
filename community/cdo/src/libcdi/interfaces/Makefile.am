## Process this file with automake to produce Makefile.in
#
SWIGCDILIB = cdilib.i
SWIGCDIOBJ = cdiobj.i
BINDINGS   =

EXTRA_DIST = $(SWIGCDIOBJ) $(SWIGCDILIB) \
	\
	ruby/cdilib_wrap.c ruby/cdiobj_wrap.cpp ruby/testLib.rb \
	ruby/testObj.rb \
	\
	python/cdilib_wrap.c python/cdiobj_wrap.cpp python/testLib.py \
	python/testObj.py python/setupObj.py python/setupLib.py \
	python/setup.py python/table.py python/CdiLib.py \
	python/CdiObj.py \
	\
	CdiInfo.cpp \
	\
	f2003/bindGen.rb

if CREATE_INTERFACES
BINDINGS += help
endif
#
bin_PROGRAMS =
noinst_PROGRAMS = 
noinst_LTLIBRARIES = 
if CREATE_INTERFACES
noinst_PROGRAMS += CdiInfo
noinst_LTLIBRARIES += libcdipp.la
endif
#
locallibs = $(top_builddir)/src/.libs/libcdi.so @LTLIBOBJS@ @NETCDF_LIBS@ @HDF5_LIBS@ @SZLIB_LIBS@ @ZLIB_LIBS@

libcdipp_la_SOURCES = cdi.cpp cdi.hpp
libcdipp_la_LIBADD  = $(top_builddir)/src/libcdi.la
libcdipp_la_LDFLAGS = @NETCDF_LIBS@ @HDF5_LIBS@ @SZLIB_LIBS@ @ZLIB_LIBS@ $(LDFLAGS)

CdiInfo_SOURCES     = CdiInfo.cpp
CdiInfo_LDADD       = libcdipp.la
CdiInfo_LDFLAGS     = @NETCDF_LIBS@ @HDF5_LIBS@ @SZLIB_LIBS@ @ZLIB_LIBS@ $(LDFLAGS)
#
#
if ENABLE_RUBY
BINDINGS += rubyObj
# Ruby ====================================================
rubyLibMakefile:
	cd ruby; CFLAGS="$(CFLAGS)" LIBS="$(LIBS)" LDFLAGS="$(locallibs) $(LDFLAGS)" INCFLAGS="-I../../src" ruby extconfLib.rb

rubyLibWrapper: cdilib.i
if ENABLE_SWIG
	$(SWIG) -ruby -module CdiLib -o ruby/cdilib_wrap.c $(AM_CPPFLAGS) cdilib.i
endif

rubyLib: rubyLibWrapper rubyLibMakefile
	@cd ruby; make

rubyLibTest: rubyLib
	@cd ruby; ruby testLib.rb
# =========================================================
rubyObjMakefile:
	cd ruby; CFLAGS="$(CFLAGS)" LIBS="$(LIBS)" LDFLAGS="$(LDFLAGS)" INCFLAGS="-I../" ruby extconfObj.rb

rubyObjWrapper: cdiobj.i
if ENABLE_SWIG
	$(SWIG) -ruby -c++ -module CdiObj -o ruby/cdiobj_wrap.cpp $(AM_CPPFLAGS) cdiobj.i
endif

rubyObj: rubyObjWrapper rubyObjMakefile
	@cd ruby; make

rubyObjTest: rubyObj
	@cd ruby; ruby testObj.rb

rubyObjInstall: rubyObj
	@cd ruby; make install
endif
# =========================================================
if ENABLE_PYTHON
BINDINGS += pythonObj
# Python ==================================================
pythonLibWrapper: cdilib.i
if ENABLE_SWIG
	$(SWIG) -python -module CdiLib -o python/cdilib_wrap.c $(AM_CPPFLAGS) cdilib.i
endif

pythonLib: pythonLibWrapper
	@cd python; CFLAGS="$(CFLAGS)" LIBS="$(LIBS)" LDFLAGS="$(LDFLAGS)" INCFLAGS="-I../../src" python setupLib.py build_ext --inplace

pythonLibTest: pythonLib
	@cd python; python testLib.py
# =========================================================
python/cdiobj_wrap.cpp: cdiobj.i
if ENABLE_SWIG
	$(SWIG) -python -c++ -module CdiObj -o python/cdiobj_wrap.cpp $(AM_CPPFLAGS) cdiobj.i
endif

pythonObj: python/cdiobj_wrap.cpp
	cd python; CC="$(CC)" CXX="$(CXX)" BUILDLIBDIR="$(top_builddir)/src/.libs" LIBDIR=$(libdir) CFLAGS="$(CFLAGS)" LIBS="$(LIBS)" LDFLAGS="-L../../src/.libs $(LDFLAGS)" INCFLAGS="-I../" python setupObj.py build_ext --inplace

pythonObjTest: pythonObj
	@cd python; python testObj.py

pythonObjInstall: pythonObj
	cd python; CC="$(CC)" CXX="$(CXX)" BUILDLIBDIR="$(top_builddir)/src/.libs" LIBDIR=$(libdir) CFLAGS="$(CFLAGS)" LIBS="$(LIBS)" LDFLAGS="-L../../src/.libs $(LDFLAGS)" INCFLAGS="-I../" python setupObj.py install --prefix=$(prefix)
endif
# =========================================================
help:
	@echo "#=====================================================#"
	@echo "# EXPERIMENTAL FEATURE ===============================#"
	@echo "#=====================================================#"
	@echo "Bindings for: python ruby"
	@echo "(please exchange <lang> with the corresponding language)"
	@echo "#=====================================================#"
	@echo "  help             - this message"
	@echo "  <lang>Obj        - create object interface for Ruby"
	@echo "  <lang>Lib        - create function interface for Ruby"
	@echo "  <lang>ObjWrapper - create the wrapper files using swig (object interface)"
	@echo "  <lang>LibWrapper - create the wrapper files using swig (function interface)"
	@echo "  <lang>ObjTest    - run object interface test (./<lang>/testObj.rb). Accepts input file argument."
	@echo "  <lang>LibTest    - run function library test (./<lang>/testLib.rb). Accepts input file argument."
	@echo "  <lang>ObjInstall - install shared libraries and optional language files for <lang>"
	@echo "  <lang>LibInstall - see <lang>ObjInstall"
#
AM_CPPFLAGS = -I$(top_srcdir)/src
#
#all-local: $(BINDINGS)
#
CLEANFILES  = `ls *~`
