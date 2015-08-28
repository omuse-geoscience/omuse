## Process this file with automake to produce Makefile.in
#
#TESTS = test_grib.sh
TESTS = cksum_verify \
	test_cksum_grib test_cksum_nc test_cksum_extra \
	test_cksum_service test_cksum_nc2 test_cksum_nc4 test_cksum_ieg \
	test_chunk_cksum \
	pio_write_run pio_cksum_mpinonb pio_cksum_fpguard \
	pio_cksum_asynch pio_cksum_writer pio_cksum_cdf \
	test_resource_copy pio_write_deco2d_run \
	test_f2003 test_cdf_transformation
check_PROGRAMS = cksum_verify test_grib cksum_write cksum_read pio_write \
	test_resource_copy cksum_write_chunk pio_write_deco2d

if ENABLE_NETCDF
check_PROGRAMS += test_cdf_write test_cdf_read
endif

#
test_grib_SOURCES = test_grib.c
cksum_verify_SOURCES = cksum_verify.c
cksum_write_SOURCES = cksum_write.c simple_model_helper.h simple_model_helper.c
cksum_write_chunk_SOURCES = cksum_write_chunk.c simple_model_helper.h simple_model_helper.c
cksum_read_SOURCES = cksum_read.c \
	var_cksum.c var_cksum.h \
	stream_cksum.c stream_cksum.h \
	ensure_array_size.h ensure_array_size.c
pio_write_SOURCES = pio_write.c pio_write.h simple_model.c \
	simple_model_helper.h simple_model_helper.c
pio_write_deco2d_SOURCES = pio_write.c pio_write.h deco2d_model.c \
	simple_model_helper.h simple_model_helper.c
test_resource_copy_SOURCES = test_resource_copy.c
test_resource_copy_LDADD = $(UUID_C_LIB) ../src/libcdiresunpack.la $(LDADD)
test_resource_copy_mpi_SOURCES = test_resource_copy.c
test_cdf_write_SOURCES = test_cdf_write.c
test_cdf_read_SOURCES = test_cdf_read.c
#
AM_CFLAGS = $(PPM_CORE_CFLAGS) $(YAXT_CFLAGS) $(MPI_C_INCLUDE)
if USE_MPI
pio_write_LDADD = ../src/libcdipio.la $(UUID_C_LIB) $(MPI_C_LIB)
pio_write_deco2d_LDADD = ../src/libcdipio.la $(UUID_C_LIB) $(MPI_C_LIB)
TESTS +=  test_resource_copy_mpi_run
check_PROGRAMS += test_resource_copy_mpi
test_resource_copy_mpi_LDADD = ../src/libcdipio.la $(UUID_C_LIB)
test_resource_copy_mpi_CFLAGS = $(AM_CFLAGS) $(CFLAGS) \
	-DMPI_MARSHALLING
else
pio_write_LDADD = $(LDADD) $(UUID_C_LIB)
pio_write_deco2d_LDADD = $(LDADD) $(UUID_C_LIB)
test_resource_copy_mpi_LDADD = $(LDADD) $(UUID_C_LIB)
endif

LDADD = ../src/libcdi.la -lm
AM_CPPFLAGS = -I$(top_srcdir)/src
#
EXTRA_DIST = test_cdf_write.c test_cdf_read.c
#
CLEANFILES  = `ls *~ *.grb *.nc *.srv *.ext example_*.cksum`
#
AUTOMAKE_OPTIONS = color-tests
