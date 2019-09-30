
include make.inc

all: libfishpack testfishpack

libfishpack:
	mkdir -p ./lib
	mkdir -p ./objs
	( cd ./src; $(MAKE) clean; $(MAKE) )

testfishpack:
	( cd ./test; $(MAKE) clean; $(MAKE) )

clean:
	( cd ./src; $(MAKE) clean; cd ../test; $(MAKE) clean )
