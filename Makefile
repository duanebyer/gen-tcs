CC              = g++ -std=c++11 -Wall
CC_OBJ_FLAGS    = -c -fPIC
CC_Shared_FLAGS = -shared -Wl,-soname,libTCSGen.so
ROOT_CFLAGS     = $(shell ${ROOTSYS}/bin/root-config --cflags)
ROOT_LIBS       = $(shell ${ROOTSYS}/bin/root-config --libs)
libTCSGEN	= libTCSGen

all:	    TCSGen.cc TTCSKine.o KinFuncs.o CrsFuncs.o TTCSCrs.o GPDs.o
	    mkdir -p lib ; rm -f lib/*.so
	    $(CC) $(CC_Shared_FLAGS) -o lib/${libTCSGEN}.so.1.0.1 TTCSKine.o KinFuncs.o CrsFuncs.o TTCSCrs.o GPDs.o
	    cd lib;\
	    ln -sf ${libTCSGEN}.so.1.0.1 ${libTCSGEN}.so.1; ln -sf ${libTCSGEN}.so.1.0.1 ${libTCSGEN}.so
	    cd ../;
	    $(CC) -o TCSGen.exe TCSGen.cc -I ./include -L./lib -lTCSGen $(ROOT_CFLAGS) $(ROOT_LIBS)
	
TTCSKine.o: src/TTCSKine.cc include/TTCSKine.h
	    $(CC) $(CC_OBJ_FLAGS) src/TTCSKine.cc -o $@ $(ROOT_CFLAGS) -I ./include
	
GPDs.o:	    src/GPDs.cc include/GPDs.h
	    $(CC) $(CC_OBJ_FLAGS) src/GPDs.cc -o $@ $(ROOT_CFLAGS) -I ./include

TTCSCrs.o:  src/TTCSCrs.cc include/TTCSCrs.h GPDs.o
	    $(CC) $(CC_OBJ_FLAGS) src/TTCSCrs.cc -o $@ $(ROOT_CFLAGS) -I ./include 
	
KinFuncs.o: src/KinFunctions.cc include/KinFunctions.h
	    $(CC) $(CC_OBJ_FLAGS) src/KinFunctions.cc -o $@ -I ./include

CrsFuncs.o: src/CrsFunctions.cc include/CrsFunctions.h
	    $(CC) $(CC_OBJ_FLAGS) src/CrsFunctions.cc -o $@ -I ./include



clean:	    
	    rm -f TCSGen.exe *.o lib/*.so.* lib/*.so
