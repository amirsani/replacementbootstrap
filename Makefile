# ----------------------------------------------------------------------------
# plain simple rules to make and cleanup the library:
# make default;   compiles the library
# make test;      compiles and executes test. O.K. message marks success.
# make clean;     removes temporary files
# make cleanall;  removes temporary files, the library, and programs
# ----------------------------------------------------------------------------

CXX      = g++
CPPFLAGS = -O3 -I. -I/usr/local/include -std=c++0x -Wall
LDFLAGS  = -L. -lgzstream -lz
LDMPFLAGS = -fopenmp
AR       = ar cr

default: main

functions.o : functions.cpp functions.h
	${CXX} ${CPPFLAGS} -c -o functions.o functions.cpp

bootstrap_methods.o : bootstrap_methods.cpp bootstrap_methods.h functions.h
	${CXX} ${CPPFLAGS} -c -o bootstrap_methods.o bootstrap_methods.cpp ${LDFLAGS}

main.o : main.cpp bootstrap_methods.h functions.h
	${CXX} ${CPPFLAGS} -c -o main.o main.cpp

main : main.o bootstrap_methods.o functions.o
	${CXX} ${CPPFLAGS} -o main main.o bootstrap_methods.o functions.o ${LDFLAGS}

gzstream.o : gzstream.C gzstream.h
	${CXX} ${CPPFLAGS} -c -o gzstream.o gzstream.C

libgzstream.a : gzstream.o
	${AR} libgzstream.a gzstream.o

clean :
	rm *.o gen_ryabko_process main functions ryabko_process

cleanall :
	rm *.o libgzstream.a test_gzip test_gunzip

# ============================================================================
# EOF
