UNAME_S := $(shell uname -s)

#include for testing 
TEST_INCLUDE = -I./external/gmock-1.7.0/gtest/include -I./external/gmock-1.7.0/include
INCLUDE = -I./include $(TEST_INCLUDE)
#Testing libraries
GTEST_DIR=./external/gmock-1.7.0/gtest
GMOCK_DIR=./external/gmock-1.7.0/

AR=ar -rv
ifeq ($(UNAME_S),Linux)
    CPP = g++
    CFLAGS = -std=c++11 -Wall -Wextra 
	EXE_NAME = _linux
endif

ifeq ($(UNAME_S),Darwin)
    CPP = clang++
    CFLAGS=-std=c++11 -stdlib=libc++ -Wall -Wextra -w
	EXE_NAME = _darwin
endif

# I use cygwin :)
# removed -pedantic because it is crazy at the moment
ifeq ($(UNAME_S),CYGWIN_NT-6.1)
	CPP = g++
	CFLAGS = -std=c++11 -Wall -Wextra -w
	EXE_NAME = _cygwin
endif

all: main

main: copl_linalg.o copl_core.o main.o copl_algorithm.o copl_newton.o
	$(CPP) $(CFLAGS) $(INCLUDE) -std=c++11 ./bin/main.o ./bin/copl_linalg.o ./bin/copl_core.o ./bin/copl_algorithm.o ./bin/copl_newton.o \
	-o ./bin/main$(EXE_NAME).exe

main.o: 	
	$(CPP) $(CFLAGS) $(INCLUDE) -c ./src/main.cpp -o ./bin/main.o

copl_linalg.o:	
	$(CPP) $(CFLAGS) $(INCLUDE) -c ./src/copl_linalg.cpp -o ./bin/copl_linalg.o

copl_newton.o:	
	$(CPP) $(CFLAGS) $(INCLUDE) -c ./src/copl_newton.cpp -o ./bin/copl_newton.o

copl_core.o:	
	$(CPP) $(CFLAGS) $(INCLUDE) -c ./src/copl_core.cpp -o ./bin/copl_core.o

copl_algorithm.o:	
	$(CPP) $(CFLAGS) $(INCLUDE) -c ./src/copl_algorithm.cpp -o ./bin/copl_algorithm.o

gmock.a: 
	$(CPP) $(CFLAGS) -isystem ${GTEST_DIR}/include -I${GTEST_DIR} \
           -isystem ${GMOCK_DIR}/include -I${GMOCK_DIR} \
           -pthread -c ${GTEST_DIR}/src/gtest_main.cc -o ./bin/gtest_main.o
	$(CPP) $(CFLAGS) -isystem ${GTEST_DIR}/include -I${GTEST_DIR} \
      	   -isystem ${GMOCK_DIR}/include -I${GMOCK_DIR} \
           -pthread -c ${GMOCK_DIR}/src/gmock-all.cc -o ./bin/gmock-all.o
	$(AR) ./bin/libgmock.a ./bin/gtest_main.o ./bin/gmock-all.o 

unittest_linalg.o: 	
	$(CPP) $(CFLAGS) $(INCLUDE) $(TEST_INCLUDE) -c ./test/unittest_linalg.cpp -o ./bin/unittest_linalg.o

test: unittest_linalg.o gmock.a
	$(CPP) $(CFLAGS) $(TEST_LIB) ./bin/unittest_linalg.o ./bin/copl_linalg.o ./bin/copl_core.o ./bin/libgmock.a  -o ./bin/unittest_linalg.exe

clean:
	cd ./bin
	rm *.o
	rm *.exe
	cd ./external/gtest-1.7.0; make clean
