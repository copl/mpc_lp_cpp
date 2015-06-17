UNAME_S := $(shell uname -s)

#include for testing 
TEST_INCLUDE = -I./external/gmock-1.7.0/gtest/include -I./external/gmock-1.7.0/include
INCLUDE = -I./include $(TEST_INCLUDE)
#Testing libraries
GTEST_DIR=./external/gmock-1.7.0/gtest
GMOCK_DIR=./external/gmock-1.7.0/

AR=ar -rv

COPL_CPP_FILES := $(wildcard src/*.cpp) 
#remove src/main.cpp 
COPL_CPP_FILES := $(filter-out src/main.cpp, $(COPL_CPP_FILES))
TEST_CPP_FILES := $(wildcard test/*.cpp) 
COPL_OBJ_FILES := $(addprefix bin/,$(notdir $(COPL_CPP_FILES:.cpp=.o)))
TEST_OBJ_FILES := $(addprefix testbin/,$(notdir $(TEST_CPP_FILES:.cpp=.o)))

ifeq ($(UNAME_S),Linux)
    CPP = g++
    CFLAGS = -std=c++11 -Wall -Wextra -fPIC
	EXE_NAME = _linux
	SHARED_EXTENSION = so
endif

ifeq ($(UNAME_S),Darwin)
    CPP = clang++
    CFLAGS = -std=c++11 -stdlib=libc++ -Wall -Wextra -w -g -fPIC
	EXE_NAME = _darwin
	SHARED_EXTENSION = dylib
endif

# I use cygwin :)
# removed -pedantic because it is crazy at the moment
ifeq ($(UNAME_S),CYGWIN_NT-6.1)
	CPP = g++
	CFLAGS = -std=c++11 -Wall -Wextra -w -fPIC
	EXE_NAME = _cygwin
	SHARED_EXTENSION = dll
endif

all: main

#We treat main differently to all the other files so that when we test we 
#dont link two files with a main function
main.o:	
	$(CPP) $(CFLAGS) $(INCLUDE) $(TEST_INCLUDE) -c ./src/main.cpp -o ./bin/main.o

main: $(COPL_OBJ_FILES) main.o
	mkdir -p ./bin
	$(CPP) -stdlib=libc++ -std=c++11 $(COPL_OBJ_FILES) ./bin/main.o -o ./bin/main$(EXE_NAME).exe

#compile gtest and gmock
bin/libgmock.a: 
	$(CPP) $(CFLAGS) -isystem ${GTEST_DIR}/include -I${GTEST_DIR} \
           -isystem ${GMOCK_DIR}/include -I${GMOCK_DIR} \
           -pthread -c ${GTEST_DIR}/src/gtest_main.cc -o ./bin/gtest_main.o
	$(CPP) $(CFLAGS) -isystem ${GTEST_DIR}/include -I${GTEST_DIR} \
           -isystem ${GMOCK_DIR}/include -I${GMOCK_DIR} \
           -pthread -c ${GTEST_DIR}/src/gtest-all.cc -o ./bin/gtest-all.o           
	$(CPP) $(CFLAGS) -isystem ${GTEST_DIR}/include -I${GTEST_DIR} \
      	   -isystem ${GMOCK_DIR}/include -I${GMOCK_DIR} \
           -pthread -c ${GMOCK_DIR}/src/gmock-all.cc -o ./bin/gmock-all.o
	$(AR) ./bin/libgmock.a ./bin/gtest_main.o ./bin/gtest-all.o ./bin/gmock-all.o 

#Compile all test objects
testbin/%.o: test/%.cpp
	$(CPP) $(CFLAGS) $(INCLUDE) $(TEST_INCLUDE) -c $< -o $@

#Compile all algorithm objects
bin/%.o: src/%.cpp 	
	$(CPP) $(CFLAGS) $(INCLUDE) $(TEST_INCLUDE) -c $< -o $@

lib: $(COPL_OBJ_FILES)
	mkdir -p ./lib
	$(CPP)  -stdlib=libc++ -std=c++11 $(COPL_OBJ_FILES) -o ./lib/copl_sovler.$(SHARED_EXTENSION) --shared -fPIC

test: $(TEST_OBJ_FILES) $(COPL_OBJ_FILES) bin/libgmock.a
	$(CPP) $(CFLAGS) $(TEST_OBJ_FILES) $(COPL_OBJ_FILES) bin/libgmock.a -o ./bin/unittest.exe

clean:
	rm ./bin/*; rm ./testbin/*; rm ./lib/*
