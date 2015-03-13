UNAME_S := $(shell uname -s)

#include for testing 
TEST_INCLUDE = -I./external/gmock-1.7.0/gtest/include -I./external/gmock-1.7.0/include
INCLUDE = -I./include $(TEST_INCLUDE)
#Testing libraries
GTEST_DIR=./external/gmock-1.7.0/gtest
GMOCK_DIR=./external/gmock-1.7.0/

AR=ar -rv

COPL_CPP_FILES := $(wildcard src/*.cpp) 
TEST_CPP_FILES := $(wildcard test/*.cpp) 
COPL_OBJ_FILES := $(addprefix bin/,$(notdir $(COPL_CPP_FILES:.cpp=.o)))
TEST_OBJ_FILES := $(addprefix testbin/,$(notdir $(TEST_CPP_FILES:.cpp=.o)))

ifeq ($(UNAME_S),Linux)
    CPP = g++
    CFLAGS = -std=c++11 -Wall -Wextra 
	EXE_NAME = _linux
endif

ifeq ($(UNAME_S),Darwin)
    CPP = clang++
    CFLAGS = -std=c++11 -stdlib=libc++ -Wall -Wextra -w
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

main: $(COPL_OBJ_FILES)
	$(CPP) -stdlib=libc++ -std=c++11 $(COPL_OBJ_FILES) -o ./bin/main$(EXE_NAME).exe

#compile gtest and gmock
gmock.a: 
	$(CPP) $(CFLAGS) -isystem ${GTEST_DIR}/include -I${GTEST_DIR} \
           -isystem ${GMOCK_DIR}/include -I${GMOCK_DIR} \
           -pthread -c ${GTEST_DIR}/src/gtest_main.cc -o ./bin/gtest_main.o
	$(CPP) $(CFLAGS) -isystem ${GTEST_DIR}/include -I${GTEST_DIR} \
      	   -isystem ${GMOCK_DIR}/include -I${GMOCK_DIR} \
           -pthread -c ${GMOCK_DIR}/src/gmock-all.cc -o ./bin/gmock-all.o
	$(AR) ./bin/libgmock.a ./bin/gtest_main.o ./bin/gmock-all.o 

#Compile all test objects
testbin/%.o: test/%.cpp
	$(CPP) $(CFLAGS) $(INCLUDE) $(TEST_INCLUDE) -c $< -o $@

#Compile all algorithm objects
bin/%.o: src/%.cpp 	
	$(CPP) $(CFLAGS) $(INCLUDE) $(TEST_INCLUDE) -c $< -o $@

test: $(TEST_OBJ_FILES) $(COPL_OBJ_FILES) gmock.a
	$(CPP) $(CFLAGS) $(TEST_LIB) ./bin/libgmock.a $(TEST_OBJ_FILES) $(COPL_OBJ_FILES) -o ./bin/unittest.exe

help:
	echo $(COPL_CPP_FILES)	
	echo $(COPL_OBJ_FILES)

clean:
	rm ./bin/*.o	
	rm ./testbin/*.o
	rm ./bin/*.exe
	cd ./external/gtest-1.7.0; make clean
