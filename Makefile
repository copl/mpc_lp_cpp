UNAME_S := $(shell uname -s)
INCLUDE = -I./include
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
	
main: copl_linalg.o copl_core.o main.o copl_algorithm.o
	$(CPP) $(CFLAGS) $(INCLUDE) -std=c++11 ./bin/main.o ./bin/copl_linalg.o ./bin/copl_core.o ./bin/copl_algorithm.o \
	-o ./bin/main$(EXE_NAME).exe

main.o: 	
	$(CPP) $(CFLAGS) $(INCLUDE) -c ./src/main.cpp -o ./bin/main.o

copl_linalg.o:	
	$(CPP) $(CFLAGS) $(INCLUDE) -c ./src/copl_linalg.cpp -o ./bin/copl_linalg.o

copl_core.o:	
	$(CPP) $(CFLAGS) $(INCLUDE) -c ./src/copl_core.cpp -o ./bin/copl_core.o

copl_algorithm.o:	
	$(CPP) $(CFLAGS) $(INCLUDE) -c ./src/copl_algorithm.cpp -o ./bin/copl_algorithm.o

eigen_example:
	$(CPP) $(CFLAGS) $(INCLUDE) ./test/eigen_example.cpp -o ./bin/eigen_example.exe

copl_linalg_test: copl_linalg.o
	$(CPP) $(CFLAGS) $(INCLUDE) ./test/copl_linalg_test.cpp ./bin/copl_linalg.o -o ./bin/copl_linalg_example.exe

clean:
	cd ./bin
	rm *.o
	rm *.exe