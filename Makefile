UNAME_S := $(shell uname -s)
INCLUDE = -I./include
ifeq ($(UNAME_S),Linux)
    CPP = g++
    CFLAGS = -std=c++11 -Wall -Wextra -pedantic
endif

ifeq ($(UNAME_S),Darwin)
    CPP = clang++
    CFLAGS = -std=c++11 -stdlib=libc++ -Wall -Wextra -pedantic
endif

# I use cygwin :)
ifeq ($(UNAME_S),CYGWIN_NT-6.1-WOW)
	CPP = g++
	CFLAGS = -std=c++11 -Wall -Wextra -pedantic
endif

all: main
	
main:
	$(CPP) $(CFLAGS) $(INCLUDE) -std=c++11 ./src/main.cpp -o ./bin/main.exe
	

copl_linalg.o:	
	$(CPP) $(CFLAGS) $(INCLUDE) -c ./src/copl_linalg.cpp -o ./bin/copl_linalg.o

eigen_example:
	$(CPP) $(CFLAGS) $(INCLUDE) ./test/eigen_example.cpp -o ./bin/eigen_example.exe

copl_linalg_test: copl_linalg.o
	$(CPP) $(CFLAGS) $(INCLUDE) ./test/copl_linalg_test.cpp ./bin/copl_linalg.o -o ./bin/copl_linalg_example.exe

clean:
	cd ./bin
	rm *.o
	rm *.exe