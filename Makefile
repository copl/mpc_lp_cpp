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

all: stream_sampler

stream_sampler:
	$(CPP) $(CFLAGS) $(INCLUDE) ./src/eigen_example.cpp -o ./bin/eigen_example.exe
