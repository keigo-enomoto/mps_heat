all: a.out

a.out: param.cpp ../../src/param.hpp
	g++ -Wall -Wextra -std=c++11 param.cpp

.PHONY: clean

clean:
	rm -f a.out