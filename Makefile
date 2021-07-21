
# C++ compiler
CC = g++

# compiler flags
#	CFLAGS = -c -I $(HOME)/include/ -std=c++17
# CFLAGS = -Wshadow -Wall -c -I $(HOME)/include/ -std=c++17
CFLAGS = -DMANIFEM_COLLECT_CM -Wshadow -Wall -c -I $(HOME)/include/ -std=c++17
#	CFLAGS = -DNDEBUG -c -O4 -I $(HOME)/include/ -std=c++17
# CFLAGS = -DNDEBUG -c -I $(HOME)/include/ -std=c++17

%.o: %.cpp
	$(CC) $(CFLAGS) $^

main-%.o: examples-manual/parag-%.cpp
	$(CC) $(CFLAGS) -o $@ $^

manifem-exe-%: iterator.o field.o finite-elem.o function.o global.o manifold.o mesh.o main-%.o progressive.o
	$(CC) -o $@ -std=c++17 $^

manifem-exe-test-%: iterator.o field.o finite-elem.o function.o global.o manifold.o mesh.o test-%.o progressive.o
	$(CC) -o $@ -std=c++17 $^

manifem-exe-12.10: main-12.10.o
	$(CC) -o $@ -std=c++17 $^

run-%: manifem-exe-%
	./$<

clean:
	rm *.o manifem-exe-*

clean-all: clean
	rm *.eps *.msh

.SECONDARY:

.PHONY: clean clean-all

