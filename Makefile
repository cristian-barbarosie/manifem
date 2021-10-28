
# C++ compiler
CC = g++

# compiler flags
#	CFLAGS = -c -I . -I $(HOME)/include/ -std=c++17
CFLAGS = -Wshadow -Wall -c -I . -I $(HOME)/include/ -std=c++17
# CFLAGS = -DMANIFEM_COLLECT_CM -Wshadow -Wall -c -I . -I $(HOME)/include/ -std=c++17
# CFLAGS = -DMANIFEM_COLLECT_CM -DNDEBUG -O4 -c -I . -I $(HOME)/include/ -std=c++17
# CFLAGS = -DNDEBUG -c -O4 -I . -I $(HOME)/include/ -std=c++17
# CFLAGS = -DNDEBUG -c -I . -I $(HOME)/include/ -std=c++17

manifem_objects = iterator.o field.o finite-elem.o function.o global.o manifold.o mesh.o progressive.o

%.o: %.cpp
	$(CC) $(CFLAGS) $^

parag-%.o: examples-manual/parag-%.cpp
	$(CC) $(CFLAGS) -o $@ $^

manifem-exe-%: parag-%.o $(manifem_objects)
	$(CC) -o $@ -std=c++17 $^

# manifem-exe-test-%: test-%.o $(test_objects)
# 	$(CC) -o $@ -std=c++17 $^

parag-12.10.o: parag-12.10.cpp
	$(CC) $(CFLAGS) -o $@ $^

manifem-exe-12.10: parag-12.10.o
	$(CC) -o $@ -std=c++17 $^

# manifem-exe-test-12.10: test-12.10.o
# 	$(CC) -o $@ -std=c++17 $^

run-%: manifem-exe-%
	./$<

clean:
	rm *.o manifem-exe-*

clean-all: clean
	rm *.eps *.msh

.SECONDARY:

.PHONY: clean clean-all

