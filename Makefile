
# C++ compiler
CC = g++

# compiler flags
# CFLAGS = -c -I . -I src/ -I $(HOME)/include/ -std=c++17
CFLAGS = -c -Wshadow -Wall -I . -I src/ -I $(HOME)/include/ -std=c++17
# CFLAGS = -c -DMANIFEM_COLLECT_CM -Wshadow -Wall -I . -I src/ -I $(HOME)/include/ -std=c++17
# CFLAGS = -c -DNDEBUG -O4 -I . -I src/ -I $(HOME)/include/ -std=c++17
# CFLAGS = -c -DMANIFEM_COLLECT_CM -DNDEBUG -O4 -I . -I src/ -I $(HOME)/include/ -std=c++17
# CFLAGS = -c -DNDEBUG -I . -I src/ -I $(HOME)/include/ -std=c++17

manifem_objects = iterator.o field.o finite-elem.o frontal.o function.o global.o manifold.o mesh.o

%.o: src/%.cpp
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

manifem-exe-user-%: user-%.o $(manifem_objects)
	$(CC) -o $@ -std=c++17 $^

user-%.o: user-%.cpp
	$(CC) $(CFLAGS) -o $@ $^

run-%: manifem-exe-%
	./$<

static-lib: $(manifem_objects)
	ar cr libmaniFEM.a $^

clean:
	rm -f *.o manifem-exe-* libmaniFEM.a

clean-all: clean
	rm -f *.eps *.msh libmaniFEM*.a

.SECONDARY:

.PHONY: clean clean-all

