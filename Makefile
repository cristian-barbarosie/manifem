

%.o: %.cpp
#	g++ -c -I $(HOME)/include/ -std=c++17 $^
	g++ -Wshadow -Wall -c -I $(HOME)/include/ -std=c++17 $^
#	g++ -DNDEBUG -c -I $(HOME)/include/ -std=c++17 $^

manifem-exe-%: field.o finite-elem.o function.o global.o iterator.o manifold.o mesh.o main-%.o progressive.o
	g++ -o $@ -std=c++17 $^

manifem-exe-9.15: main-9.15.o
	g++ -o $@ -std=c++17 $^

run-%: manifem-exe-%
	./$<

run-test: manifem-exe-test
	./$<

manifem-exe-test: test.o finite-elem.o function.o mesh.o field.o global.o iterator.o manifold.o mesh.o progressive.o
	g++ -o $@ -std=c++17 $^

clean:
	rm *.o manifem-exe-*

clean-all: clean
	rm *.eps *.msh

.SECONDARY:

.PHONY: clean clean-all run-test run-%

