# maniFEM
ManiFEM is a C++ library for solving partial differential equations
through the finite element method.
The name comes from "finite elements on manifolds". 
ManiFEM was designed with the goal of coping with very general meshes,
in particular meshes on Riemannian manifolds, even manifolds which cannot be embedded in R^3,
like the torus R^2/Z^2.
Also, maniFEM was written with the goal of being conceptually clear and easy to read.
We hope it will be particularly useful for people who want fine control over the mesh, 
e.g. for implementing their own meshing or remeshing algorithms.

ManiFEM is just a collection of C++ classes.
It has no user-friendly interface nor graphic capabilities. 
The user should have some understanding of programming and of C++. 
However, maniFEM can be used at a basic level by people with no deep knowledge of C++.

Similar products (competitors) : [FreeFem++](http://www3.freefem.org/),
[FEniCS/Dolphin](https://fenicsproject.org/), [deal.II](https://dealii.org/)

In its current version, release 21.11, maniFEM works well for mesh generation,
including meshes on quotient manifolds. 
Anisotropic Riemann metrics (paragraph 3.24 in the manual) and
variational formulations are not yet implemented. 
Finite elements (section 6 in the manual) are implemented in a rather rudimentary manner for now.
A changelog is available at the end of the
[manual](http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf).
To check which version of maniFEM is installed in your computer, see at the beginning
of the file `maniFEM.h`.

ManiFEM is being developed by [Cristian Barbarosie](mailto:cristian.barbarosie@gmail.com),
[Sérgio Lopes](mailto:slopes@adm.isel.pt) and
[Anca-Maria Toader](mailto:anca.maria.toader@gmail.com);
see its [homepage](http://manifem.rd.ciencias.ulisboa.pt).

To learn maniFEM, you should read the
[manual](http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf) (version 21.11).

To use maniFEM, choose a [release](https://github.com/cristian-barbarosie/manifem/releases)
and download all files to some directory in your computer.
Latest code might be unstable; releases are stable.
To check which version of maniFEM is installed in your computer,
see at the beginning of the file `maniFEM.h`.
You can then run the examples in the manual :
just `make run-1.1` for the example in paragraph 1.1, 
`make run-2.6` for the example in paragraph 2.6, and so on.
Paragraph 11.16 in the [manual](http://manifem.rd.ciencias.ulisboa.pt/manual-manifem.pdf)
gives more details.
You will need a recent C++ compiler (we use `g++`) and the `make` utility. 
Under linux it should be easy to install them. 
It is not that easy to install and use them under Windows, but it is certainly possible,
for instance by using [cygwin](https://cygwin.org).
Some examples require the [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) library; 
just copy its source tree somewhere in your computer and be sure that path is mentioned in your 
`Makefile` under the `-I` flag of your compiler.
You may want to use [gmsh](http://gmsh.info/) for visualization purposes. 

A component of maniFEM, [MetricTree](https://github.com/cristian-barbarosie/MetricTree),
can be used independently.

This work is supported by National Funding from FCT - Fundação para a Ciência e a Tecnologia
(Portugal), through Faculdade de Ciências da Universidade de Lisboa and 
Centro de Matemática, Aplicações Fundamentais e Investigação Operacional,
project UID/MAT/04561/2020.

Copyright 2019, 2020, 2021 Cristian Barbarosie cristian.barbarosie@gmail.com

ManiFEM is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ManiFEM is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

Full text of the GNU Lesser General Public License is available 
in files [COPYING](src/COPYING) and [COPYING.LESSER](src/COPYING.LESSER).
It can also be found at <https://www.gnu.org/licenses/>.
