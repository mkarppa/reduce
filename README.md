# reduce
Reduce - A tool for symmetry breaking

## Overview

This git repository contains the C++ source code for `reduce', an
experimental software implementation of adaptive prefix-assignment
symmetry reduction; cf.

  T. Junttila, M. Karppa, P. Kaski, J. Kohonen,
  "An adaptive prefix-assignment technique for symmetry reduction".

This experimental software is supplied to accompany the aforementioned
manuscript.

##Licence

The source code is subject to the MIT Licence, see 'LICENSE' for details.

##Requirements

The software requires the following:
- UNIX-like environment (GNU/Linux and MacOS X should work)
- Sufficiently recent C++ compiler (tested to work with CLang 9.0.0,
G++ 5.4.0, and G++ 6.4.0) and other relevant build tools (especially make)
- nauty (version 2.6r10 should work)
- GMP (version 6.1.2 should work)
- TCLAP (version 1.2.1 should work)
- For parallel execution, an MPI implementation is required
(recommended: OpenMPI)

##Build instructions

1. Download 'nauty' from http://pallini.di.uniroma1.it/ and build it
under nauty/ following instructions from the website (usually simply
./configure && make)

2. Download 'TCLAP' from http://tclap.sourceforge.net/ and build it
under 'tclap' following instructions from the website (usually simply
./configure && make)

3. Download 'GMP' from https://gmplib.org/ and build it under 'gmp'
Note: use '--enable-cxx' flag when configuring. On MacOS X, you should
also enable '--with-pic'.

4. If you do not plan to use MPI, edit the Makefile and replace
   'CXX=mpic++' with 'CXX=c++' (simply uncomment the lower line and
   comment out the first line); also remove '-DWITH_OPEN_MPI' from
   CXXFLAGS. If you have OpenMPI installed, you can simply ignore this
   step.

5. If the library versions differ from those assumed, please replace
   the corresponding path variables in Makefile with correct
   versions. We apologize for this inconvenience!

6. Run 'make'

## Testing
You may use the provided files 'A000088_8.g', 'mm2227.cnf+g',
'mm2227.cnf', and 'mm2227.g' for testing as follows.

```
$ ./reduce -ng A000088_8.g -p "1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
17 18 19 20 21 22 23 24 25 26 27 28"
$ ./reduce -g -p "1 2 3 4 5 6 7" mm2227.cnf+g
$ ./reduce -p "1 2 3 4 5" bpm65.cnf
```

The first one should produce 12346 partial assignments. The second
should provide a satisfiable CNF instance with 8 partial
assignments. The last one should provide an unsatisfiable CNF instance
with 6 partial assignments.

## Usage

The command
```
$ ./reduce -h
```
prints an usage display of the available command-line options. 

### Usage without an explicit symmetry graph

By default, 'reduce' expects as input a CNF formula in DIMACS format 
and a prefix supplied from the command line. Input can be read from
standard input by specifying '-' as the input file.

### Usage with CNF and an explicit symmetry graph

In most cases one wants to use 'reduce' so that the input symmetries
are described by an explicit graph. When the '-g' switch is specified,
the graph is expected to be appended at the end of the CNF input file.

Here is an example that contains a system of clauses in CNF, followed by a graph
representation of the symmetry in the system:

```
p cnf 6 3
1 2 0
1 3 5 0
2 4 6 0
p edge 8 4
e 1 3
e 1 5
e 2 4
e 2 6
c 1 1
c 2 1
c 3 1
c 4 1
c 5 1
c 6 1
c 7 2
c 8 3
p variable 6
v 1 1
v 2 2
v 3 3
v 4 4
v 5 5
v 6 6
p value 2
r 7 false
r 8 true
```

Let us look at the description of the symmetry graph in parts.
First, we describe the vertices and edges of the graph:
```
p edge 8 4
e 1 3
e 1 5
e 2 4
e 2 6
```
That is, the graph has 8 vertices and 4 edges. Each edge is described by
a unique 'e'-line of the format 'e &lt;i&gt; &lt;j&gt;' for an edge joining the
vertex &lt;i&gt; to the vertex &lt;j&gt;. The vertices are numbered from 1 to N, where
N is the number of vertices in the graph. Next, we describe the colors
of the vertices using 'c'-lines:
```
c 1 1
c 2 1
c 3 1
c 4 1
c 5 1
c 6 1
c 7 2
c 8 3
```
The format of a 'c'-line is 'c &lt;i&gt; &lt;k&gt;' that indicates that the vertex &lt;i&gt;
has color &lt;k&gt;. The color &lt;k&gt; must be a positive integer. 

Finally, we define which vertices of the symmetry graph correspond to which
variables and which values of the CNF instance. This is done using 
'v'-lines and 'r'-lines as follows:
```
p variable 6
v 1 1
v 2 2
v 3 3
v 4 4
v 5 5
v 6 6
p value 2
r 7 false
r 8 true
```
These lines declare that there are six variable-vertices in the symmetry
graph, corresponding to the six variables in the CNF instance. Each 'v'-line
gives one such correspondence, with the format 'v &lt;i&gt; &lt;q&gt;' indicating that
the vertex &lt;i&gt; in the graph corresponds to the variable &lt;q&gt; in the CNF 
instance above. Furthermore, there are two value-vertices in the graph,
each declared with an 'r'-line, with the format 'r &lt;i&gt; &lt;v&gt;' indicating
that the vertex &lt;i&gt; in the graph corresponds to the value &lt;v&gt; in the CNF
instance. Both values 'false' and 'true' must be present.

### Usage with the symmetry graph only

We can also instruct 'reduce' not to expect CNF input at all by using
the '-n' switch (for 'no CNF'). In that case, symmetry graph should be
provided along with the '-g' switch.

### Parallel execution

Parallel execution differs from sequential execution only by the fact
that reduce must be run using 'mpiexec' (or similar). For further
details, please refer to the documentation of your MPI implementation.
