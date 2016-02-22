Implementation of Single-level and Multilevel sparsification. See paper here: http://arxiv.org/abs/1601.05527
Note: This project is still in development.

To build:

	make

To run:

For single-level:
	./msparse -i ~/test_graphs/fb-uf.edges -n 35111 --zero-index --gname fb-uf -o out.edges -s 0.4 --single-level --weak

For multi-level
	./msparse -i ~/test_graphs/fb-uf.edges -n 35111 --zero-index --gname fb-uf -o out.edges --param-list 0.1,0.1,0.1 --weak

To set number of openmp threads (set to 8 by default)
	export OMP_NUM_THREADS=16

To disable openmp, remove '-fopenmp' from LFLAGS in Makefile

Notes:
======
==> Param list starts from coarsest level

==> c++11 required (version g++ 4.8 or over)

==> Graph file must already be undirected, node numbering has to be continuous and need only contain one direction of relationship e.g

	Good graph file
		1 2
		3 4
		5 3

	Bad graph file
		1 2
		2 1
		3 4
		4 3
		5 3
		3 5
	graph files (edge list) can also contain 3rd column for weights.

==> Space delimeter is used by default. To change to tab use --delim tab

==> Weak + Strong edges sparsified by default.

Options
=======

--weak: retain weak edges
--strong: Retain strong edges
--normalize: normalize algebraic distance by product of degrees.

-i - path of input graph
-n - number of nodes
--gname - shortname for graph.
--zero-index - index of first node starts at 0.
--permute - random shuffle node ordering
-r -  interpolation order. Default 1.
--delim - delimeter of input file. Only space or tab allowed.

 See Main.cpp for all available options.

 Output:

 sparsified graph
