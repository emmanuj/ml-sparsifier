#include "solver.h"

#include <iostream>
#include "node.h"
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include "globals.h"
#include "OptionParser.h"

int main(int argc, char * argv[]) {
	time_t rtime = time(0);
	srand (rtime);

	using namespace optparse;

	const std::string usage ="Usage: usage...";//TODO: explain usage
	OptionParser parser = OptionParser()
	.usage(usage);

	parser.set_defaults("k1", "6");
	parser.set_defaults("k2", "6");
	parser.set_defaults("k3", "4");

	parser.set_defaults("mu", "2");
	parser.set_defaults("q", "0.4");
	parser.set_defaults("a", "0.4");
	parser.set_defaults("s", "0.1"); //size of graph sample. 50% by default
	parser.set_defaults("g", "1");
	parser.set_defaults("r", "2");
	parser.set_defaults("c", "20");//#clusters
	parser.set_defaults("t", "2");
	parser.set_defaults("b","0.01");
	parser.set_defaults("f", "edgelist");
	parser.set_defaults("winsize", "20");
	parser.set_defaults("permute", "0");
	parser.set_defaults("debug", "0");
	parser.set_defaults("log-info", "0");
	parser.set_defaults("zero-Index", "0");
	parser.set_defaults("weak","0");
	parser.set_defaults("strong","0");
	parser.set_defaults("normalize","0");
	parser.set_defaults("level-span","3");
	parser.set_defaults("sparse-level","0");

	parser.set_defaults("d", "space");

	parser.add_option("--zero-index") .action("store_true").dest("zero-Index") .help("Indicate whether node ID's are zero indexed");
	parser.add_option("--single-level") .action("store_true").dest("single-level") .help("Sparsify only at the finest level");
	parser.add_option("--permute") .action("store_true").dest("permute").help("Permute initial node ordering");
	parser.add_option("--weak").action("store_true").dest("weak").help("Prefer weak edges in sparsification");
	parser.add_option("--strong").action("store_true").dest("strong").help("Prefer strong edges in sparsification");
	parser.add_option("--normalize").action("store_true").dest("normalize").help("Normalize algebraic distance computation");
	parser.add_option("--debug") .action("store_true").dest("debug").help("Print debugging information");
	parser.add_option("--log-info") .action("store_true").dest("debug").help("Print log messages");
	parser.add_option("--digraph").action("store_true").dest("digraph").help("input graph is directed");
	parser.add_option("--gname").dest("gname").type("string").help("short name for this graph");

	parser.add_option("--param-list").dest("param-list").type("string").help("list of params");

	parser.add_option("-i","--input").dest("i").type("string").help("Input file");
	parser.add_option("-o","--output").dest("o").type("string").help("Output file");
	parser.add_option("-d","--delimeter").dest("d").type("string").help("Delimeter of the input graph. Default is space");
	parser.add_option("-f","--format").dest("f").type("string").help("Graph file format( e.g edgelist or graph)");
	parser.add_option("-n","--num-nodes").dest("n").type("long").help("number of nodes in graph");
	parser.add_option("-c","--clust").dest("c").type("int").help("number of clusters to find in the graph");
	parser.add_option("-s","--sample").dest("s").type("double").help("Sampling parameter");

	parser.add_option("--level-span").dest("level-span").type("double").help("Level span of msparse parameter");
	parser.add_option("--sparse-level").dest("sparse-level").type("double").help("specify which levels to sparsify");


	parser.add_option("--k1").dest("k1").type("int").help("#Compatible relaxations");
	parser.add_option("--k2").dest("k2").type("int").help("#Gauss-Seidel relaxations");
	parser.add_option("--k3").dest("k3").type("int").help("# of minimization sweeps");
	parser.add_option("-m", "--mu").type("int") .dest("mu").help("be verbose!");
	parser.add_option("-w", "--winsize").type("int") .dest("winsize").help("Size of the window");
	parser.add_option("-q").type("double").dest("q") .help("be silent!");
	parser.add_option("-a").type("double").dest("a").help("number of files (default: %default)");
	parser.add_option("-g") .type("int").dest("g").help("alternative help");
	parser.add_option("-b") .type("double").dest("b").help("Increment constant");
	parser.add_option("-r").type("int").dest("r").help("alternative version");
	parser.add_option("-t").type("int").dest("t").help("default:");

	Values& options = parser.parse_args(argc, argv);

	if(!options.is_set("n") || !options.is_set("i")){
		std::cout<<parser.usage()<<std::endl;
		exit(0);
	}

	MinLASolver solver(options);
	solver.start();
}
