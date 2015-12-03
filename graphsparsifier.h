#ifndef GRAPHSPARSIFY_H
#define GRAPHSPARSIFY_H
#include <iostream>
#include <vector>
#include "edge.h"
#include "neighbor.h"
#include "node.h"
#include "graph.h"
#include "OptionParser.h"

class GraphSparsifier{
	public:
		static Graph sparsifyRandom(std::vector<Edge>& edges, long numnodes);
		static void sparsifyByStrata(Graph& g, int edgeCount, int numClus, float frac, std::string& infile);
		static void writeClusterToFile(Graph& g, int numClus, std::string& outfile);
		static Graph sparsifyBySpectral(std::vector<Edge>& edges, long numnodes);
		static void sparsifyByBinning(Graph& g, float frac, long bin_min);
		static void sparsifyByBinning(Graph& g, float frac);
		static void sparsifyByBinning2(Graph& g, optparse::Values& options);
};

#endif
