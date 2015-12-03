#ifndef MINLASOLVER_H
#define MINLASOLVER_H

#include "graph.h"
#include <stdint.h>
#include "OptionParser.h"
#include "coarseedge.h"
class MinLASolver{
public:
	MinLASolver(optparse::Values& op);
	void start();
private:
	MinLASolver(const MinLASolver&);
	MinLASolver& operator=(const MinLASolver&);

	void ML(Graph& g);
	void uncoarsen(Graph& fineG, Graph& coarseG);
	void removeEdgesFromCoarse(Graph& fineG, Graph& coarseG);
	void markEdgesForDeletion(Graph& g);

	optparse::Values& options;
	std::vector<double> paramlist;
};
#endif
