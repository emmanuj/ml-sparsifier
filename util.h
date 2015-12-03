#ifndef UTIL__H
#define UTIL__H

#include <vector>
#include <string>
#include "node.h"
#include "graph.h"

class SortByVPrime{
public:
	SortByVPrime(Graph& graph):g(graph){}
	bool operator()(const Node& lhs, const Node& rhs){
		int lprime = 0;
		int rprime = 0;

		std::vector<Neighbor> lnbs = g.neighbors(lhs.getNodeId());
		std::vector<Neighbor> rnbs = g.neighbors(rhs.getNodeId());
		for(int i=0;i< lnbs.size();i++){
			if(g.getNode(lnbs[i].nodeId).inVPrime()){
				++lprime;
			}
		}
		for(int i=0;i< rnbs.size();i++){
			if(g.getNode(rnbs[i].nodeId).inVPrime()){
				++rprime;
			}
		}
		return lprime > rprime;
	}
private:
	Graph& g;
};

class SortByYi
{
	public:
	    bool operator() (const Node& lhs, const Node& rhs){
	       return lhs.getY() < rhs.getY();
	    }
};


class Util{
public:
    static std::vector<std::string> splitv(const std::string &s, char delim, std::vector<std::string> &elems);
    static std::vector<std::string> split(const std::string &s, char delim);
    static std::string trim(std::string &s);
    static std::vector<double> computeXi(std::vector<Node> &nodes);
    static void printCoarseGraph(Graph& g);
    static void doPumping(Graph& g, std::vector<Node>& nodes);
private:

};
#endif
