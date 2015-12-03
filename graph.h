#ifndef GRAPH_H
#define GRAPH_H

#include <vector>
#include "globals.h"
#include "node.h"
#include <map>
#include <utility>
#include <vector>
#include <list>
#include "OptionParser.h"
#include "coarseedge.h"
class Graph{
private:
	std::vector<Node> nodes;
	count numNodes;
	int level;
	long edge_count;
	std::vector<long> degs; //node degrees
public:
	Graph();
	Graph(const Graph& g);
	Graph(const count& n); //default
	//Graph(const count& n, std::vector<Node>& vertices,  std::vector<std::pair<Index, Index> >& edges, std::vector<double>& weights); //default
	~Graph(); //destructor
	Graph& operator=(const Graph& g){
		if(this != &g){
			nodes=g.nodes;
			numNodes = g.numNodes;
			level = g.level;
			edge_count =g.edge_count;
			degs = g.degs;
		}
		return *this;
	}
	Graph(std::vector<Index>&cnodes, std::vector<std::vector< CoarseEdge >>&cedges,
		Graph& fineGraph,std::vector<Index>& mapper, int level);

	void print();
	void printGraph(const std::string& fname);
	Graph sparsifyGraph();
	void printStats();
	std::vector<Neighbor>& neighbors(const Index idx);
	std::vector<Neighbor>& coarseNeighbors(const Index idx);
	void addEdge(const Index &node1, const Index &node2, const edgeweight& weight);
	void addEdge(const Index &node1, const Index &node2, const edgeweight& weight, const std::vector<std::pair<Index, Index> >& fineEdges );
	void addEdges(std::map<Index, std::map<Index, edgeweight> >& e );
	edgeweight weight(Index idx,Index idx2) const;
	edgeweight weightedDegree(const Index idx);
	count degree(const Index& idx);
	const std::vector<Node>& getNodes()const {return nodes;}
	
	long getEdgeCount() const{ return edge_count/2;}
	Node& getNode(const Index& i ){ return nodes[i]; }
	void setCoarse(const Index& idx, const bool isCoarse){ //sets the coarse flag for a node
		nodes[idx].setCoarse(isCoarse);
		//update coarse neighborhood of neighbors
		std::vector<Neighbor> nbs = nodes[idx].neighbors();
		for(unsigned i=0;i<nbs.size();i++){
			Neighbor nb(idx, nbs[i].weight);
			nb.algebraicDist = nbs[i].algebraicDist;
			nodes[nbs[i].nodeId].addCNeighbor(nb);
		}
	}
	void setFutureVolume(const Index& idx, const double& volume){ nodes[idx].setFutureVolume(volume);}
	void setVolume(const Index& idx, const long double& volume){ nodes[idx].setVolume(volume);}

	/**
	* @return count - number of nodes in the graph
	*/
	count getSize() const { return numNodes; }
	int getLevel()const{ return level; }
	void setLevel(const int _level){ level =_level; }
	void computeAlgebraicDistances(optparse::Values& options);
	void printUnweightedGraph(std::ofstream& outfile);
	void removeMarkedEdges(Index from);
	void markEdgeForDeletion(Index u, Neighbor& v);


};
#endif
