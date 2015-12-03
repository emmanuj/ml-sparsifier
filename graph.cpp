#include "graph.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include <iomanip>
#include "util.h"
#include "edge.h"
#include <fstream>
#include <cmath>


class SortByAlgebraicDistance
{
	public:
	    bool operator() (const Edge& lhs, const Edge& rhs){
	       return lhs.weight < rhs.weight;
	    }
};

Graph::~Graph(){}

Graph::Graph():
nodes(),
numNodes(),
level(0),
edge_count(0),
degs(){}

Graph::Graph(const Graph& g):
nodes(g.nodes),
numNodes(g.numNodes),
level(g.level),
edge_count(g.edge_count),
degs(g.degs)
{
}

Graph::Graph(const count& n):
nodes(n),
numNodes(n),
level(0),
edge_count(0),
degs(n)
{
	for(Index i=0; i<numNodes;i++){
		nodes[i] = Node(i);
	}
}

Graph::Graph(std::vector<Index> &nds, std::vector<std::vector<CoarseEdge> >& coarseEdges,
	Graph& fineGraph, std::vector<Index>& mapper, int l):
nodes(nds.size()),
numNodes(nds.size()),
level(l),
edge_count(0),
degs(nds.size())
{
	for(Index i=0;i<nds.size();i++){
		//copy node properties
		nodes[i].setNodeId(i);
		nodes[i].setPrevId(nds[i]);
	}
	//add edges
	for(Index i=0;i<coarseEdges.size();i++){
		if(coarseEdges[i].size() == 0) continue;
		//sort them by index
		sort(coarseEdges[i].begin(), coarseEdges[i].end(),[=](CoarseEdge a, CoarseEdge b){
            return a.nodeId > b.nodeId;
        });

		int j =0;
		Index cur_idx = coarseEdges[i][j].nodeId;
		edgeweight sum = 0;

		while(j < coarseEdges[i].size()){
			std::vector< std::pair<Index, Index> > fineEdges;
			while( j< coarseEdges[i].size() && cur_idx == coarseEdges[i].at(j).nodeId ){
				sum += coarseEdges[i][j].weight;
				fineEdges.push_back(coarseEdges[i][j].fineEdge);
				j++;
			}
			addEdge(i, cur_idx, sum, fineEdges);
			//std::cout<<"======= Edge ["<<nodes[i].getPrevId()<<" : "<<nodes[cur_idx].getPrevId()<<"] ======="<<std::endl;
			//for(std::pair<Index,Index>& f: fineEdges){
			//	std::cout<<f.first<<" - "<<f.second<<std::endl;
			//}
			sum=0;
			if(j != coarseEdges[i].size() ){
				cur_idx = coarseEdges[i].at(j).nodeId;
			}
		}
	}

	long double volSum =0.0;
	for(int i=0;i<nds.size();i++){
		setVolume(mapper[nds[i]], fineGraph.getNode(nds[i]).V());
	}
	//compute volumes
	for(int i=0;i<fineGraph.getSize();i++){
		if(fineGraph.getNode(i).isC()) continue;
		Index idx = fineGraph.getNode(i).getNodeId();
		std::vector<Neighbor>& cnbs = fineGraph.coarseNeighbors(idx);
		for(unsigned j=0;j<cnbs.size();j++){
			double p = cnbs.at(j).pij;
			double vol = (p * fineGraph.getNode(idx).V());
			setVolume(mapper[cnbs[j].nodeId], nodes[mapper[cnbs[j].nodeId]].V() + vol);
		}
	}
	//if(level ==2) print();
	for(Node& n: nodes){
		volSum+=n.V();
	}
	//if((int)ceil(volSum) != 100 ) print();
	std::cout<<"Total volume: "<<volSum<<std::endl;
	//print();
}
double getRandom(){
	return (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) - 0.5;
}

void Graph::computeAlgebraicDistances(optparse::Values& options){
	//srand(time(0));
    const double alfa = 0.5;
    const double epsilon = 0.001;
    const int Ri = 10;
    for(int h=0;h<Ri;h++){
        std::vector<double> x_vec(nodes.size());
        std::generate(x_vec.begin(), x_vec.end(), getRandom);//generate the random number between -0.5 to 0.5
        for(int k=0;k<40;k++){
            std::vector<double> c_vec(nodes.size());
            for(int i=0;i<nodes.size();i++){
                std::vector<Neighbor>& nbs = neighbors(nodes[i].getNodeId());
                double prod_sum =0;
                double sum = 0;
                for(Neighbor& n: nbs){
                    if(!n.isdelete){
                    	sum +=n.weight;
                    	prod_sum+= (n.weight * x_vec[n.nodeId]);
                    }
                }

                if(sum == 0){ //Avoid divide by zero error
                	c_vec[i] = (alfa * x_vec[i]);
                }else{
                	c_vec[i] = (alfa * x_vec[i]) + ( (1 - alfa) * (prod_sum/sum));
                }

            }
            x_vec =  c_vec;
        }

        //rescale
        double mn = *std::min_element(x_vec.begin(), x_vec.end());
        double mx = *std::max_element(x_vec.begin(), x_vec.end());
        for(int i=0;i<x_vec.size();i++){
            const double a = x_vec[i] - mn;
            const double b = mx - x_vec[i];
            double scaledX = (0.5 * ( a - b))/(a + b);
            x_vec[i] = scaledX;
        }
        for(int i=0;i<nodes.size();i++){
            for(int j=0;j<nodes[i].neighbors().size();j++){
                Index nid = nodes[i].neighbors()[j].nodeId;
                nodes[i].neighbors()[j].Ro.push_back(fabs(x_vec[i] - x_vec[nid]));
            }
        }
    }

    for(int i=0;i<nodes.size();i++){
        for(int j=0;j<nodes[i].neighbors().size();j++){
            double mx = *std::max_element(nodes[i].neighbors()[j].Ro.begin(), nodes[i].neighbors()[j].Ro.end());
            nodes[i].neighbors()[j].algebraicDist = (1/(mx + epsilon));
            if((bool)options.get("normalize")){
	            double ad = nodes[i].neighbors()[j].algebraicDist;
	            Index nbid = nodes[i].neighbors()[j].nodeId;
	            nodes[i].neighbors()[j].algebraicDist = ad /(sqrt(degree(nodes[i].getNodeId())* degree(nbid)));
	        }
	        std::vector<double>().swap(nodes[i].neighbors()[j].Ro);
        }
    }
}


void Graph::printGraph(const std::string& fname){
	std::ofstream outfile(std::string(fname+"_"+std::to_string(edge_count)).c_str());
	for(int i=0;i<nodes.size();i++){
		std::vector<Neighbor> nbs = neighbors(nodes[i].getNodeId());
		for(int j=0; j<nbs.size();j++){
			outfile<<nodes[i].getNodeId()<<" "<<nbs[j].nodeId<<" "<<nbs[j].weight<<std::endl;
		}
	}
}

void Graph::printUnweightedGraph(std::ofstream& outfile){
	for(int i=0;i<nodes.size();i++){
		std::vector<Neighbor> nbs = neighbors(nodes[i].getNodeId());
		for(int j=0; j<nbs.size();j++){
			if(nodes[i].getNodeId() > nbs[j].nodeId) continue;
			outfile<<nodes[i].getNodeId()<<" "<<nbs[j].nodeId<<std::endl;
		}
	}
}
void Graph::print(){
	std::cout<<"============== Level: "<<level<<" =============="<<std::endl;
	for(int i=0;i<nodes.size();i++){
		std::vector<Neighbor>& nbs = neighbors(nodes[i].getNodeId());
		if(nbs.size() == 0 ){
			std::cout<<nodes[i].getNodeId()<<" has no neighbors. Size: "<<getSize()<<std::endl;
			exit(0);
		}
		for(int j=0; j<nbs.size();j++){
			if(i > nbs[j].nodeId) std::cout<<nodes[i].getNodeId()<<" : "<<nbs[j].nodeId<<std::endl;
			//std::cout<<nodes[i].getNodeId()<<" : "<<nbs.at(j).nodeId<<" : "<<nbs.at(j).weight<<std::endl;
		}
	}
	
	std::cout<<"#edges: "<<getEdgeCount()<<std::endl;
}
void Graph::printStats(){
	int edgeCount =0;
	int minDegree = numNodes;
	int maxDegree =0;
	//double avgDegree=0.0;
	int sumDegree = 0;
	for(int i=0;i<nodes.size();i++){
		std::vector<Neighbor> nbs = nodes[i].neighbors();
		int deg = degree(nodes[i].getNodeId());
		if(deg > maxDegree){
			maxDegree = deg;
		}
		if(deg<minDegree){
			minDegree = deg;
		}
		sumDegree+=deg;
		for(int j=0; j<nbs.size();j++){
			edgeCount++;
		}
	}
	std::cout<<"[INPUT STATS: ]#nodes: "<<numNodes<<", #edges: "<<edgeCount<<", degree stat: ("<<"min: "<<minDegree<<" avg: "<<(sumDegree/(numNodes*1.0))<<" max: "<<maxDegree<<")"<<std::endl;
}
std::vector<Neighbor>& Graph::neighbors(const Index idx){
	return nodes.at(idx).neighbors();
}
std::vector<Neighbor>& Graph::coarseNeighbors(const Index idx){
	return nodes.at(idx).cNeighbors();
}

void Graph::addEdge(const Index& n1, const Index& n2, const edgeweight& weight){
	if(weight == 0){
		std::cout<<"Error: Edge weight can't be 0."<<std::endl;
		exit(1);
	}
	long pos1 = nodes.at(n1).addNeighbor(n2, weight);
	long pos2 = nodes.at(n2).addNeighbor(n1, weight);
	//std::cout<<pos1<<" "<<pos2<<std::endl;
	nodes.at(n1).getNeighborAt(pos1).invPos = pos2;
	nodes.at(n2).getNeighborAt(pos2).invPos = pos1;

	nodes.at(n1).setWeightedDegree(nodes.at(n1).getWeightedDegree() + weight );
	nodes.at(n2).setWeightedDegree(nodes.at(n2).getWeightedDegree() + weight );

	degs[n1] = degs[n1] + 1;
	degs[n2] = degs[n2] + 1;

	edge_count+=2;
}

void Graph::addEdge(const Index& n1, const Index& n2, const edgeweight& weight, const std::vector<std::pair<Index, Index>>& fineEdges){
	if(weight == 0){
		std::cout<<"Error: Edge weight can't be 0."<<std::endl;
		exit(1);
	}
	long pos1 = nodes.at(n1).addNeighbor(n2, weight, fineEdges);
	long pos2 = nodes.at(n2).addNeighbor(n1, weight, fineEdges);
	//std::cout<<pos1<<" "<<pos2<<std::endl;
	nodes.at(n1).getNeighborAt(pos1).invPos = pos2;
	nodes.at(n2).getNeighborAt(pos2).invPos = pos1;

	nodes.at(n1).getNeighborAt(pos1).fineEdges = fineEdges;
	nodes.at(n2).getNeighborAt(pos2).fineEdges = fineEdges;

	nodes.at(n1).setWeightedDegree(nodes.at(n1).getWeightedDegree() + weight );
	nodes.at(n2).setWeightedDegree(nodes.at(n2).getWeightedDegree() + weight );
	degs[n1] = degs[n1] + 1;
	degs[n2] = degs[n2] + 1;
	edge_count+=2;
}

void Graph::removeMarkedEdges(Index u){
	nodes.at(u).removeMarkedEdges();
}
void Graph::markEdgeForDeletion(Index u, Neighbor& v){
	if(degs[u] !=0 && !v.isdelete){
		nodes.at(u).neighbors()[v.pos].isdelete = true;
		nodes.at(u).setWeightedDegree(nodes.at(u).getWeightedDegree() - v.weight );
		edge_count--;
		degs.at(u) = degs.at(u) - 1;
	}
}

edgeweight Graph::weightedDegree(const Index idx){
	return nodes.at(idx).getWeightedDegree();
}

count Graph::degree(const Index& idx){
	return degs.at(idx);
}
