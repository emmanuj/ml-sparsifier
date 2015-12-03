#include <algorithm>
#include "time.h"
#include <string>
#include <fstream>
#include <vector>
#include <cmath>
#include "node.h"
#include "graph.h"
#include "globals.h"
#include "solver.h"
#include <map>
#include <stack>
#include "util.h"
#include "edge.h"
#include <iomanip>
#include "OptionParser.h"
#include "graphutil.h"
#include "reader.h"
#include "edgelistreader.h"
#include "graphreader.h"
#include "graphsparsifier.h"
#include "etimer.h"
#include "log.h"

class SortNeighbors{
public:
	SortNeighbors(optparse::Values& op): options(op){}
	bool operator() (const Neighbor& lhs, const Neighbor& rhs){
       	if(!lhs.isdelete){
       		if(rhs.isdelete){
        		return true;
        	}else{
        		if((bool)options.get("weak")){
        			return lhs.algebraicDist < rhs.algebraicDist;
        		}else{
        			return lhs.algebraicDist > rhs.algebraicDist;
        		}

        	}
        }
        return false;
    }
private:
	optparse::Values& options;
};

MinLASolver::MinLASolver(optparse::Values& op):
options(op), paramlist(){
}

void printMarkedStats(Graph& graph){
	long c = 0;
	long t =0;
	for(int i=0;i<graph.getSize();i++){
		const std::vector<Neighbor> nbs = graph.neighbors(i);
		t+= nbs.size();
		for(const Neighbor& n: nbs){
			if(n.isdelete){
				if(i<n.nodeId) continue;
				c++;
			}
		}
	}
	std::cout<<"Marked: "<<c<<", Total: "<<t/2<<std::endl;
}
void MinLASolver::markEdgesForDeletion(Graph& g){
	int count=0;
	double param = paramlist[g.getLevel()];
	std::cout<<"Level: "<<g.getLevel()<<" Params: "<<param<<std::endl;
	if(param > -1){
		std::vector<std::stack<long> > deleted(g.getSize());
		for(unsigned i=0;i<g.getSize();i++){

			//mark previously stored edges
			while(deleted[i].size() > 0 ){
				g.markEdgeForDeletion(i, g.getNode(i).getNeighborAt(deleted[i].top()));
				deleted[i].pop();
			}
			std::vector<Neighbor> nbs = g.getNode(i).neighbors();

			if(g.degree(i) == 0) continue;

			int nh = round(pow(g.degree(i), param)); //TODO: Can change this to either a floor, a ceil and a round. Explore the different options.

			if(nh == 0) continue;
			if((bool)options.get("weak") || (bool)options.get("strong")){
				partial_sort(nbs.begin(),nbs.begin()+nh, nbs.end(), SortNeighbors(options));
				//mark. remember to take into account all the time edges already marked.
				int c = 0;
				while(c < nh && !nbs[c].isdelete){
					g.markEdgeForDeletion(i,  nbs[c]);
					deleted[nbs[c].nodeId].push(nbs[c].invPos);//stores the position of the second edge that needs to be removed
					c++;
				}
			}else{ //a little bit of both
				//divide into bins and select a few

				//sort and group the data
				std::sort(nbs.begin(), nbs.end(), [=](Neighbor a, Neighbor b){
		            if(!a.isdelete){
			       		if(b.isdelete){
			        		return true;
			        	}else{
			        		return a.algebraicDist < b.algebraicDist;
			        	}
			        }
			        return false;
		        });
		        double stdev =0;
				std::vector<double> v(g.degree(i));
				for(int j=0;j<g.degree(i);j++){
					v[j] = nbs[j].algebraicDist;
				}

				double sum = std::accumulate(v.begin(), v.end(), 0.0);
				double mean = sum / v.size();

				std::vector<double> diff(v.size());
				std::transform(v.begin(), v.end(), diff.begin(),
				               std::bind2nd(std::minus<double>(), mean));
				double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
				stdev = std::sqrt(sq_sum / v.size());
				//clean up
				std::vector<double>().swap(v);
				std::vector<double>().swap(diff);

				float h = (3.5 * stdev)/cbrt(g.degree(i)); //bin width
				int k = h==0 ? 1 :ceil((nbs[g.degree(i)-1].algebraicDist - nbs[0].algebraicDist)/h) + 1; //number of bins
				//std::cout<<"Node "<<id<<std::endl;
				//std::cout<<"#bins: "<<k<<" "<<" bin width: "<<h<<" # to remove :"<<nh<<" #data points "<<g.degree(i)<<std::endl;

				int l =0;
				int idx = 0;
				double min_value = nbs[0].algebraicDist;
				double max_value =min_value + h;

				if(k!=1){
					std::vector<std::vector<Neighbor>> bins(k);
					while( l < k && idx < nbs.size() && !nbs[idx].isdelete){
						if(nbs[idx].algebraicDist >= min_value && nbs[idx].algebraicDist < max_value){
							bins[l].push_back(nbs[idx]);
							idx++;
						}

						if(nbs[idx].algebraicDist >= max_value){
							l++;
							min_value = max_value;
							max_value = min_value + h;
						}
					}

					//now filter by bins
					long deg = g.degree(i);
					std::random_shuffle ( bins.begin(), bins.end() );//random shuffle the bins
					for(int p=0;p<bins.size();p++){

						if(bins[p].size() == 0) continue;
						int nh_bin =ceil((bins[p].size()*1.0 * nh)/deg);
						int c = 0;
						int end = bins[p].size();
						while(c < nh_bin && end > 0){
							int pos = rand() % end;
							g.markEdgeForDeletion(i,  bins[p][pos]);
							deleted[bins[p][pos].nodeId].push(bins[p][pos].invPos);//stores the position of the second edge that needs to be removed
							//move the neigbor to back and reduce range for random
							if(pos != (end - 1)){
								Neighbor n = bins[p][pos];
								bins[p][pos] = bins[p][end - 1];
								bins[p][end - 1] = n;
							}
							end--;
							c++;
						}
					}

				}else{
					int c = 0;
					while(c < nh && !nbs[c].isdelete){
						g.markEdgeForDeletion(i,  nbs[c]);
						deleted[nbs[c].nodeId].push(nbs[c].invPos);//stores the position of the second edge that needs to be removed
						c++;
					}
				}

			}

		}
	}
}

void MinLASolver::removeEdgesFromCoarse(Graph& fineG, Graph& coarseG){
	//std::cout<<"Coarse edges: "<<std::endl;
	int count =0;
	for(unsigned i=0;i<coarseG.getSize();i++){
		std::vector<Neighbor>& nbs = coarseG.neighbors(i);
		for(Neighbor& n: nbs){
			if(n.isdelete){
				if(i < n.nodeId) continue;
				//std::cout<<i<<" "<<n.nodeId<<" "<<n.fineEdges.size()<<std::endl;
				for(std::pair<Index, Index>& f_edge: n.fineEdges){
					Neighbor nb = fineG.getNode(f_edge.first).getNeighborAt(f_edge.second);
					//std::cout<<"\t"<<f_edge.first<<" "<<nb.nodeId<<std::endl;
					fineG.markEdgeForDeletion(f_edge.first, fineG.getNode(f_edge.first).getNeighborAt(f_edge.second));
					fineG.markEdgeForDeletion(nb.nodeId, fineG.getNode(nb.nodeId).getNeighborAt(nb.invPos));
					count++;
				}
			}
		}
	}
}

void MinLASolver::uncoarsen(Graph& fineG, Graph& coarseG){
	std::cout<<"Uncoarsening... Level: "<<fineG.getLevel()<<std::endl;
	removeEdgesFromCoarse(fineG, coarseG);
	fineG.computeAlgebraicDistances(options);
	markEdgesForDeletion(fineG);
	//printMarkedStats(fineG);
}

void MinLASolver::ML(Graph& g){
	int l = g.getLevel();
	std::cout<<"==================== level: "<<l<<" ======================"<<std::endl;
	std::cout<<"********************Graph stats ****************************"<<std::endl;
	std::cout<<" # nodes: "<<g.getNodes().size()<<std::endl;
	std::cout<<" #Edges "<< g.getEdgeCount()<<std::endl;

	long max = 0;
	long min = 9999999999;
	for(int i=0;i<g.getSize();i++){
		if(max < g.degree(i)){
			max = g.degree(i);
		}
		if(g.degree(i)< min){
			min= g.degree(i);
		}
	}

	std::cout<<" Max Degree: "<<max<<std::endl;
	std::cout<<" Min Degree: "<<min<<std::endl;
	std::cout<<"********************end graph stats ************************"<<std::endl;

	if(g.getSize()==10){
		//g.print();
		std::cout<<"In coarsest level. computing solution..."<<std::endl;
		ETimer timer;
		timer.start();
		g.computeAlgebraicDistances(options);

		//create parameter list for sparsifying for each level
		int levelSpan = (int) options.get("level-span");
		int part = (int) options.get("sparse-level");
		double param = (double) options.get("s");

		paramlist = std::vector<double>(g.getLevel()+1,-1.0);
		if(part ==0){ //sparsify only at coarsest levels
			int startIdx = paramlist.size() -1;
			for(int i= startIdx; i > (startIdx - levelSpan) && i > 0; i--){
				paramlist[i] = param;
			}
		}else if(part == 1){// sparsify at the mid levels
			int startIdx = floor(paramlist.size()/2) - ceil(levelSpan/2);
			for(int i=startIdx;i < (startIdx+levelSpan) && i < paramlist.size(); i++){
				paramlist[i] = param;
			}
		}else if(part == 2){ //sparsify at the finest levels
			int startIdx = 0;
			for(int i= startIdx; i < (startIdx + levelSpan) && i < paramlist.size(); i++){
				paramlist[i] = param;
			}
		}
		std::cout<<"Got here "<<std::endl;
		markEdgesForDeletion(g);
		printMarkedStats(g);
		timer.stop("Completed in ");
	}else{
		Graph coarseG = GraphUtil::coarsenGraph(g, options);
		for(int i=0;i<(int)options.get("g");i++){
			ML(coarseG);
			uncoarsen(g, coarseG);
		}
	}
}

void MinLASolver::start(){
	ETimer timer;
	timer.start();
	//read in graph
	//std::cout<<"Reading graph.... "<<std::endl;
	//TODO: check if file exist
	//TODO: do some error checking on the graph
	std::ifstream listfile((std::string)options.get("i"));

	EdgeListReader elist(options);
	GraphReader gReader(options);

	Graph graph;

	if((std::string)options.get("f") == "graph"){
		graph = gReader.read(listfile, (long)options.get("n"));
	}else if((std::string)options.get("f") == "edgelist"){
		graph = elist.read(listfile, (long) options.get("n"));
	}

	timer.stop("Finished reading graph in ");

	std::vector<Node> solution;
	std::cout<<"Started coarsening..."<<std::endl;
	timer.start();

	int edge0 = graph.getEdgeCount();

	ML(graph);
	timer.stop("Algorithm completed in");

	//Lets print the sparsified graph
	std::ofstream fout((std::string)options.get("o"));
	for(int i=0;i<graph.getSize();i++){
		const std::vector<Neighbor> nbs = graph.neighbors(i);
		for(const Neighbor& n: nbs){
			if(i<n.nodeId) continue;
			if((!n.isdelete)||(graph.degree(n.nodeId) <= 1 || graph.degree(i) <=1)){
				fout<<i<<" "<< n.nodeId<<std::endl;
			}
		}
	}

	int edge1 = graph.getEdgeCount();
	std::cout<<"Orig # edges: "<<edge0<<" new # edge: "<<edge1<<std::endl;
	//printMarkedStats(graph);
}
