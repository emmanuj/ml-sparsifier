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
#include "etimer.h"
#include "log.h"
#include <omp.h>

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
options(op), paramlist(), param(-1.0){
}

void printMarkedStats(Graph& g){
	long c = 0;
	long t =0;
	for(int i=0;i<g.getSize();i++){
		std::vector<Neighbor> nbs = g.neighbors(i);
		t+=nbs.size();
		for(Neighbor& n: nbs){
			if(n.isdelete){
				if(i<n.nodeId){
					c++;
				}
			}
		}
	}
	std::cout<<"Marked: "<<c<<", Total: "<<t/2<<std::endl;
}

void MinLASolver::markEdges(Graph& g){
	int count=0;
	double param = paramlist[g.getLevel()];
	std::cout<<"Level: "<<g.getLevel()<<" Params: "<<param<<std::endl;
	if(param > -1){
		std::vector<std::stack<long> > marked(g.getSize());
		for(unsigned i=0;i<g.getSize();i++){

			//mark previously stored edges
			while(marked[i].size() > 0 ){
				g.markEdge(i, g.getNode(i).getNeighborAt(marked[i].top()));
				marked[i].pop();
			}
			std::vector<Neighbor> nbs = g.getNode(i).neighbors();

			if(g.degree(i) == 0) continue;

			int nh = ceil(pow(g.degree(i), param));

			if(nh == 0) continue;
			if((bool)options.get("weak") || (bool)options.get("strong")){
				partial_sort(nbs.begin(),nbs.begin()+nh, nbs.end(), SortNeighbors(options));
				//mark. remember to take into account all the time edges already marked.
				int c = 0;
				while(c < nh && !nbs[c].isdelete){
					g.markEdge(i,  nbs[c]);
					marked[nbs[c].nodeId].push(nbs[c].invPos);//stores the position of the second edge that needs to be removed
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
				////std::cout<<"Node "<<id<<std::endl;
				////std::cout<<"#bins: "<<k<<" "<<" bin width: "<<h<<" # to remove :"<<nh<<" #data points "<<g.degree(i)<<std::endl;

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
							g.markEdge(i,  bins[p][pos]);
							marked[bins[p][pos].nodeId].push(bins[p][pos].invPos);//stores the position of the second edge that needs to be removed
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
						g.markEdge(i,  nbs[c]);
						marked[nbs[c].nodeId].push(nbs[c].invPos);//stores the position of the second edge that needs to be removed
						c++;
					}
				}

			}

		}
	}
}

void MinLASolver::inheritMarkedEdges(Graph& fineG, Graph& coarseG){
	////std::cout<<"Coarse edges: "<<std::endl;
	int count =0;
	for(unsigned i=0;i<coarseG.getSize();i++){
		std::vector<Neighbor>& nbs = coarseG.neighbors(i);
		for(Neighbor& n: nbs){
			if(n.isdelete){
				if(i < n.nodeId) continue;
				////std::cout<<i<<" "<<n.nodeId<<" "<<n.fineEdges.size()<<std::endl;
				for(std::pair<Index, Index>& f_edge: n.fineEdges){
					Neighbor nb = fineG.getNode(f_edge.first).getNeighborAt(f_edge.second);
					////std::cout<<"\t"<<f_edge.first<<" "<<nb.nodeId<<std::endl;
					fineG.markEdge(f_edge.first, fineG.getNode(f_edge.first).getNeighborAt(f_edge.second));
					fineG.markEdge(nb.nodeId, fineG.getNode(nb.nodeId).getNeighborAt(nb.invPos));
					count++;
				}
			}
		}
	}
}

void MinLASolver::uncoarsen(Graph& fineG, Graph& coarseG){
	std::cout<<"Uncoarsening... Level: "<<fineG.getLevel()<<std::endl;
	inheritMarkedEdges(fineG, coarseG);
	fineG.computeAlgebraicDistances(options);
	markEdges(fineG);
	printMarkedStats(fineG);
}

void MinLASolver::ML(Graph& g){
	int l = g.getLevel();
	std::cout<<"==================== level: "<<l<<" ======================"<<std::endl;
	//std::cout<<"********************Graph stats ****************************"<<std::endl;
	//std::cout<<" # nodes: "<<g.getNodes().size()<<std::endl;
	//std::cout<<" #Edges "<< g.getEdgeCount()<<std::endl;

	//long max = 0;
	/*long min = 9999999999;
	for(int i=0;i<g.getSize();i++){
		if(max < g.degree(i)){
			max = g.degree(i);
		}
		if(g.degree(i)< min){
			min= g.degree(i);
		}
	}

	//std::cout<<" Max Degree: "<<max<<std::endl;
	//std::cout<<" Min Degree: "<<min<<std::endl;
	//std::cout<<"********************end graph stats ************************"<<std::endl;*/

	if(g.getSize()==10){
		//g.print();
		std::cout<<"In coarsest level. computing solution..."<<std::endl;
		ETimer timer;
		timer.start();

		g.computeAlgebraicDistances(options);

		paramlist = std::vector<double>(g.getLevel()+1, -1.0);

		std::string plist = (std::string) options.get("param-list");
		std::vector<std::string> plist_st = Util::split(plist,',');

		int k = paramlist.size();
		for(int i=0;i< plist_st.size() && k>0;i++){
			paramlist[k-1] = std::stod(plist_st[i]);
			k--;
		}

		//Benchmark code
		//create parameter list for sparsifying for each level
		/*int levelSpan = (int) options.get("level-span");
		int part = (int) options.get("sparse-level");
		//double param = (double) options.get("s");
		//std::cout<<"here "<<std::endl;
		//paramlist = std::vector<double>(g.getLevel()+1,param);
		paramlist = std::vector<double>(g.getLevel()+1,-1.0);
		if(part ==0){ //sparsify only at coarsest levels
			//double p =0;
			int startIdx = paramlist.size() -1;
			for(int i= startIdx; i > (startIdx - levelSpan) && i > 0; i--){
				paramlist[i] = param;
				//p+=0.1;
			}
		}else if(part == 1){// sparsify at the mid levels
			int startIdx = floor(paramlist.size()/2) - ceil(levelSpan/2);
			//double p =param + 0.3;
			for(int i=startIdx;i < (startIdx+levelSpan) && i < paramlist.size(); i++){
				paramlist[i] = param;
				//p-=0.1;
			}
		}else if(part == 2){ //sparsify at the finest levels
			int startIdx = 0;
			//double p =param + 0.3;
			for(int i= startIdx; i < (startIdx + levelSpan) && i < paramlist.size(); i++){
				paramlist[i] = param;
				//p-=0.1;
			}
		}
		//std::cout<<"Got here "<<std::endl;*/
		markEdges(g);
		//printMarkedStats(g);
		timer.stop("Completed in ");
	}else{
		Graph coarseG = GraphUtil::coarsenGraph(g, options);
		for(int i=0;i<(int)options.get("g");i++){
			ML(coarseG);
			uncoarsen(g, coarseG);
		}
	}
}
inline double mround( double val )
{
    if( val < 0 ) return ceil(val - 0.5);
    return floor(val + 0.5);
}
//parameter search for benchmarking
void MinLASolver::binarySearchParameter(Graph& g, double target_lower, double target_upper){

	double lo =0.0;
	double hi = 1.0;
	Graph sparse;
	int i =0;
	int niter = 20;
	while(i < niter){
		double mid = lo + ((hi-lo)*0.5);
		param = mid;
		sparse =g;
		ML(sparse);

		double actual_ratio = mround((sparse.getEdgeCount()/(g.getEdgeCount() * 1.0)) * 100 ) /100.0;
		std::cout<<"Run: "<<i<<"Orig: "<<g.getEdgeCount()<<" Sparse count: "<<sparse.getEdgeCount()<<" actual ratio: "<<actual_ratio<<std::endl;
		if(actual_ratio >= target_lower && actual_ratio <= target_upper){
			break;
		}else if(actual_ratio > target_upper){
			lo = mid;
		}else{//} if(actual_ratio < target_lower){
			hi = mid;
		}
	}

	g = sparse;
}
void MinLASolver::start(){
	ETimer timer;
	timer.start();
	//read in graph
	std::cout<<"Reading graph.... "<<std::endl;
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

	timer.stop("Finished reading graph in");

	std::vector<Node> solution;
	//std::cout<<"Started coarsening..."<<std::endl;
	//timer.start();

	int edge0 = graph.getEdgeCount();
	timer.start();

	std::cout<<"Sparsifying graph..."<<std::endl;
	if((bool)options.get("single-level")){

		paramlist = std::vector<double>(1);
		paramlist[0] = std::stod((std::string)options.get("s"));
		markEdges(graph);
	}else{
		ML(graph);
	}

	timer.stop("Algorithm completed in");

	//lets try to find the parameter
	//binarySearchParameter(graph, 0.2, 0.4);
	//timer.stop((std::string)options.get("gname"));
	

	//Lets print the sparsified graph
	std::cout<<"Writing sparsified graph"<<std::endl;
	std::ofstream fout((std::string)options.get("gname")+"_"+(std::string)options.get("o"));
	std::vector<long> degs = graph.getDegrees();
	int edge1 =0;
	for(int i=0;i<graph.getSize();i++){
		const std::vector<Neighbor> nbs = graph.neighbors(i);
		//ensure at least one edge per node
		if(!((bool)options.get("single-level"))){
			if(graph.degree(i) == 0 && degs[i] == 0 ){
				Neighbor n = nbs[rand() % nbs.size() ];
				fout<<i<<" "<< n.nodeId<<std::endl;
				degs[i]++;
				edge1++;
				continue;
			}
		}
		for(const Neighbor& n: nbs){
			if(i<n.nodeId) continue;
			if(n.isdelete){
				edge1++;
				fout<<i<<" "<< n.nodeId<<std::endl;
			}
		}
	}

	fout.flush();
	fout.close();


	
	std::cout<<"Orig # edges: "<<edge0<<" new edge count: "<<edge1<<std::endl;
	//printMarkedStats(graph);
}
