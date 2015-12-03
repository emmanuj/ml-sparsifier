#include "graphsparsifier.h"
#include <algorithm>
#include <functional>
#include <cmath>
#include <fstream>
#include "etimer.h"
#include "util.h"
#include <sstream>
#include <stack>

/*
 	There are three steps to sparsification
	1. Weight the edges
	2. Filter the edges based on some heuristic
*/

Graph GraphSparsifier::sparsifyRandom(std::vector<Edge>& edges, long numnodes){
	Graph sGraph(numnodes);
	//Random edge selection
	int limit = ceil(0.2 * edges.size());
	int count =0;

	while(count != limit && edges.size() > 0){
		int idx = rand() % edges.size();

		Index first = edges[idx].first;
		Index second = edges[idx].second;
		edgeweight w = edges[idx].weight;
		sGraph.addEdge(first,second, w);

		edges.erase(edges.begin() + idx); //todo see whether to use a list here instead
		count++;
	}

	return sGraph;
}
inline double mround( double val )
{
    if( val < 0 ) return ceil(val - 0.5);
    return floor(val + 0.5);
}
inline long sparsifyGraph(Graph& orig, float e,std::ostream& fout, bool write, optparse::Values& options){
	Graph g = orig;
	g.computeAlgebraicDistances(options);
	long count = 0;
	std::vector<std::stack<long> > deleted(g.getSize());
	std::vector<double> v;
	double stdev =0;
	if(!(bool)options.get("weak") && ! (bool)options.get("strong")){
		for(int i=0;i<g.getSize();i++){
			std::vector<Neighbor> nbs = g.neighbors(i);
			for(int j=0;j<nbs.size();j++){
				//std::cout<<nbs[j].algebraicDist<<std::endl;
				v.push_back(nbs[j].algebraicDist);
			}
		}
		double sum = std::accumulate(v.begin(), v.end(), 0.0);
		double mean = sum / v.size();

		std::vector<double> diff(v.size());
		std::transform(v.begin(), v.end(), diff.begin(),
		               std::bind2nd(std::minus<double>(), mean));
		double sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
		stdev = std::sqrt(sq_sum / v.size());
	}
	

	for(int i=0;i<g.getSize();i++){
		const Index id = g.getNode(i).getNodeId();
		//delete previously stored edges
		bool rem = false;
		if(deleted[id].size() > 0){
			rem = true;
		}
		while(deleted[id].size() > 0 ){
			g.markEdgeForDeletion(id, g.getNode(id).getNeighborAt(deleted[id].top()));
			deleted[id].pop();
		}
		if(rem){
			g.removeMarkedEdges(id);
		}

		std::vector<Neighbor> nbs = g.neighbors(id);
		if(nbs.size() == 0) continue;
		//return top deg(u)^e
		int nh = ceil(pow(g.degree(id), e));
		if((bool)options.get("weak")){
			partial_sort(nbs.begin(), nbs.begin()+nh, nbs.end(), [=](Neighbor a, Neighbor b){
	            return a.algebraicDist < b.algebraicDist;
	        });
		}else if((bool) options.get("strong")){
			partial_sort(nbs.begin(), nbs.begin()+nh, nbs.end(), [=](Neighbor a, Neighbor b){
	            return a.algebraicDist > b.algebraicDist;
	        });
		}else{
			//divide into bins and select a few

			//sort and group the data
			std::sort(nbs.begin(), nbs.end(), [=](Neighbor a, Neighbor b){
	            return a.algebraicDist < b.algebraicDist;
	        });

			float h = (3.5 * stdev)/cbrt(nbs.size()); //bin width
			int k = ceil((nbs[nbs.size()-1].algebraicDist - nbs[0].algebraicDist)/h) + 1; //number of bins
			//std::cout<<"Node "<<id<<std::endl;
			//std::cout<<"#bins: "<<k<<" "<<" bin width: "<<h<<" # to remove :"<<nh<<" #data points "<<nbs.size()<<std::endl;

			std::vector<std::vector<Neighbor>> bins(k);

			int l =0;
			int idx = 0;
			float min_value = 0;
			float max_value = h;


			while( l < k && idx < nbs.size()){
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
			int x = nh;
			for(int p=0;p<bins.size();p++){
				int nh_bin =round((bins[p].size()*1.0)/(nbs.size()) * nh);
				int c = 0;
				int end = bins[p].size();
				while(c < nh_bin && x>0){
					int e = rand() % end;
					count++;
					if(write) fout<<id<<" "<< bins[p][e].nodeId<<std::endl;//" "<<ed.weight<<std::endl;
					g.markEdgeForDeletion(id,  bins[p][e]);
					deleted[bins[p][e].nodeId].push(bins[p][e].invPos);//stores the position of the second edge that needs to be removed
					//move the neigbor to back and reduce range for random
					Neighbor n = bins[p][e];
					bins[p][e] = bins[p][bins[p].size() - 1];
					bins[p][bins[p].size() - 1] = n;
					end--;

					c++;
					x--;
				}
				if(nh !=0) g.removeMarkedEdges(id);
			}

		}

		if((bool)options.get("weak") || (bool)options.get("strong")){
			int c = 0;
			while(c < nh){
				count++;
				if(write) fout<<id<<" "<< nbs[c].nodeId<<std::endl;//" "<<ed.weight<<std::endl;
				g.markEdgeForDeletion(id,  nbs[c]);
				deleted[nbs[c].nodeId].push(nbs[c].invPos);//stores the position of the second edge that needs to be removed
				c++;
			}
			if(nh !=0) g.removeMarkedEdges(id);
		}

	}
	return count;
}
//find the best paramter with binary search
/*inline float parameterize( Graph& g, float edgeRatio, float l, float u, int max_steps, bool findLowest, float& minEdgeRatio){
	std::cout<<"edge ratio: "<<edgeRatio<<std::endl;
	float upper_bound = u;
	float lower_bound = l;
	float bestParameter = lower_bound;
	float estimate = lower_bound;
	float minDistance = upper_bound;

	const float _ABS_ZERO = 1e-3;
	int step = 0;
	float prevEdgeRatio = 0;
	while(step < max_steps){
		std::cout<<"l: "<<lower_bound<<" u: "<<upper_bound<<std::endl;
		estimate = (lower_bound+upper_bound) /2.0;
		float sparseEdgeCount = (float) sparsifyGraph(g, estimate, std::cout, false);
		//the estimate has stoped improving the result.
		//In other words, we've achieved the minimum sparsification possible
		float currentEdgeRatio = sparseEdgeCount/g.getEdgeCount();
		
		currentEdgeRatio = mround(currentEdgeRatio * 1000.0) /1000.0;
		std::cout<<sparseEdgeCount<<" "<<currentEdgeRatio<<std::endl;

		if(fabs(currentEdgeRatio - prevEdgeRatio) < _ABS_ZERO && findLowest){
			minEdgeRatio = currentEdgeRatio;
			return estimate;
		}


		float distance = fabs(currentEdgeRatio - edgeRatio);
		if(distance < minDistance && fabs(currentEdgeRatio) > _ABS_ZERO){
			minDistance = distance;
			bestParameter = estimate;

			//"Exact" hit?
			if (fabs(currentEdgeRatio - edgeRatio) < _ABS_ZERO){
				std::cout<<"Break: "<<fabs(currentEdgeRatio - edgeRatio)<<std::endl;
				break;
			}
		}

		bool increase = (currentEdgeRatio < edgeRatio);

		if(increase){
			lower_bound = estimate;
		}
		else{
			upper_bound = estimate;
		}

		prevEdgeRatio = currentEdgeRatio;
		step++;
	}

	return bestParameter;
}
inline float parameterize2( Graph& g, float edgeRatio, float l, float u){
	float upper_bound = u;
	float lower_bound = l;
	float estimate = lower_bound;

	float currentEdgeRatio = 0.0;
	
	const float _ABS_ZERO = 1e-3;
	while(fabs(currentEdgeRatio - edgeRatio) > _ABS_ZERO){
		std::cout<<"Got here"<<std::endl;
		estimate = (lower_bound+upper_bound) * 0.5;
		std::cout<<"l: "<<lower_bound<<" u: "<<upper_bound<<std::endl;
		float sparseEdgeCount = (float) sparsifyGraph(g, estimate, std::cout, false);
		std::cout<<"Got here1"<<std::endl;
		currentEdgeRatio = sparseEdgeCount/g.getEdgeCount();
		std::cout<<"Got here2"<<std::endl;
		currentEdgeRatio = mround(currentEdgeRatio * 1000.0) /1000.0;
		std::cout<<"Got here3 "<<std::endl;
		std::cout<<sparseEdgeCount<<" "<<currentEdgeRatio<<std::endl;

		bool increase = (currentEdgeRatio < edgeRatio);
		if(increase){
			lower_bound = estimate;
		}else{
			upper_bound = estimate;
		}
	}

	return estimate;
}*/
void GraphSparsifier::sparsifyByBinning2(Graph& orig, optparse::Values& options){
	std::cout<<"Sparisifying graph... with edge count "<<orig.getEdgeCount()<<std::endl;
	std::string name = (std::string) options.get("gname");
	float frac = (float) options.get("s");
	long numrun = (long) options.get("b");

	float r = 0.1;	
	float lower_bound =0;
	const float upper_bound = 1;
	ETimer timer;
	//float e = parameterize(orig, r, lower_bound, upper_bound, 20, true, r);
	while(r < 1.0){
		std::ostringstream ss;
		ss << r;
		std::string perc = ss.str();
		for(int m=0;m<numrun;m++){
			std::vector<std::stack	<long> > deleted(orig.getSize());
			std::ofstream fout(("out/"+name+"/"+name+"_"+perc+"_"+std::to_string(m)+".txt").c_str());
			sparsifyGraph(orig, r, fout, true, options);
			fout.flush();
			fout.close();
		}
		r+=0.1;
	}
	std::cout<<"Done sparsifying all graphs"<<std::endl;
}
