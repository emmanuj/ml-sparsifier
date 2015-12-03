#ifndef GRAPHREADER_H
#define GRAPHREADER_H

#include "reader.h"
#include "graph.h"
#include <iostream>
#include "util.h"
#include <fstream>
#include <string>
#include <algorithm>
class GraphReader: Reader
{
    public:
    	GraphReader(const optparse::Values& op):options(op), graph(0){}

        Graph& read(std::ifstream& in, long n){

        	graph = Graph(n);
			//std::srand ( unsigned ( std::time(0) ) );
			std::vector<Index> perm(n);
			std::map<Index, Index> mapped_id;

			if((bool)options.get("permute")){
				std::iota(perm.begin(), perm.end(), 0);
				std::random_shuffle ( perm.begin(), perm.end() );
				for(int i=0;i<perm.size();i++){
					mapped_id[perm[i]] = i;
				}
			}

			std::string line;
			int start =0;
			Index from_node =0;
			//create graph
			while(std::getline(in,line)){
				line = Util::trim(line);
				if(line[0]=='#' || line[0]=='%' || line[0]=='p') continue;
				if(start ==0){
					start++; 
					continue;
				}
				std::vector<std::string> vertices = Util::split(line, ' ');
				edgeweight w(1.0);
				for(std::string& v: vertices){
					Index n = std::stoi(v);
					if(!(bool)options.get("zero-Index")){
						n = std::stoi(v) -1;
					}
					if((bool)options.get("permute")){
						graph.addEdge(mapped_id[from_node],mapped_id[n], w);
						if((bool)options.get("digraph")){graph.addEdge(mapped_id[from_node],mapped_id[n], w);}
					}else{
						graph.addEdge(from_node,n, w);
						if((bool)options.get("digraph")){graph.addEdge(from_node,n, w);}
					}
				}
				from_node++;
			}

			return graph;
        }
    private:
    	Graph graph;
    	optparse::Values options;
};
#endif
