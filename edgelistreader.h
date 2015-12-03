#ifndef EDGELISTREADER_H
#define EDGELISTREADER_H

#include "reader.h"
#include "graph.h"
#include <iostream>
#include "util.h"
#include <fstream>
#include <string>
#include <algorithm>
class EdgeListReader: Reader
{
    public:
    	EdgeListReader(const optparse::Values& op):options(op), graph(0){}

        Graph& read(std::ifstream& in, long n){
        	graph = Graph(n);
			//std::srand ( unsigned ( std::time(0) ) );
			std::vector<Index> perm(n);
			std::map<Index, Index> mapped_id;
			std::string dm = (std::string)options.get("d");
			char delim = ' ';
			if(dm == "tab"){
				delim = '\t';
			}

			if((bool)options.get("permute")){
				std::iota(perm.begin(), perm.end(), 0);
				std::random_shuffle ( perm.begin(), perm.end() );
				for(int i=0;i<perm.size();i++){
					mapped_id[perm[i]] = i;
				}
			}

			std::string line;
			//create graph
			while(std::getline(in,line)){
				line = Util::trim(line);
				if(line[0]=='#' || line[0]=='%' || line[0]=='p') continue;
				std::vector<std::string> vertices = Util::split(line, delim);
				Index n1;
				Index n2;
				edgeweight w(1.0);

				if(line[0]=='e'){
					n1 = (std::stoi(vertices[1]));
					n2 =(std::stoi(vertices[2]));
					if(vertices.size() == 4 ){
						w = std::stof(vertices[3]);
					}
				}else{
					n1 = (std::stoi(vertices[0]));
					n2 =(std::stoi(vertices[1]));
					if(vertices.size() == 3 ){
						w = std::stof(vertices[2]);
					}
				}
				if(n1 == n2) continue;//no self loops allowed
				if(!(bool)options.get("zero-Index")){
					n1-=1;
					n2-=1;
				}

				if((bool)options.get("permute")){
					graph.addEdge(mapped_id[n1],mapped_id[n2], w);
				}else{
					graph.addEdge(n1,n2, w);
				}
			}

			return graph;
        }
    private:
    	Graph graph;
    	optparse::Values options;
};
#endif
