#ifndef EDGE_H
#define EDGE_H

#include "node.h"
#include "globals.h"
#include <iostream>
class Edge{
public:
	Edge():
    first(),
    second(),
    weight(),
    algDist(), 
    secondPos()
    {}

	Edge(const Index& node1, const Index& node2, const edgeweight& w):
	first(node1),
	second(node2),
	weight(w),
	algDist(0),
	secondPos()
	{}

	Edge(const Index& node1, const Index& node2, const edgeweight& w, const edgeweight& alg):
	first(node1),
	second(node2),
	weight(w),
	algDist(alg),
	secondPos()
	{}

	Edge(const Edge& e):
	first(e.first),
	second(e.second),
	weight(e.weight), algDist(e.algDist),
	secondPos(e.secondPos){

	}
	Edge& operator=(const Edge& edge){
		if(this != &edge){
			first = edge.first;
			second = edge.second;
			weight = edge.weight;
			algDist = edge.algDist;
			secondPos = edge.secondPos;
		}
		return *this;
	}
	
	friend std::ostream& operator<<(std::ostream& out, const Edge& e){
		out<<e.first<<" "<<e.second<<" "<<e.weight<<" "<<e.algDist;
		return out;
	}

	Index first;
	Index second;
	edgeweight weight;
	edgeweight algDist;
	bool isdelete;
	long secondPos; //index of the second node in the neighbor list of the first

	friend bool operator<(const Edge& lhs, const Edge& rhs);
	friend bool operator>(const Edge& lhs, const Edge& rhs);
	friend bool operator<=(const Edge& lhs, const Edge& rhs);
	friend bool operator>=(const Edge& lhs, const Edge& rhs);
	friend bool operator==(const Edge& lhs, const Edge& rhs);

};

#endif
