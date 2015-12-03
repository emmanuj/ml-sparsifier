#ifndef NEIGHBOR_H
#define NEIGHBOR_H

#include "globals.h"

class Neighbor{
public:
	Neighbor(const Index& nodeid, const edgeweight& w):
	nodeId(nodeid),
	weight(w),
	Ro(),
	pij(0),
	algebraicDist(0),
	isdelete(0),
	pos(-1),
	invPos(-1),
	fineEdges()
	{}

	Neighbor(const Neighbor& n):
	nodeId(n.nodeId),
	weight(n.weight),
	Ro(n.Ro),
	pij(n.pij),
	algebraicDist(n.algebraicDist),
	isdelete(n.isdelete),
	pos(n.pos),
	invPos(n.invPos),
	fineEdges(n.fineEdges)
	{
		//std::cout<<"Original value: "<<n.isdelete<<std::endl;
	}

	Neighbor& operator=(const Neighbor& n){
		if(this != &n){
			nodeId = n.nodeId;
			weight = n.weight;
			Ro =n.Ro;
			algebraicDist = n.algebraicDist;
			pij = n.pij;
			isdelete = n.isdelete;
			pos = n.pos;
			invPos = n.invPos;
			fineEdges = n.fineEdges;
		}
		return *this;
	}

	bool operator==(const Neighbor& n) const {
		return nodeId == n.nodeId;
	}

	Index nodeId; //tail id
	edgeweight weight; //weight of the connection.
	std::vector<double> Ro; //list of Roij for each R in computing algebraic distances.
	double pij;
	double algebraicDist;
	bool isdelete;
	long pos;
	long invPos;
	std::vector<std::pair<Index, Index> > fineEdges;
};

#endif
