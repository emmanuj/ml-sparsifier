#ifndef NODE_H
#define NODE_H

#include "globals.h"
#include <iostream>
#include <vector>
#include <algorithm>
#include "neighbor.h"


class Node {
public:
	Node() : volume(1.0), weighted_deg(0), futureVol(1.0), isCoarse(false), nodeId(-1), _neighbors(),
	x(-1), y(-1), v_prime(false), prevId(-1), coarseNeighborhood(), c_weights()
	{ }

	Node(Index n) : volume(1.0), weighted_deg(0), futureVol(1.0), isCoarse(false), nodeId(n),
	 _neighbors(), x(-1), y(-1), v_prime(false), prevId(n), coarseNeighborhood(), c_weights()
	 {}

	Node(const Node& a) : volume(a.volume), 
	weighted_deg(a.weighted_deg),
	futureVol(a.futureVol), 
	isCoarse(a.isCoarse),nodeId(a.nodeId), 
	_neighbors(a._neighbors),x(a.x), y(a.y),
	v_prime(a.v_prime), prevId(a.prevId),
	coarseNeighborhood(a.coarseNeighborhood),
	c_weights(a.c_weights)
	{ }
	
	bool isC()const{return isCoarse;} // return a bool whether the node is coarse
	double V() const{return volume;}//return the Volume
	void setVolume(const double v){ volume = v;}
	void setCoarse(const bool c){ isCoarse = c;}
	void setFutureVolume(const double fvol){ futureVol = fvol;}
	void setNodeId(const Index _id){
		nodeId = _id;
	}
	Index getNodeId() const {return nodeId; }

	void setPrevId(const Index _id) {
		prevId = _id;
	}
	Index getPrevId() const{
		return prevId;
	}

	double futureVolume()const{
		return futureVol;	
	}

	void setX(const double i){ x=i; }
	void setY(const double i){ y=i; }
	void setVPrime(const bool v){ v_prime = v; }

	double getX() const {
		return x;
	}
	double getY()const{
		return y;
	}

	bool inVPrime()const{
		return v_prime;
	}

	void setWeightedDegree(const double wd){
		weighted_deg = wd;
	}

	double getWeightedDegree() const{
		return weighted_deg;
	}

	bool operator<(const Node& n) const {
		return nodeId > n.nodeId;
	}
	bool operator==(const Node& n) const {
		return nodeId == n.nodeId;
	}

	Node& operator=(const Node& node){
		if(this != &node){
			futureVol = node.futureVol;
			volume = node.volume;
			isCoarse = node.isCoarse;
			nodeId = node.nodeId;
			x = node.x;
			y = node.y;
			v_prime = node.v_prime;
			_neighbors = node._neighbors;
			prevId = node.prevId;
			weighted_deg = node.weighted_deg;
			coarseNeighborhood = node.coarseNeighborhood;
			c_weights = node.c_weights;
		}
		return *this;
	}

	long addNeighbor(const Index idx, const double wgt) {
		Neighbor nb(idx, wgt);
		nb.pos = _neighbors.size();
		_neighbors.push_back(nb);
		return nb.pos;
	}

	long addNeighbor(const Index idx, const double wgt, const std::vector<std::pair<Index, Index>> fineEdges) {
		Neighbor nb(idx, wgt);
		nb.pos = _neighbors.size();
		nb.fineEdges = fineEdges;
		_neighbors.push_back(nb);
		return nb.pos;
	}
	long addNeighbor(Neighbor& n) {
		n.pos = _neighbors.size();
		_neighbors.push_back(n);
		return n.pos;
	}

	Neighbor& getNeighborAt(Index pos){
		return _neighbors[pos];
	}

	void removeMarkedEdges(){
		_neighbors.erase(std::remove_if(_neighbors.begin(), _neighbors.end(),[=](Neighbor b){
			//std::cout<<b.isdelete<<std::endl;
            return b.isdelete == true;
        }), _neighbors.end());
	}

	std::vector<Neighbor>& neighbors(){
		return _neighbors;
	}
	std::vector<Neighbor>& cNeighbors(){
		return coarseNeighborhood;
	}

	void addCNeighbor(const Neighbor& n){
		coarseNeighborhood.push_back(n);
		//std::cout<<n.weight<<std::endl;
		c_weights+=n.weight;
	}

	double getWeightedCDegree() const{
		return c_weights;
	}

private:
	double volume;
	double weighted_deg;
	double futureVol;
	bool isCoarse;
	Index nodeId;
	std::vector<Neighbor> _neighbors;
	double x;
	double y;
	bool v_prime;
	Index prevId;
	std::vector<Neighbor> coarseNeighborhood;
	double c_weights;
};

#endif
