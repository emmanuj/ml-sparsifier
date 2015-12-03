#ifndef COARSEEDGE_H
#define COARSEEDGE_H
class CoarseEdge{
public:
	CoarseEdge():nodeId(-1), weight(0), fineEdge(){}
	CoarseEdge(const Index n, const double w, const std::pair<Index, Index>& fine_e):
		nodeId(n), weight(w), fineEdge(fine_e){}
	CoarseEdge(const CoarseEdge& c): nodeId(c.nodeId),weight(c.weight), fineEdge(c.fineEdge){}
	CoarseEdge& operator=(const CoarseEdge& c){
		if(this !=  &c){
			nodeId = c.nodeId;
			weight = c.weight;
			fineEdge = c.fineEdge;
		}
		return *this;
	}
	Index nodeId;
	double weight;
	std::pair<Index, Index> fineEdge; //first = first nodeId, second = position of the neighbor in the neighbor list
};

#endif
