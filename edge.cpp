#include "edge.h"

bool operator<(const Edge& lhs, const Edge& rhs){
	if(lhs.first < rhs.first)
		return true;
	else if(lhs.first == rhs.first)
		if(lhs.second < rhs.second)
			return true;
	return false;
}
//bool operator>(const Edge& lhs, const Edge& rhs){return rhs < lhs;}
//bool operator<=(const Edge& lhs, const Edge& rhs){return !(lhs > rhs);}
//bool operator>=(const Edge& lhs, const Edge& rhs){return !(lhs < rhs);}

bool operator==(const Edge& lhs, const Edge& rhs){
	return (lhs.first == rhs.first && lhs.second == rhs.second)
		|| (lhs.first == rhs.second && lhs.second == rhs.first );
}
