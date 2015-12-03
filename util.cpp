#include "util.h"
#include <sstream>
#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>
#include "graph.h"
#include "node.h"
#include <iomanip>

//FUcntions taken from stack overflow
std::vector<std::string> Util::splitv(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}

std::vector<std::string> Util::split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    splitv(s, delim, elems);
    return elems;
}

// trim from start
static inline std::string &ltrim(std::string &s) {
        s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
        return s;
}

// trim from end
static inline std::string &rtrim(std::string &s) {
        s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
        return s;
}

// trim from both ends
std::string Util::trim(std::string &s) {
        return ltrim(rtrim(s));
}

std::vector<double> Util::computeXi(std::vector<Node> &nodes){
    std::vector<double> xiVec(nodes.size());
    for(int i=0; i < nodes.size(); i++){
        //compute vol of idxi
        double xi = 0.5 * nodes[i].V();
        if(i>0){
            xi+= (xiVec[i-1] + (0.5 * nodes[i-1].V()));
        }
        //std::cout<<xi<<" : "<<nodes[i].V()<<std::endl;
        xiVec[i] = xi;
    }
    
    return xiVec;
}

void Util::printCoarseGraph(Graph& g){
   std::cout<<std::endl;
    std::cout<<"Coarse Graph"<<std::endl;
    std::cout<<"From Node\tTo Node\t Weight"<<std::endl;
    int countEdges=0;
    //create coarse graph
    std::vector<Node>::const_iterator itr = g.getNodes().begin();
    for(itr = g.getNodes().begin(); itr!= g.getNodes().end(); ++itr ){
        std::vector<Neighbor> nbs = g.neighbors((*itr).getNodeId());
        for(int j=0;j< nbs.size();j++){
            if((*itr).getNodeId() < nbs[j].nodeId){
                std::cout<<(*itr).getNodeId()<<"\t\t"<<nbs[j].nodeId<<"\t " <<nbs[j].weight<<std::endl;
                countEdges++;
            }
        }
    }
    std::cout<<"Volumes: "<<std::endl;
    for(itr = g.getNodes().begin(); itr!= g.getNodes().end(); ++itr ){
        std::cout<<"Node: "<<(*itr).getNodeId()<<": Vol: "<<(*itr).V()<<std::endl;
    }
    std::cout<<"Coarse Graph {# Edges "<<countEdges<<", # Nodes: "<<g.getNodes().size()<<" }"<<std::endl;
}
void Util::doPumping(Graph& g, std::vector<Node>& nodes){
    sort(nodes.begin(), nodes.end(), SortByYi());
    for(int i=0;i<nodes.size();i++){
        
        double xi = 0.5 * nodes[i].V();
        if(i>0){
            xi+= (nodes[i-1].getX() + (0.5 * nodes[i-1].V()));
        }
        nodes[i].setX(xi);
        g.getNode(nodes[i].getNodeId()).setX(xi);
    }
}
