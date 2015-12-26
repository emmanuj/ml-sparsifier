#include "util.h"
#include <sstream>
#include <algorithm>
#include <functional>
#include <cctype>
#include <locale>
#include "graph.h"
#include "node.h"
#include "OptionParser.h"
#include <iomanip>
#include "graphutil.h"
#include <ctime>
#include "neighbor.h"
#include <list>
#include "etimer.h"
#include <omp.h>

class SortByFutureVolume
{
    public:
        bool operator() (const Node& lhs, const Node& rhs){
           return lhs.futureVolume() > rhs.futureVolume();
        }
    private:
};

std::vector<Index> GraphUtil::createSeeds(Graph &g, optparse::Values& options){
    //std::cout<<"\tCreating seed nodes"<<std::endl;
	std::vector<Index> c_nodes;
    std::list<Node> _nodes(g.getNodes().begin(), g.getNodes().end());
    //calculate future volumes
    //initially F = V
    std::list<Node>::iterator itr = _nodes.begin();
    #pragma omp parallel
    {
        #pragma omp for
        for(itr = _nodes.begin(); itr!=_nodes.end(); itr++ ){
            const Index nid = (*itr).getNodeId();
            std::vector<Neighbor>& neighbors = g.neighbors(nid);
            const double vol_i = g.getNode(nid).V();
            double j_term=0;
            for(unsigned j = 0;j<neighbors.size();j++){
                const double vol_j = g.getNode(neighbors[j].nodeId).V();
                j_term+=(vol_j * std::min(1.0, neighbors[j].weight/g.weightedDegree(neighbors[j].nodeId)));
            }
            g.setFutureVolume(nid, j_term + vol_i);
            (*itr).setFutureVolume(j_term + vol_i);
        }
    }

    //moves graph nodes to C nodes and remove those nodes from F
    const double thresh = (int)options.get("mu") * averageVolume(g);

    while(itr != _nodes.end()){
        const Index nid = (*itr).getNodeId();
        if(g.getNode(nid).futureVolume() > thresh){
            c_nodes.push_back(nid);
            g.setCoarse(nid, true);
            itr = _nodes.erase(itr);
        }else{
            ++itr;
        }
    }
    

    //std::cout<<"# c nodes after first FV computation: "<<c_nodes.size()<<std::endl;
    //recompute future volumes of the remaining fine graph.getNodes()
    #pragma omp parallel
    {   
        #pragma omp for
        for(itr = _nodes.begin(); itr!=_nodes.end(); ++itr ){
            const Index nid = (*itr).getNodeId();
            const std::vector<Neighbor>& neighbors = g.neighbors(nid);
            const double vol_i = g.getNode(nid).V();
            double j_term=0;
            for(unsigned j = 0;j<neighbors.size();j++){
                const double vol_j = g.getNode(neighbors[j].nodeId).V();
                double sumFWeights = g.weightedDegree(neighbors[j].nodeId) - g.getNode(neighbors[j].nodeId).getWeightedCDegree();
                j_term+=(vol_j * std::min(1.0, neighbors[j].weight/sumFWeights));
            }
            g.setFutureVolume(nid, j_term+vol_i);
            (*itr).setFutureVolume(j_term + vol_i);
        }
    }

    //sort f in descending order of V. move cnodes to the end
    _nodes.sort(SortByFutureVolume());
    itr = _nodes.begin();
    while(itr != _nodes.end()){
        const double T = g.getNode((*itr).getNodeId()).getWeightedCDegree()/(*itr).getWeightedDegree();
        if( T <= (double)options.get("q")){
            c_nodes.push_back((*itr).getNodeId());
            g.setCoarse((*itr).getNodeId(), true);
            itr = _nodes.erase(itr);
        }else{
            ++itr;
        }
    }
    //std::cout<<"# c nodes after second FV computation: "<<c_nodes.size()<<std::endl;
    ////std::cout<<"# Seeds: "<<c_nodes.size()<<std::endl;
    if(c_nodes.size() <10){
        while(c_nodes.size()<10 && _nodes.size() > 0){
            Index n = (*_nodes.begin()).getNodeId();
            c_nodes.push_back(n);
            g.setCoarse(n, true);
            _nodes.erase(_nodes.begin());
        }
    }

    //std::cout<<"Total # of C nodes found: "<<c_nodes.size()<<std::endl;

    return c_nodes;
}

void GraphUtil::filterConnectionsAD(Node& n, optparse::Values& options, int level){
    std::vector<Neighbor>& cnbs = n.cNeighbors();
    if(cnbs.size()>(int)options.get("r")){
        sort(cnbs.begin(), cnbs.end(),[=](Neighbor a, Neighbor b){
            if(a.pij == 0 && b.pij == 0) return false;
            if(a.pij > 0 && b.pij == 0) return true;
            if(a.pij == 0 && b.pij > 0) return false;
            return a.algebraicDist > b.algebraicDist;
        });

        for(int i=(int) options.get("r");i<cnbs.size();i++){
            cnbs[i].pij = 0.0;
        }

    }

    //normalize it
    double sum =0;
    for(int i=0;i<(int) options.get("r") && i< cnbs.size();i++){
        sum+=cnbs.at(i).pij;
    }

    if(sum != 0){
        for(int i=0;i<(int) options.get("r") && i< cnbs.size();i++){
            cnbs.at(i).pij = cnbs.at(i).pij/sum;
        }
    }
}
void GraphUtil::filterConnections(Node& n, optparse::Values& options){

    std::vector<Neighbor>& cnbs = n.cNeighbors();
    if(cnbs.size()>(int)options.get("r")){
        sort(cnbs.begin(), cnbs.end(),[=](Neighbor a, Neighbor b){
            return a.pij > b.pij;
        });

        for(int i=(int) options.get("r");i<cnbs.size();i++){
            cnbs[i].pij = 0.0;
        }

    }
    //normalize it
    double sum =0;
    for(int i=0;i<(int) options.get("r") && i< cnbs.size();i++){
        sum+=cnbs[i].pij;
    }

    if(sum != 0){
        for(int i=0;i<(int) options.get("r") && i< cnbs.size();i++){
            cnbs[i].pij = cnbs[i].pij/sum;
        }
    }
}

Graph GraphUtil::coarsenGraph(Graph &g, optparse::Values& options){
    //ETimer timer, m_timer;
    //if(g.getLevel() == 2)exit(0);
    //timer.start();
    //m_timer.start();
    //std::cout<<"Computing algebraic distance "<<std::endl;
    g.computeAlgebraicDistances(options); //This should always come first before seed creation.
    //TODO: consider using a list of pointers for cNeighbors vector.
    //timer.stop("Completed alg distance");
    std::vector<Index> c_nodes = createSeeds(g, options);

    if(c_nodes.size()==0){
        //std::cout<<"Error: Could not coarsen graph. Unable to create seed nodes. Graph size: "<<g.getSize()<<std::endl;
        return g;
    }

    //compute interpolation
    //std::cout<<"Computing Pij values "<<std::endl;
    //timer.start();
    #pragma omp parallel
    {
        #pragma omp for
        for(int i=0;i<g.getSize();i++){
            if(g.getNode(i).isC()) continue;
            Index idx = g.getNode(i).getNodeId();
            for(int j=0;j<g.getNode(i).cNeighbors().size();j++){
                Neighbor& nb = g.getNode(i).cNeighbors()[j];
                nb.pij = nb.weight/g.getNode(i).getWeightedCDegree();
            }

            //filterConnections(g.getNode(i), options);
            filterConnectionsAD(g.getNode(i), options, g.getLevel());
        }
    }
    
    //timer.stop("PIj computation completed");


    //std::cout<<"Creating list of graph edges"<<std::endl;
    //timer.start();
    std::list<Edge> edges; //all edges in the graph
    std::vector<std::vector<CoarseEdge>> cEdges(c_nodes.size());
    std::vector<Index> mapper(g.getSize());

    #pragma omp parallel
    {
        #pragma omp for
        for(int i =0;i<g.getSize();i++){
            std::vector<Neighbor>& nbs = g.neighbors(g.getNode(i).getNodeId());
            for(unsigned j=0;j<nbs.size();j++){
                if(g.getNode(i).getNodeId() > nbs[j].nodeId){
                    Edge ed(g.getNode(i).getNodeId(), nbs[j].nodeId, nbs[j].weight);
                    ed.secondPos = nbs[j].pos;
                    edges.push_back(ed);
                }
            }
        }
        //timer.stop("edge list creation ");

        //std::cout<<"Creating coarse edges (using tripples)"<<std::endl;
        //timer.start();
        
        #pragma omp for
        for(int i=0;i<c_nodes.size();i++){
            mapper[c_nodes[i]] = i;
        }

        //std::vector<std::vector< std::pair<Index, edgeweight> >> cEdges(c_nodes.size());
        
        //loop through all edges and create the connections
        std::list<Edge>::iterator it = edges.begin();

        #pragma omp for
        for(it = edges.begin(); it!= edges.end();it++){
            Edge& e = *it;
            Node& from = g.getNode(e.first);
            Node& to = g.getNode(e.second);
            Index idx1 = e.first;
            Index idx2 = e.second;

            double wij = e.weight;
            if(from.isC() && to.isC()){ //C-C
    			////std::cout<<"coarse to coarse"<<std::endl;
                if(idx1 > idx2){ //since the graph is symmetric only update one time.
                    cEdges[mapper[idx1]].push_back(CoarseEdge(mapper[idx2], wij, std::make_pair(idx1, e.secondPos)));
                }else{
                    cEdges[mapper[idx2]].push_back(CoarseEdge(mapper[idx1], wij, std::make_pair(idx1, e.secondPos)));
                }
            }else if((!from.isC()) && (to.isC())){//F-C

                std::vector<Neighbor>& ni = g.coarseNeighbors(idx1); //get coarse neigbors of the fine node
                for(int i=0;i<ni.size();i++){
                    if(ni[i].nodeId != idx2){
                        double w = ni[i].pij * wij * 1;
                        if(w == 0) continue;
    					////std::cout<<"fine to coarse "<<ni[i].pij<<std::endl;
                        if(ni[i].nodeId > idx2){
                            cEdges[mapper[ni[i].nodeId]].push_back(CoarseEdge(mapper[idx2], w, std::make_pair(idx1, e.secondPos)));
                        }else{
                            cEdges[mapper[idx2]].push_back(CoarseEdge(mapper[ni[i].nodeId], w, std::make_pair(idx1, e.secondPos)));
                        }
                    }
                }
            }else if((from.isC()) && (!to.isC())){//C-F
                std::vector<Neighbor>& ni = g.coarseNeighbors(idx2);
                for(int i=0;i<ni.size();i++){
                    if(ni[i].nodeId != idx1){
                        double w = ni[i].pij * wij * 1;
                        if(w == 0) continue;
    					////std::cout<<"coarse to fine "<< ni[i].pij<<std::endl;
                        if(ni[i].nodeId > idx1){
                            cEdges[mapper[ni[i].nodeId]].push_back(CoarseEdge(mapper[idx1], w, std::make_pair(idx1, e.secondPos)));
                        }else{
                            cEdges[mapper[idx1]].push_back(CoarseEdge(mapper[ni[i].nodeId], w, std::make_pair(idx1, e.secondPos)));
                        }
                    }
                }
            }else{//F-F
                std::vector<Neighbor>& ni = g.coarseNeighbors(idx1);
                std::vector<Neighbor>& nj = g.coarseNeighbors(idx2);
                for(int i=0;i<ni.size();i++){
                    for(int j=0;j<nj.size();j++){
                        if(ni[i].nodeId != nj[j].nodeId){
                            double w = ni[i].pij * wij * nj[j].pij;
                            if(w == 0) continue;
    						////std::cout<<"fine to fine"<<std::endl;
                            if(ni[i].nodeId > nj[j].nodeId){
                                cEdges[mapper[ni[i].nodeId]].push_back(CoarseEdge(mapper[nj[j].nodeId], w, std::make_pair(idx1, e.secondPos)));
                            }else{
                                cEdges[mapper[nj[j].nodeId]].push_back(CoarseEdge(mapper[ni[i].nodeId], w, std::make_pair(idx1, e.secondPos)));
                            }
                        }
                    }
                }
            }
        }
    }
    //timer.stop("Coarse edges creation completed in ");
    //std::cout<<"Creating coarse graph..."<<std::endl;
    //timer.start();
    Graph coarseGraph(c_nodes, cEdges, g, mapper, g.getLevel()+1);
    //timer.stop("Coarse graph creation completed in");
    //m_timer.stop("coarsening from level "+ std::to_string(g.getLevel()) +" to level "+ std::to_string(g.getLevel() +1));
    return coarseGraph;
}


//computes average future volume of graph
double GraphUtil::averageVolume(Graph& g){
    double sum = 0.0;
    std::vector<Node> nodes = g.getNodes();
    for(int i = 0; i<nodes.size(); i++){
        sum += nodes[i].futureVolume();
    }
    return sum/(nodes.size()*1.0);
}

std::vector<Node> GraphUtil::neighborsInSet(Index i, Graph &g, std::vector<Node> v_set) {
    std::vector<Neighbor> neighbors = g.neighbors(i);
    std::vector<Node> found;
    for(unsigned j=0;j<neighbors.size();j++){
        if(std::find(v_set.begin(), v_set.end(), g.getNode(neighbors[j].nodeId)) != v_set.end()) {
            found.push_back(g.getNode(neighbors[j].nodeId));
        }
    }

    return found;
}
