#ifndef GRAPHUTIL__H
#define GRAPHUTIL__H

#include "globals.h"
#include <vector>
#include <string>
#include "node.h"
#include "graph.h"
#include "edge.h"
#include "coarseedge.h"

class GraphUtil{
public:
    static std::vector<Index> createSeeds(Graph &g, optparse::Values& options);
    //static void filterConnections(dense_hash_map<Index, double>& p_row, optparse::Values& options);
    static Graph coarsenGraph(Graph &g, optparse::Values& options);
    static double averageVolume(Graph& g);
    static std::vector<Node> neighborsInSet(Index i, Graph &g, std::vector<Node> v_set);
    //static void filterConnectionsAD(dense_hash_map<Index, double>& p_row, optparse::Values& options, std::map<Index, double>& a_dist);
    static void filterConnectionsAD(Node& n, optparse::Values& options, int l);
    static void filterConnections(Node& n, optparse::Values& options);
private:

};
#endif
