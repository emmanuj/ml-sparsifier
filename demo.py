#!/usr/bin/python
class Graph:
    def __init__(self, n):
        self.numnodes = n
    def size(self):
        return self.numnodes
    def coarsen(self, frac):
        print "coarsening... Graph reduced to", (self.numnodes * frac), "nodes"
        return Graph(self.numnodes * frac)
    def __str__(self):
        return "graph [numnodes = "+str(self.numnodes)+"]"
    __repr__=__str__

def uncoarsen(G, g):
    print "Uncoarsening", G, g
    #delete from G edges that formed edges marked in g
    #if not fine graph
        #compute algebraic distance
        #mark edges for deletion
        
def ML(G):
    if(G.size() <= 20):
        print "reached coarsest level"
        print G
    else:
        g = G.coarsen(0.5)
        ML(g)
        uncoarsen(G,g)

if __name__ =='__main__':
    g = Graph(100)
    ml(g)
        
    