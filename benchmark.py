#!/usr/bin/env python3

import matplotlib
#matplotlib.use("Agg")

import subprocess #import check_output
import operator
from os import mkdir
from shutil import rmtree
from networkit import *
from pylab import *
import matplotlib.pyplot as plt
from spearman import spearman_rank
import os.path
import sys, traceback
import time
import json

_properties = ["Graph","#edges","#nodes","clust_coef", "clust_dist", "diameter", "betweenness", "pagerank", "deg_centr","modularity", "components"]
def getModularity(G):
	plm = community.PLM(G).run()
	return community.Modularity().getQuality(plm.getPartition(), G)
def numberOfComponents(G):
	return properties.ConnectedComponents(G).run().numberOfComponents()
def clustering(G):
	return properties.ClusteringCoefficient().avgLocal(G)

def averageBetweennessPositions(G):
	positions =[]
	count = 5
	for i in range(0, G.numberOfNodes()):
		positions.append(0)

	for i in range(0,count):
		bt = centrality.ApproxBetweenness(G, 0.2,0.1).run().scores()

		ids = range(0, len(bt))
		nodes_bt = dict(zip(ids, bt))
		sorted_bt = sorted(nodes_bt.items(), key=operator.itemgetter(1))
		pos = G.numberOfNodes()
		for (nodeid, betweennes) in sorted_bt:
		#if pos == 1: print(nodeid, betweennes)
			positions[nodeid] = positions[nodeid] + pos
			pos-=1

	for i in range(0, len(positions)):
		positions[i] = positions[i]/count
	#print(positions[107])
	return positions

def main():

	if(len(sys.argv) < 3): #minimum to run
		print("Invalid number of parameter: Usage: ./benchmark.py name path")
		sys.exit(0)

	setLogLevel("FATAL")

	name = sys.argv[1]
	path = sys.argv[2]

	print("Running ",name)
	start_time = time.time()

	outfile = "reports/_"+name+".txt"

	G = graphio.EdgeListReader(' ',0,"#",True, False).read(path)
	orig_diameter = properties.Diameter.exactDiameter(G)
	orig_clustC = clustering(G)
	numComp = numberOfComponents(G)

	loc_clust_dist = centrality.LocalClusteringCoefficient(G).run().scores()
	#print(" computing degeree cent")
	deg_centr = centrality.DegreeCentrality(G).run().scores()
	#print(" computing page rank")
	page_rank = centrality.PageRank(G).run().scores()
	#print(" computing Betweenness")
	betw = averageBetweennessPositions(G)
	#print(" computing modularity")
	modularity = getModularity(G)
	# fitter = ParameterFitter(G, 20)

	orig_result = ["Original G0", G.numberOfEdges(), G.numberOfNodes(), orig_clustC, 1.0, orig_diameter, 1.0, 1.0, 1.0, modularity, numComp]
	all_g = [orig_result]

	files = ["_coarsest.txt", "_mid.txt", "_finest.txt"]
	titles = ["Coarse G1", "Mid G2", "Fine G3"]
	gnames = dict(zip(files, titles))
	#print(G.numberOfNodes())
	# read other graphs
	for f in files:
		sg = graphio.EdgeListReader(' ',0,"#",True, False).read("out/"+name+f)
		#print(sg.numberOfNodes())
		sg_result = [
			gnames[f], sg.numberOfEdges(), sg.numberOfNodes(), clustering(sg),
			spearman_rank(loc_clust_dist, centrality.LocalClusteringCoefficient(sg).run().scores()), #clust dist
			properties.Diameter.exactDiameter(sg),
			spearman_rank(betw,averageBetweennessPositions(sg)), # rho Betweenness
			spearman_rank(page_rank,centrality.PageRank(sg).run().scores()), # rho page rank
			spearman_rank(deg_centr,centrality.DegreeCentrality(sg).run().scores()), # rhow degree
			getModularity(sg),
			numberOfComponents(sg)
		]

		all_g.append(sg_result)

	# bind properties and write to file

	data = []
	for g in all_g:
		data.append(dict(zip(_properties, g)))

	json.dump(data, open("reports/"+name+".txt","w"))

	print("Completed run for graph: "+name+" in "+ str(time.time() - start_time))

if __name__ == '__main__':
	main()
