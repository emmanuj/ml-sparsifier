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
from scipy.stats import spearmanr
import os.path
import sys, traceback
import time
from param_fitter import ParameterFitter

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

	if(len(sys.argv) < 4): #minimum to run
		print("Invalid number of parameter: Usage: ./generateSparseGraphs name path numnodes")
		sys.exit(0)

	setLogLevel("FATAL")

	name = sys.argv[1]
	path = sys.argv[2]
	num_nodes =sys.argv[3]
	run_args = sys.argv[4:]
	numrun= 2

	print("Running ",name)
	start_time = time.time()
	if os.path.exists("reports/"+name):
		rmtree("reports/"+name)

	mkdir("reports/"+name)

	if os.path.exists("out/"+name):
		rmtree("out/"+name)

	mkdir("out/"+name)

	d_file = "reports/"+name+"/_diameter.txt"
	c_file = "reports/"+name+"/_clust_coef.txt"
	c_file2 = "reports/"+name+"/_clust_coef2.txt"
	comp_file = "reports/"+name+"/_components.txt"
	rho_deg_file = "reports/"+name+"/_degCentrality.txt"
	rho_pag_file = "reports/"+name+"/_pagerank.txt"
	rho_bet_file = "reports/"+name+"/_betweenness.txt"
	mod_file = "reports/"+name+"/_modularity.txt"
	clust_dist_file = "reports/"+name+"/_clust_dist.txt"


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
	r = 0.1
	#read the graphs and run computations
	while(round(r,1) <= 0.9):
		# e = fitter.binarySearchParameter(round(r,1))
		try:
			print("Generating graphs for "+ str(round(r,1)))
			run_cmd = ["./minla", "-i", path, "-n", num_nodes, "--zero-index", "-b", str(numrun), "-s", str(round(r,1)), "--gname", name, ]
			if(len(run_args)):
				for a in run_args:
					run_cmd.append(a)

			#output = check_output(run_cmd)
			subprocess.call(run_cmd)
		except Exception as e:
			print("Process execution failed.",e)
			traceback.print_exc(file=sys.stdout)
			sys.exit(0)

		sparse_g = []
		for i in range(numrun):
			try:
				#print(" Reading graph file for ", round(r,1)," at ",i)
				s_file = "out/"+name+"/"+name+"_"+str(round(r,1))+"_"+str(i)+".txt"
				sparse_g.append(graphio.EdgeListReader(' ',0,"#",True, False).read(s_file))
			except Exception as e:
				print("Failed to read the graph file at "+str(r)+" "+str(i),e)
				sys.exit(0)

		#print(" Computing properties at ", round(r,1))
		#compute diameter
		sumD = 0.0
		sumC = 0.0
		connComp = 0.0
		rho_deg =0
		rho_bet=0
		rho_pag = 0
		rho_clust = 0
		mod =0

		edge_avg = .0
		for g in sparse_g:
			sumD = sumD + properties.Diameter.exactDiameter(g)
			sumC = sumC + clustering(g)
			connComp = connComp + numberOfComponents(g)
			mod += getModularity(g)
			edge_avg+= g.numberOfEdges()

		edge_avg = edge_avg/float(numrun)

		for q in range(numrun):
			sg = sparse_g[q]
			rho_clust += spearmanr(loc_clust_dist, centrality.LocalClusteringCoefficient(sg).run().scores())[0]
			rho_deg+= spearmanr(deg_centr,centrality.DegreeCentrality(sg).run().scores())[0]
			rho_bet+= spearmanr(betw,averageBetweennessPositions(sg))[0]
			rho_pag+= spearmanr(page_rank,centrality.PageRank(sg).run().scores())[0]

		edgeRatio = edge_avg / float(G.numberOfEdges())
		edgeRatio = round(edgeRatio,2)

		avgD = sumD/len(sparse_g) #average Diameter
		avgC = sumC/len(sparse_g) #average Clustering Coefficient
		connComp = (connComp / len(sparse_g))/float(numComp)

		rho_deg = rho_deg/len(sparse_g)
		rho_bet = rho_bet/len(sparse_g)
		rho_pag = rho_pag/len(sparse_g)
		rho_clust = rho_clust/len(sparse_g)

		mod = mod/len(sparse_g)
		mod =mod/modularity

		#print(" Writing to file ", round(r,1))
		with open(d_file,"a") as f:
			f.write(str(round(orig_diameter/avgD,4)) +" "+ str(edgeRatio) +"\n")
		with open(c_file,"a") as f:
			f.write(str(round(avgC - orig_clustC,4)) +" "+ str(edgeRatio) +"\n")
		with open(c_file2,"a") as f:
			f.write(str(round(avgC,4)) +" "+ str(edgeRatio) +"\n")
		with open(comp_file,"a") as f:
			f.write(str(round(connComp,2)) +" "+ str(edgeRatio) +"\n")
		with open(rho_deg_file,"a") as f:
			f.write(str(round(rho_deg,4)) +" "+ str(edgeRatio) +"\n")
		with open(rho_bet_file,"a") as f:
			f.write(str(round(rho_bet,4)) +" "+ str(edgeRatio) +"\n")
		with open(rho_pag_file,"a") as f:
			f.write(str(round(rho_pag,4)) +" "+ str(edgeRatio) +"\n")
		with open(mod_file,"a") as f:
			f.write(str(round(mod,4)) +" "+ str(edgeRatio) +"\n")
		with open(clust_dist_file,"a") as f:
			f.write(str(round(rho_clust,4)) +" "+ str(edgeRatio) +"\n")

		r =round(r,1) + 0.1
		#remove output files
		if os.path.exists("out/"+name):
			rmtree("out/"+name)

		mkdir("out/"+name)


	print("Finalizing ", name)
	#write properties of the full graph
	with open(d_file,"a") as f:
		f.write("1.0 1.0\n")
	with open(c_file,"a") as f:
		f.write("0.0 1.0\n")
	with open(c_file2,"a") as f:
		f.write( str(round(orig_clustC,2)) +" 1.0\n")
	with open(comp_file,"a") as f:
		f.write("1.0 1.0\n")
	with open(mod_file,"a") as f:
		f.write("1.0 1.0\n")
	with open(rho_deg_file,"a") as f:
		f.write("1.0 1.0\n")
	with open(rho_bet_file,"a") as f:
		f.write("1.0 1.0\n")
	with open(rho_pag_file,"a") as f:
		f.write("1.0 1.0\n")
	with open(clust_dist_file,"a") as f:
		f.write("1.0 1.0\n")

	print("Completed run for graph: "+name+" in "+ str(time.time() - start_time))

if __name__ == '__main__':
	main()
