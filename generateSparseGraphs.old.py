#!/usr/bin/python3
import subprocess #import check_output

from os import mkdir
from shutil import rmtree
from networkit import *
from pylab import *
import matplotlib.pyplot as plt
from scipy.stats import spearmanr
import os.path

setLogLevel("FATAL")
error = 0.01
nSamples = math.ceil(math.log(10) / (error**2)) # fixed confidence of 90%

def getModularity(G):
	plm = community.PLM(G).run()
	return community.Modularity().getQuality(plm.getPartition(), G)
def numberOfComponents(G):
	return properties.ConnectedComponents(G).run().numberOfComponents()
def clustering(G):
	return properties.ClusteringCoefficient().avgLocal(G)
'''

graphs = [
	["fb","/home/emmanuj/graphs/facebook_combined.txt","4039","r"],
	["fb","/home/emmanuj/graphs/facebook_combined.txt","4039","r"],
	["yt","/home/emmanuj/graphs/social_new/youtube.txt","1134890","r"],
	["amazon","/home/emmanuj/graphs/social_new/amazon.txt","334863","b"],
	["dblp","/home/emmanuj/graphs/social_new/dblp.txt","317080","g"]
]
'''
group_name ="fb_amazon_dblp"
levels = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0]
lev = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9]

for gr in graphs:
	name = gr[0]
	path = gr[1]
	num_nodes =gr[2]
	print("Running ",name)
	r = 0.1
	if os.path.exists("reports/"+name):
		rmtree("reports/"+name)

	mkdir("reports/"+name)
	
	if os.path.exists("out/"+name):
		rmtree("out/"+name)
	
	mkdir("out/"+name)

	d_file = "reports/"+name+"/_diameter.txt"
	c_file = "reports/"+name+"/_clust_coef.txt"
	comp_file = "reports/"+name+"/_components.txt"
	rho_deg_file = "reports/"+name+"/_degCentrality.txt"
	rho_pag_file = "reports/"+name+"/_pagerank.txt"
	rho_bet_file = "reports/"+name+"/_betweenness.txt"
	mod_file = "reports/"+name+"/_modularity.txt"

	G = graphio.EdgeListReader(' ',0,"#",True, False).read(path)
	orig_diameter = properties.Diameter.exactDiameter(G)
	orig_clustC = clustering(G)
	numComp = numberOfComponents(G)
	
	print("here: computing degeree cent")
	deg_centr = centrality.DegreeCentrality(G).run().scores()
	print("here: computing page rank")
	page_rank = centrality.PageRank(G).run().scores()
	print("here: computing Betweenness")
	betw = centrality.ApproxBetweenness(G, 0.05, 0.1).run().scores()
	print("here: computing modularity")
	modularity = getModularity(G)
	
	
	s_dia = []
	s_con = []
	s_clust= []
	s_rho_deg =[]
	s_rho_bet =[]
	s_rho_pag = []
	s_mod =[]

	#generate the graphs
	try:
		print("here: Generating graphs")
		run_cmd = ["./minla", "-i", path, "-n", num_nodes, "--digraph", "--zero-index", "-b", "10", "-s", str(round(r,1)), "--gname", name]
		#output = check_output(run_cmd)
		subprocess.call(run_cmd)
	except Exception as e:
		print("Process execution failed.",e)

	#read the graphs and run computations
	while(round(r,1) <= 0.9):
		lstr = []
		for i in range(10):
			try:
				print("here: Reading graph file for ", round(r,1)," at ",i)
				s_file = "out/"+name+"/"+name+"_"+str(round(r,1))+"_"+str(i)+".txt"
				lstr.append(graphio.EdgeListReader(' ',0,"#",True, False).read(s_file))
			except Exception as e:
				print("Failed to read the graph file.",e)

		print("here: Computing properties at ", round(r,1))
		#compute diameter
		sumD = 0.0
		sumC = 0.0
		connComp = 0.0
		rho_deg =0
		rho_bet=0
		rho_pag = 0
		mod =0
		for g in lstr:
			sumD = sumD + properties.Diameter.exactDiameter(g)
			sumC = sumC + clustering(g)
			connComp = connComp + numberOfComponents(g)
			mod += getModularity(g)

		for q in range(5):
			rho_deg+= spearmanr(deg_centr,centrality.DegreeCentrality(g).run().scores())[0]
			rho_bet+= spearmanr(betw,centrality.ApproxBetweenness(g, 0.05, 0.1).run().scores())[0]
			rho_pag+= spearmanr(page_rank,centrality.PageRank(g).run().scores())[0]
		
		avgD = sumD/len(lstr) #average Diameter
		avgC = sumC/len(lstr) #average Clustering Coefficient
		connComp = (connComp / len(lstr))/numComp
		
		rho_deg = rho_deg/5
		rho_bet = rho_bet/5
		rho_pag = rho_pag/5
		mod = mod/len(lstr)
		mod =mod/modularity

		#TODO: add modularity computation
		print("here: Writing to file ", round(r,1))
		with open(d_file,"a") as f:
			f.write(str(round(orig_diameter/avgD,2)) +" "+ str(round(r,1)) +"\n")
		with open(c_file,"a") as f:
			f.write(str(round(avgC - orig_clustC,2)) +" "+ str(round(r,1)) +"\n")
		with open(comp_file,"a") as f:
			f.write(str(round(connComp,1)) +" "+ str(round(r,1)) +"\n")
		with open(rho_deg_file,"a") as f:
			f.write(str(round(rho_deg,3)) +" "+ str(round(r,1)) +"\n")
		with open(rho_bet_file,"a") as f:
			f.write(str(round(rho_bet,3)) +" "+ str(round(r,1)) +"\n")
		with open(rho_pag_file,"a") as f:
			f.write(str(round(rho_pag,3)) +" "+ str(round(r,1)) +"\n")
		with open(mod_file,"a") as f:
			f.write(str(round(mod,3)) +" "+ str(round(r,1)) +"\n")

		r =round(r,1) + 0.1

	print("here: Plotting for ", name)
	#write properties of the full graph
	with open(d_file,"a") as f:
		f.write("1.0 1.0\n")
	with open(c_file,"a") as f:
		f.write("0.0 1.0\n")
	with open(comp_file,"a") as f:
		f.write(str(numComp)+" 1.0\n")
	with open(mod_file,"a") as f:
		f.write("1.0 1.0\n")

	s_dia.append(1.0)
	s_con.append(1.0)
	s_clust.append(0.0)
	s_mod.append(1.0)

	figure(1)
	plot(levels, s_dia,gr[3],label=name)
	title("Diameter")
	xlabel("percentage sparsification")
	ylabel("Orig Diameter/diameter")
	legend()

	figure(2)
	plot(levels, s_clust,gr[3],label=name)
	title("Clustering Coefficient")
	xlabel("percentage sparsification")
	ylabel("cc - orig. cc")
	legend()

	figure(3)
	plot(levels, s_con,gr[3], label=name)
	title("Connected Components")
	xlabel("percentage sparsification")
	ylabel("#connected components/orig #conn components")
	legend()

	figure(4)
	plot(lev, s_rho_bet,gr[3],label=name)
	title("Betweenness Centrality")
	xlabel("percentage sparsification")
	ylabel("rho")
	legend()

	figure(5)
	plot(lev, s_rho_deg,gr[3],label=name)
	title("Degree Centrality")
	xlabel("percentage sparsification")
	ylabel("rho")
	legend()

	figure(6)
	plot(lev, s_rho_pag,gr[3],label=name)
	title("PageRank")
	xlabel("percentage sparsification")
	ylabel("rho")
	legend()

	figure(7)
	plot(levels, s_mod, gr[3],label=name)
	title("Modularity")
	xlabel("percentage sparsification")
	ylabel("modularity/orig modularity")
	legend()
	rmtree("out/"+name)

print("Saving all plots")

figure(1)
savefig("reports/"+group_name+"_diameter.png")
figure(2)
savefig("reports/"+group_name+"_clust.png")
figure(3)
savefig("reports/"+group_name+"_components.png")
figure(4)
savefig("reports/"+group_name+"_betweenness.png")
figure(5)
savefig("reports/"+group_name+"_degCentrality.png")
figure(6)
savefig("reports/"+group_name+"_pagerank.png")
figure(7)
savefig("reports/"+group_name+"_modularity.png")

