#!/usr/bin/env python3
# Author: Emmanuel John (emmanuj (a) clemson.edu)
from networkit import *
import math

class ParameterFitter:
	def __init__(self, G, niter):
		self.G = G
		self.niter = niter

	def get_edge_ratio(self, e):
		edgeCount = .0
		for u in range(self.G.numberOfNodes()):
			edgeCount += math.ceil(math.pow(self.G.degree(u), e))
		return float(edgeCount/self.G.numberOfEdges())

	def binarySearchParameter(self, target_ratio):
		lo = .0
		hi = 1.0
		for i in range(self.niter):
			mid = float(lo + ((hi-lo)*.5))
			actual_ratio = round(self.get_edge_ratio(mid),2)
			#print(actual_ratio, target_ratio)
			if actual_ratio > target_ratio:
				hi = mid
			else:
				lo = mid
		return lo
