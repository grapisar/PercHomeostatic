import networkx as nx
import numpy as np
import random

def in_strength(G,n):
	eout = G.in_edges(n)
	s = 0
	for e in eout:
		s += G[e[0]][e[1]]['weight']
	return s

def out_strength(G,n):
	eout = G.in_edges(n)
	s = 0
	for e in eout:
		s += G[e[0]][e[1]]['weight']
	return s


def max_out_deg(G):
	deg = [G.out_degree(n) for n in G.nodes]
	return max(deg)

def max_in_deg(G):
	deg = [G.in_degree(n) for n in G.nodes]
	return max(deg)


def assign_weights(G,dist):
	if(dist == 'uniform'):
		for e in G.edges():
			G[e[0]][e[1]]['weight'] = random.random()
	if(dist == 'gauss'):
		for e in G.edges():
			G[e[0]][e[1]]['weight'] = random.normalvariate(0.5,0.1)

		
def remove_thresh(G,f,a):
	for n in G.nodes():
		nn = [m for m in G.predecessors(n)]
		to_remove = []
		Delta = 0
		for m in nn:
			if(G[m][n]['weight'] < f):
				to_remove.append([m,n])
				Delta += G[m][n]['weight']
		for e in to_remove:
			G.remove_edge(e[0],e[1])
		nn_new = [m for m in G.predecessors(n)]
		for m in nn_new:
			G[m][n]['weight'] += (Delta*a)/len(nn_new)
			

