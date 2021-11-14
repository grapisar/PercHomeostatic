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

def w_k_corr(G):
	w_out = []
	k_out = []

	for e in G.edges():
		w_out.append(G[e[0]][e[1]]['weight'])
		k_out.append(G.out_degree(e[0]))	
	
	C = np.corrcoef(w_out,k_out)
	return C[0][1]
	
def w_kin_corr(G):
	w_in = []
	k_in = []

	for e in G.edges():
		w_in.append(G[e[0]][e[1]]['weight'])
		k_in.append(G.in_degree(e[1]))	
	
	C = np.corrcoef(w_in,k_in)
	return C[0][1]

# def adapt_after_removal(G,e):
# 	w = G[e[0]][e[1]]['weight']
# 	G.remove_edge(e[0],e[1])
# 	n = G.out_degree(e[0])
# 	if(n > 0):
# 		e_n = G.out_edges(e[0])
# 		for ee in e_n:
# 			G[ee[0]][ee[1]]['weight'] += w/n
			
# def adapt_after_removal_scale(G,e,a):
# 	w = G[e[0]][e[1]]['weight']
# 	s = out_strength(G,e[0])
# 	G.remove_edge(e[0],e[1])
# 	n = G.out_degree(e[0])
# 	if(n > 0):
# 		e_n = G.out_edges(e[0])
# 		for ee in e_n:
# 			G[ee[0]][ee[1]]['weight'] *= (s/(s- a*w))

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
			
def WCC(G):
	N = nx.number_of_nodes(G)
	out = [len(c)/N for c in sorted(nx.weakly_connected_components(G),key=len, reverse=True)]
	return out
	
def savenet(G,outstring):
	nx.write_weighted_edgelist(G,outstring)	
		
def loadnet(instring):
	G = nx.read_weighted_edgelist(instring,create_using=nx.DiGraph(),nodetype=str)
	return G

