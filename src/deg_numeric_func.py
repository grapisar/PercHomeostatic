import networkx as nx
import numpy as np
import random
import copy

def in_strength(G,n):
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

def adapt_after_removal(G,e,a):
	w = G[e[0]][e[1]]['weight']
	G.remove_edge(e[0],e[1])
	n = G.out_degree(e[0])
	if(n > 0):
		e_n = G.out_edges(e[0])
		for ee in e_n:
			G[ee[0]][ee[1]]['weight'] += a*w/n
			
def adapt_after_removal_scale(G,e,a):
	w = G[e[0]][e[1]]['weight']
	s = out_strength(G,e[0])
	G.remove_edge(e[0],e[1])
	n = G.out_degree(e[0])
	if(n > 0):
		e_n = G.out_edges(e[0])
		for ee in e_n:
			G[ee[0]][ee[1]]['weight'] *= (s/(s- a*w))

def max_out_deg(G):
	deg = [G.out_degree(n) for n in G.nodes]
	return max(deg)
	
def min_max_out_weight(G):
	w = [G[e[0]][e[1]]['weight'] for e in G.edges()]
	return min(w), max(w)

def assign_weights(G,dist):
	for e in G.edges():
		x = random.random()
		y = random.normalvariate(0.5,0.1)

		if(dist == 'uniform'):
			G[e[0]][e[1]]['weight'] = x
		if(dist == 'gauss'):
			G[e[0]][e[1]]['weight'] = y

	
	
def average_strength(G):
	if nx.is_empty(G):
		return 0
	else:
		A = nx.to_numpy_matrix(G)
		s = []
		for n in range(0,nx.number_of_nodes(G)):
			s.append(np.sum(A[n]))
		return np.mean(s)

		
def remove_thresh(G,f,a):
	for n in G.nodes():
		nn = [m for m in G.successors(n)]
		to_remove = []
		Delta = 0
		for m in nn:
			if(G[n][m]['weight'] < f):
				to_remove.append([n,m])
				Delta += G[n][m]['weight']
		for e in to_remove:
			G.remove_edge(e[0],e[1])
		nn_new = [m for m in G.successors(n)]
		k = len(nn_new)
		for m in nn_new:
			G[n][m]['weight'] += (a*Delta)/k
		
		
def remove_thresh_scale(G,f,a):
	for n in G.nodes():
		nn = [m for m in G.successors(n)]
		S = out_strength(G,n)
		to_remove = []
		Delta = 0
		for m in nn:
			if(G[n][m]['weight'] < f):
				to_remove.append([n,m])
				Delta += G[n][m]['weight']
		for e in to_remove:
			G.remove_edge(e[0],e[1])
		nn_new = [m for m in G.successors(n)]
		for m in nn_new:
			G[n][m]['weight'] *= (S - (1-a)*Delta)/(S - Delta)

		
def w_w(G):
	a = []
	w1 = []
	w2 = []
	for n in G.nodes:
		if(G.out_degree(n) > 1):
			succ = [nn for nn in G.successors(n)]
			sample = np.random.choice(succ,2,replace=False)
			w1.append(G[n][sample[0]]['weight'])
			w2.append(G[n][sample[1]]['weight'])
	a.append(w1)
	a.append(w2)
	if(len(a)>0):
		C1 = np.corrcoef(a[0],a[1])
		return C1[0][1]
		
def out_in_deg_dist(G):
	k_max = max_out_deg(G) + 1
	kk = np.arange(0,k_max,1)

	k_list = list([G.out_degree(n) for n in G.nodes()])
	pk_out = [k_list.count(i)/len(k_list) for i in kk]

	kin_list = list([G.in_degree(n) for n in G.nodes()])
	pk_in = [kin_list.count(i)/len(kin_list) for i in kk]
	
	return pk_out, pk_in
	
def WCC(G):
	N = nx.number_of_nodes(G)
	out = [len(c)/N for c in sorted(nx.weakly_connected_components(G),key=len, reverse=True)]
	return out
	
def savenet(G,outstring):
	nx.write_weighted_edgelist(G,outstring)	
		
def loadnet(instring):
	G = nx.read_weighted_edgelist(instring,create_using=nx.DiGraph(),nodetype=str)
	return G

