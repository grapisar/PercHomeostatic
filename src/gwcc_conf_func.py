import numpy as np
from numpy.polynomial.polynomial import polyval, polyval2d, polyder
from numpy import linalg as LA

def max_out_deg(G):
	max = 0
	for n in G.nodes():
		d = G.out_degree(n)
		if(d > max):
			max = d
	return max

def max_in_deg(G):
	max = 0
	for n in G.nodes():
		d = G.in_degree(n)
		if(d > max):
			max = d
	return max	



def multi_deg(G):
	M = max(max_out_deg(G),max_in_deg(G))

	U = np.zeros(shape=(M+1,M+1))

	for n in G.nodes():
		do = G.out_degree(n)
		di = G.in_degree(n)
		U[di][do] += 1
	
	U = U / np.sum(U)

	return U


def polyder2d(U,d):
	dU = np.zeros(shape=(len(U[0]),len(U[:,0])))
	if(d==1):
		for i in range(0,len(U[0])):
			der = list(polyder(U[:,i]))
			der.append(0)
			dU[:,i] = np.array(der)
	if(d==2):
		for i in range(0,len(U[:,0])):
			der = list(polyder(U[i]))
			der.append(0)
			dU[i] = np.array(der)
			
	return dU
		
def fixed_point_directed(U):

	Ui = polyder2d(U,1)
	Uo = polyder2d(U,2)
	
	Ui = Ui / np.sum(Ui)
	Uo = Uo / np.sum(Uo)
	
	F_tol = 1e-5
	
#	xi = 1
	Wi = 0.5
	Wo = 0.5
	
	k = 0
	nrm = 1
	
	while(nrm > F_tol):
		
		Wi_new = polyval2d(Wo,Wi,Ui)
		Wo_new = polyval2d(Wo,Wi,Uo)
		
		k += 1
		
		nrm = LA.norm( Wi - Wi_new ) + LA.norm( Wo - Wo_new )
		
		Wo = Wo_new
		Wi = Wi_new
		
		#if(k%10 == 0):
		#	print(nrm)
			
	W = polyval2d(Wo,Wi,U)
	
	return 1-W
		

