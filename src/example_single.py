from deg_analytic_func import *
from deg_numeric_func import *
import numpy as np
import matplotlib.pyplot as plt
from sys import argv
from scipy.special import factorial,binom
import time
import seaborn as sns
from progress.bar import Bar

ins = argv[1]
dist = argv[2]
SAMPLES = int(argv[3])
outstring = argv[4]
plot = argv[5]

G = nx.read_edgelist(ins,create_using=nx.DiGraph())
		 		
### MODEL INITIAL CONDITIONS AND PARAMETERS### 

N = nx.number_of_nodes(G)
print('N =',N)
assign_weights(G,dist)
G_noad = G.copy()

TOLLK = 1e-12

#DEFINITION OF W and K domains#

w_max = 3

dw = w_max / SAMPLES
k_max = max_out_deg(G) + 1
ww = np.arange(0,w_max,dw)
kk = np.arange(0,k_max,1)
print('Max deg =', (k_max-1))

#INITIAL DEGREE DISTRIBUTION#
k_list = list([G.out_degree(n) for n in G.nodes()])
pk = [k_list.count(i)/len(k_list) for i in kk]
pk0out = pk

kin_list = list([G.in_degree(n) for n in G.nodes()])
pk_in = [kin_list.count(i)/len(kin_list) for i in kk]
pk0in = pk_in

U = np.outer(pk,pk_in)

#INITIAL WEIGHT DISTRIBUTION#
if(dist == 'uniform'):
	pw = np.heaviside(ww,1)*np.heaviside(1-ww,1)
if(dist == 'gauss'):
	pw = norm(gaussian(ww,0.5,0.1),dw)

#INITIAL MEAN DEGREE#	
meank = mean_k(pk)

#INITIAL EXCESS DEGREE DISTRIBUTION#
pkex = PFK_EX(pk,kk)

#INITIAL JOINT W, EXCESS K DEGREE DISTRIBUTION#

pwk = []
for k in range(1,len(kk)):
	pwk.append(pw*pkex[k-1])
pwk = np.array(pwk)


ff = [0,0.25,0.50,0.75]

PW_noad_T = []
PW_ad_T = []
WW_noad = []
WW_ad = []

OUT = []

with Bar('No Adaptive',max=len(ff)) as bar:
    for f in ff:
        pw = pw*np.heaviside(ww-f,1)
        pw = norm(pw,dw)
        PW_noad_T.append(pw)

        remove_thresh(G_noad,f,0)
        w_list = [G_noad[e[0]][e[1]]['weight'] for e in G_noad.edges()]
        WW_noad.append(w_list)

        bar.next()


with Bar('Adaptive',max=len(ff)) as bar:
    for f in ff:
        wkpred = pwk
        kpred = pk
        wpred = pw

        #NEW JOINT DIST#

        pwk_ = PFKEX(wkpred,f,kk,ww,dw,kpred,TOLLK,1,1)
        #pwk_ = PFKEX_APPROX(wkpred,f,kk,ww,dw,kpred,TOLLK,1,1,dkmax)
            
        kk = np.arange(0,len(pwk_)+1,1)

        #NEW EXCESS DIST#
        pkex = []
        for k in range(0,len(pwk_)):
            pkex.append(integ(pwk_[k],dw))
        meank = np.sum(pkex)

        #NEW DEG DIST#
        pk = [0]*len(kk)
        for k in range(1,len(kk)): 
            pk[k] = pkex[k-1]/k
        pk[0] = 1 - np.sum(pk)

        #NORMALIZE JOINT AND EXCESS#
        pwk = pwk_/meank
        pkex = pkex/meank

        #NEW W DIST	
        pw = 0
        for q in pwk:
            pw += q

        PW_ad_T.append(pw)

        remove_thresh(G,f,1)
        w_list = [G[e[0]][e[1]]['weight'] for e in G.edges()]
        WW_ad.append(w_list)

        bar.next()

OUT.append(PW_noad_T)
OUT.append(WW_noad)
OUT.append(PW_ad_T)
OUT.append(WW_ad)
OUT.append(ff)
OUT.append(ww)

OUT = np.array(OUT,dtype='object')
np.save(outstring,OUT)

if plot == 'PLOT' :
    for i in range(0,len(PW_noad_T)):
        plt.figure()
        plt.hist(WW_noad[i],bins=int(np.sqrt(len(WW_noad[i]))),density=True)
        plt.title(r'$y = %.3f$ No Response' % ff[i])
        plt.xlabel(r'$x$')
        plt.ylabel(r'$w(x)$')
        plt.xlabel('x')
        plt.plot(ww,PW_noad_T[i])
        plt.show()
        plt.close()

    for i in range(0,len(PW_ad_T)):
        plt.figure()
        plt.hist(WW_ad[i],bins=int(np.sqrt(len(WW_ad[i]))),density=True)
        plt.title(r'$y = %.3f$ Homeostatic Response' % ff[i])
        plt.xlabel(r'$x$')
        plt.ylabel(r'$w(x)$')
        plt.plot(ww,PW_ad_T[i])
        plt.show()
        plt.close()

