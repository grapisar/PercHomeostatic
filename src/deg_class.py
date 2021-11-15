from deg_analytic_func import *
from deg_numeric_func import *
from gwcc_conf_func import *
import numpy as np
from progress.bar import Bar
import time

class AdaptiveDegradation:
    def __init__(self,ins,dist,SAMPLES,outstring,plot,TOLLK):
        self.instring = ins
        self.weight_distribution_type = dist
        self.samples = SAMPLES
        self.plot_result = plot
        self.tolerance = TOLLK

    def initialize_network(self):
        G = nx.read_edgelist(self.instring,create_using=nx.DiGraph())
        assign_weights(G,self.weight_distribution_type)

        return G

    def initialize_w_ax(self,G):
        w_max = max([in_strength(G,n) for n in G.nodes()])
        dw = w_max / self.samples
        ww = np.arange(0,w_max,dw)
        
        return ww

    def initialize_k_ax(self,G):

        k_max = max_in_deg(G) + 1
        kk = np.arange(0,k_max,1)
        
        return kk

    def in_degree_distribution(self,G,kk):

        kin_list = list([G.in_degree(n) for n in G.nodes()])
        pk0in = [kin_list.count(i)/len(kin_list) for i in kk]
        
        return pk0in

    def out_degree_distribution(self,G):
        kout_max = max_out_deg(G) + 1
        kk = np.arange(0,kout_max,1)
        k_list = list([G.out_degree(n) for n in G.nodes()])
        pk0out = [k_list.count(i)/len(k_list) for i in kk]
        
        return pk0out

    def weight_distribution(self,ww):
        dw = ww[1] - ww[0]
        if(self.weight_distribution_type == 'uniform'):
	        pw = np.heaviside(ww,1)*np.heaviside(1-ww,1)
        if(self.weight_distribution_type == 'gauss'):
	        pw = norm(gaussian(ww,0.5,0.1),dw)
        
        return pw

    def weight_degree_distribution(self,pk,pw,kk):
        pwk = []
        pkex = PFK_EX(pk,kk)
        for k in range(1,len(kk)):
	        pwk.append(pw*pkex[k-1])
        pwk = np.array(pwk)

        return pwk

    def degradation_no_adaptive(self,yy,G,ww,pw):
        dw = ww[1] - ww[0]
        G_noad = G.copy()
        PW_noad_T = []
        WW_noad = []

        with Bar('No Adaptive',max=len(yy)) as bar:
            for y in yy:
                pw = pw*np.heaviside(ww-y,1)
                pw = norm(pw,dw)
                PW_noad_T.append(pw)

                remove_thresh(G_noad,y,0)
                w_list = [G_noad[e[0]][e[1]]['weight'] for e in G_noad.edges()]
                WW_noad.append(w_list)

                bar.next()        
        return PW_noad_T,WW_noad

    def degradation_no_adaptive_WGCC(self,M,yy,G,ww,pw,pk0in,TOLLK):
        dw = ww[1] - ww[0]
        G_noad = G.copy()
        N = nx.number_of_nodes(G)
        WGCC_noad = []
        WGCC_noad_T = []

        with Bar('NO ADAPTIVE - MODEL:', max=len(yy)) as bar:
            t1 = time.time()
            for y in yy:
                pkin_noad = PFK_NOAD(pk0in,y,pw,ww,dw,TOLLK)
                pkout_noad  = pkin_noad
                WGCC_noad_T.append(fixed_point_directed(np.outer(pkout_noad,pkin_noad)))
                bar.next()
            Dt1 = time.time() - t1

        with Bar('NO ADAPTIVE - MC:', max=M) as bar:
            t2 = time.time()
            for i in range(M):
                G_noad2 = G_noad.copy()
                assign_weights(G_noad2,self.weight_distribution_type)
                wgcc_n = []
                for y in yy:
                    remove_thresh(G_noad2,y,0)
                    sw_noad = len(max(nx.weakly_connected_components(G_noad2), key=len))/N
                    wgcc_n.append(sw_noad)
                WGCC_noad.append(wgcc_n)
                bar.next()
            Dt2 = time.time() - t2

        print('Model Time, No Adaptive:',Dt1)
        print('MC Time, No Adaptive:',Dt2)        
        
        return WGCC_noad_T,WGCC_noad


    def degradation_adaptive_multi(self,yy,G,pwk,kk,pk,ww,pw,TOLLK):
        dw = ww[1] - ww[0]
        PW_ad_T = []
        WW_ad = []        

        with Bar('Adaptive',max=len(yy)) as bar:
            for y in yy:
                wkpred = pwk
                kpred = pk
                wpred = pw


                pwk_ = PFKEX(wkpred,y,kk,ww,dw,kpred,TOLLK,1)
                    
                kk = np.arange(0,len(pwk_)+1,1)

                pkex = []
                for k in range(0,len(pwk_)):
                    pkex.append(integ(pwk_[k],dw))
                meank = np.sum(pkex)


                pk = [0]*len(kk)
                for k in range(1,len(kk)): 
                    pk[k] = pkex[k-1]/k
                pk[0] = 1 - np.sum(pk)


                pwk = pwk_/meank
                pkex = pkex/meank


                pw = 0
                for q in pwk:
                    pw += q

                PW_ad_T.append(pw)

                remove_thresh(G,y,1)
                w_list = [G[e[0]][e[1]]['weight'] for e in G.edges()]
                WW_ad.append(w_list)

                bar.next()

        return PW_ad_T,WW_ad


    def degradation_adaptive_multi_WGCC(self,M,yy,G,pwk,kk,pk,pk0out,ww,pw,TOLLK):
        dw = ww[1] - ww[0]
        WGCC_T = []
        WGCC = []
        N = nx.number_of_nodes(G)

        with Bar('ADAPTIVE - MODEL:',max=len(yy)) as bar:
            t1 = time.time()
            for y in yy:
                wkpred = pwk
                kpred = pk
                wpred = pw

                pwk_ = PFKEX_APPROX(wkpred,y,kk,ww,dw,kpred,TOLLK,1,12)
                    
                kk = np.arange(0,len(pwk_)+1,1)

                pkex = []
                for k in range(0,len(pwk_)):
                    pkex.append(integ(pwk_[k],dw))
                meank = np.sum(pkex)

                pk = [0]*len(kk)
                for k in range(1,len(kk)): 
                    pk[k] = pkex[k-1]/k
                pk[0] = 1 - np.sum(pk)

                pkout = PFK_OUT(kpred,wpred,y,ww,dw,kk)
            
                WGCC_T.append(fixed_point_directed(np.outer(pkout,pk)))

                pwk = pwk_/meank
                pkex = pkex/meank


                pw = 0
                for q in pwk:
                    pw += q

                bar.next()
            Dt1 = time.time() - t1

        with Bar('ADAPTIVE - MC:',max=M) as bar:
            t2 = time.time()
            for i in range(M):
                G2 = G.copy()
                assign_weights(G2,self.weight_distribution_type)
                wgcc = []
                for y in yy:
                    remove_thresh(G2,y,1)

                    if not nx.is_empty(G2):	
                        sw = len(max(nx.weakly_connected_components(G2), key=len))/N
                        wgcc.append(sw)
                    else:
                        wgcc.append(0)
                WGCC.append(wgcc)
                bar.next()
            Dt2 = time.time() - t2       



        print('Model Time, Adaptive:',Dt1)
        print('MC Time, Adaptive:',Dt2)

        return WGCC_T,WGCC

    def degradation_adaptive_single(self,yy,G,ww,pw,kk,pk):

        MEANSW = []
        MEANSK = []
        RHO = []
        yy1 = []
        dw = ww[1] - ww[0]
        w_max = ww[-1]

        with Bar('MC',max=len(yy)) as bar:
            t1 = time.time()
            for y in yy:
                
                G2 = G.copy()
                
                remove_thresh(G2,y,1)
                
                if not nx.is_empty(G2):
                    w = []
                    wk = []
                    kn = []
                    for e in G2.edges():
                        w.append(G2[e[0]][e[1]]['weight'])
                        wk.append(G2.in_degree(e[1]) * G2[e[0]][e[1]]['weight'])
                        kn.append(G2.in_degree(e[1]))
                    
                    w_mean = np.mean(w)
                    MEANSW.append(w_mean)
                    
                    wk_mean = np.mean(wk)
                    kn_mean = np.mean(kn)
                
                    kin_list = [G2.in_degree(n) for n in G2.nodes()]
                    k_mean = np.mean(kin_list)
                    
                    MEANSK.append(k_mean)
                    
                    RHO.append(wk_mean - w_mean*kn_mean)
                    
                    yy1.append(y)
                bar.next()
            Dt1 = time.time() - t1

        pkex = PFK_EX(pk,kk)

        t2 = time.time()
        
        k_0 = mean_k(pk)
        kn_0 = np.sum([k*pkex[k-1] for k in range(1,max_in_deg(G) + 1)])

        beta1, beta2 = beta1_2(pw,ww,dw,w_max)
        F = primitive(pw,ww,dw)
        GF = GEN(F,pkex)

        AVGW_T = beta1 + beta2*( (F - GF) / (1 - F) )
        AVGK_T = k_0*(1 - F)
        RHO_T = beta2*(kn_0*GF - (F*(1-GF))/(1-F))
        
        Dt2 = time.time()-t2

        print('MC Time:',Dt1)
        print('Model Time:',Dt2)

        return MEANSW,MEANSK,RHO,AVGW_T,AVGK_T,RHO_T,yy1




