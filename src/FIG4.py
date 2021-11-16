import matplotlib.pyplot as plt
from sys import argv
import numpy as np
from deg_class import AdaptiveDegradation
from get_input_params import input_params
import seaborn as sns

def plotFIG4(ins):

    instring = ins + ".npy"
    data = np.load(instring,allow_pickle=True)

    MEANSW = data[0]
    AVGW_T = data[1]
    MEANSK = data[2]
    AVGK_T = data[3]
    RHO = data[4]
    RHO_T = data[5]
    yy1 = data[6]
    ww = data[7]

    sns.set_context("paper", font_scale=1.5, rc={"lines.linewidth": 2.5})

    plt.figure()	
    fig, ax1 = plt.subplots()
    ax1.set_xlabel(r'$y$', fontsize = 20)
    ax1.set_ylabel(r'$\bar{k}(y)$', fontsize = 20)
    plt.plot(yy1,MEANSK,'o',color='#ff0000',markersize=8)
    ax1.plot(ww,AVGK_T,'k-',alpha=0.95)
    ax2 = ax1.twinx() 
    ax2.set_ylabel(r'$\bar{w}(y)$', fontsize = 20)  
    ax2.plot(yy1,MEANSW,'^',color='#0099ff',markersize=8)
    ax2.plot(ww,AVGW_T,'k--',alpha=0.95)
    ax2.set_ylim(0.5,3)
    plt.xlim(0,0.99)
    plt.tight_layout()
    plt.show()
    plt.close()


    plt.figure()
    plt.plot(yy1,RHO,'ro',markersize=8)
    plt.plot(ww,RHO_T,'k-',alpha=0.95)
    plt.xlabel(r'$y$', fontsize = 20)
    plt.ylabel(r'$\langle wk \rangle _y - \langle w \rangle _y \langle k \rangle _y$', fontsize = 20)
    plt.tight_layout()
    plt.xlim(0,1)
    sns.despine()

    plt.show()

def main():

    ins,dist,SAMPLES,outstring,plot,check = input_params()   

    if(check): 
        TOLLK = 1e-12
        yy = np.linspace(0,1,50)

        adaptive_deg = AdaptiveDegradation(ins,dist,SAMPLES,outstring,plot,TOLLK)
        
        G  = adaptive_deg.initialize_network()
        ww = adaptive_deg.initialize_w_ax(G)
        kk = adaptive_deg.initialize_k_ax(G)
        pk0in = adaptive_deg.in_degree_distribution(G,kk)
        pw = adaptive_deg.weight_distribution(ww)

        MEANSW,MEANSK,RHO,AVGW_T,AVGK_T,RHO_T,yy1 = adaptive_deg.degradation_adaptive_single(yy,G,ww,pw,kk,pk0in)

        OUT = []
        OUT.append(MEANSW)
        OUT.append(AVGW_T)
        OUT.append(MEANSK)
        OUT.append(AVGK_T)
        OUT.append(RHO)
        OUT.append(RHO_T)
        OUT.append(yy1)
        OUT.append(ww)

        OUT = np.array(OUT,dtype='object')
        np.save(outstring,OUT)

        if plot == 'plot':
            plotFIG4(outstring)

if __name__ == "__main__":
    main()