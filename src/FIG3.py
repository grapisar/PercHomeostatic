import matplotlib.pyplot as plt
from sys import argv
import numpy as np
from deg_class import AdaptiveDegradation
from get_input_params import input_params
import seaborn as sns

def plotFIG3(ins):
        
        instring = ins + ".npy"
        data = np.load(instring,allow_pickle=True)

        WGCC_noad_T = data[0]
        WGCC_noad = data[1]
        WGCC_T = data[2]
        WGCC = data[3]
        yy = data[4]
        yy2 = data[5]

        low_quant = 0.25
        high_quant = 0.75

        WGCC_noad = np.array(WGCC_noad)
        WGCC = np.array(WGCC)
        W_noad_plot = []
        W_plot = []

        for i in range(0,len(yy2)):
            W_noad_plot.append([np.mean(WGCC_noad[:,i]),np.quantile(WGCC_noad[:,i],low_quant),np.quantile(WGCC_noad[:,i],high_quant)])
        for i in range(0,len(yy)):
            W_plot.append([np.mean(WGCC[:,i]),np.quantile(WGCC[:,i],low_quant),np.quantile(WGCC[:,i],high_quant)])

        W_noad_plot = np.array(W_noad_plot)
        W_plot = np.array(W_plot)

        sns.set_context("paper", font_scale=1.5, rc={"lines.linewidth": 2.5})

        lower_error = W_noad_plot[:,0] - W_noad_plot[:,1]
        upper_error = W_noad_plot[:,2] - W_noad_plot[:,0]
        asymmetric_error = [lower_error, upper_error]
        plt.errorbar(yy2, W_noad_plot[:,0], yerr=asymmetric_error, fmt='o', color='#ff0000', markersize=8,zorder=5,label=r'No Response')
        plt.plot(yy2,WGCC_noad_T,'kx',zorder=10,markersize=6,markeredgewidth=1.8)
        plt.plot(yy2,WGCC_noad_T,'k--',zorder=0,alpha=0.2)

        lower_error2 = W_plot[:,0] - W_plot[:,1]
        upper_error2 = W_plot[:,2] - W_plot[:,0]
        asymmetric_error2 = [lower_error2, upper_error2]
        plt.errorbar(yy, W_plot[:,0], yerr=asymmetric_error2, fmt='s', color='#0099ff', markersize=7,zorder=5,label=r'Homeostatic Response')
        plt.plot(yy,WGCC_T,'k+',zorder=10,markersize=8,markeredgewidth=1.8)
        plt.plot(yy,WGCC_T,'k-.',zorder=0,alpha=0.2)

        plt.xlabel(r'$y$',fontsize = 20)
        plt.ylabel(r'$S_W$',fontsize = 20)
        plt.xlim(0,2)
        plt.legend()

        plt.tight_layout()
        sns.despine()

        plt.show()




def main():

    ins,dist,SAMPLES,outstring,plot,check = input_params()   

    if(check): 

        TOLLK = 1e-12
        yy = np.linspace(0,2,20)
        yy2 = np.linspace(0,0.95,20)
        M = 100

        adaptive_deg = AdaptiveDegradation(ins,dist,SAMPLES,outstring,plot,TOLLK)

        G  = adaptive_deg.initialize_network()
        ww = adaptive_deg.initialize_w_ax(G)
        kk = adaptive_deg.initialize_k_ax(G)
        pk0in = adaptive_deg.in_degree_distribution(G,kk)
        pk0out = adaptive_deg.out_degree_distribution(G)
        pw = adaptive_deg.weight_distribution(ww)
        pwk = adaptive_deg.weight_degree_distribution(pk0in,pw,kk)

        WGCC_noad_T,WGCC_noad = adaptive_deg.degradation_no_adaptive_WGCC(M,yy2,G,ww,pw,pk0in,TOLLK)
        WGCC_T,WGCC = adaptive_deg.degradation_adaptive_multi_WGCC(M,yy,G,pwk,kk,pk0in,pk0out,ww,pw,TOLLK)

        OUT = []
        OUT.append(WGCC_noad_T)
        OUT.append(WGCC_noad)
        OUT.append(WGCC_T)
        OUT.append(WGCC)
        OUT.append(yy)
        OUT.append(yy2)

        OUT = np.array(OUT,dtype='object')
        np.save(outstring,OUT)

        if plot == 'plot':
            plotFIG3(outstring)

if __name__ == "__main__":
    main()