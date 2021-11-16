import matplotlib.pyplot as plt
from sys import argv
import numpy as np
from deg_class import AdaptiveDegradation
from get_input_params import input_params
import seaborn as sns

def plotFIG5(ins):
        
    instring = ins + ".npy"
    data = np.load(instring,allow_pickle=True)

    AVG_S_resp = data[0]
    AVG_S_noresp = data[1]
    AVG_S_resp_multi = data[2]
    AVG_S_noresp_multi = data[3]
    l_max_resp = data[4]
    l_max_noresp = data[5]
    l_max_resp_multi = data[6]
    l_max_noresp_multi = data[7]
    yy = data[8]
    yy2 = data[9]
    ww = data[10]

    low_quant = 0.25
    high_quant = 0.75

    L_max_resp = []
    L_max_noresp = []
    L_max_noresp_multi = []
    L_max_resp_multi = []

    for i in range(0,len(l_max_noresp)):
        L_max_noresp.append([np.mean(l_max_noresp[i]),np.quantile(l_max_noresp[i],low_quant),np.quantile(l_max_noresp[i],high_quant)])
        L_max_resp.append([np.mean(l_max_resp[i]),np.quantile(l_max_resp[i],low_quant),np.quantile(l_max_resp[i],high_quant)])
    for i in range(0,len(yy2)):
        L_max_noresp_multi.append([np.mean(l_max_noresp_multi[:,i]),np.quantile(l_max_noresp_multi[:,i],low_quant),np.quantile(l_max_noresp_multi[:,i],high_quant)])
        L_max_resp_multi.append([np.mean(l_max_resp_multi[:,i]),np.quantile(l_max_resp_multi[:,i],low_quant),np.quantile(l_max_resp_multi[:,i],high_quant)])

    L_max_resp = np.array(L_max_resp)
    L_max_noresp = np.array(L_max_noresp)
    L_max_resp_multi = np.array(L_max_resp_multi)
    L_max_noresp_multi = np.array(L_max_noresp_multi)

    mks = 7.5
    alph_shade = 0.15
    ln_shade = 0.95

    sns.set_context("paper", font_scale=1.5, rc={"lines.linewidth": 2.5})
    plt.figure()
    plt.plot(yy,L_max_noresp[:,0],'o',label = r'$\lambda_{max}$, No response',markersize=mks,color='#ff0000')
    plt.fill_between(yy, L_max_noresp[:,1], L_max_noresp[:,2], facecolor='#ff0000', alpha=alph_shade)
    plt.plot(yy,L_max_resp[:,0],'s',label = r'$\lambda_{max}$, Homeostatic response',markersize=mks,color='#0099ff')
    plt.fill_between(yy, L_max_resp[:,1], L_max_resp[:,2], facecolor='#0099ff', alpha=alph_shade)
    plt.plot(ww,AVG_S_resp,'k-',label = r'$\bar{w}\bar{k}$',alpha = ln_shade)
    plt.plot(ww,AVG_S_noresp,'k--',label = r'$\beta_1\bar{k}$',alpha = ln_shade)
    plt.xlabel('y',fontsize = 20)
    plt.ylabel(r'$\bar{s}$, $\lambda_{max}$',fontsize = 20)
    plt.xlim(0,1)
    plt.legend(fontsize = 10.6)
    plt.tight_layout()
    sns.despine()
    plt.show()
    plt.close()

    plt.figure()
    plt.plot(yy2,L_max_noresp_multi[:,0],'o',label = r'$\lambda_{max}$, No response',markersize=mks,color='#ff0000')
    plt.fill_between(yy2,L_max_noresp_multi[:,1],L_max_noresp_multi[:,2], facecolor='#ff0000', alpha=alph_shade)
    plt.plot(yy2,L_max_resp_multi[:,0],'s',label = r'$\lambda_{max}$, Homeostatic response',markersize=mks,color='#0099ff')
    plt.fill_between(yy2,L_max_resp_multi[:,1],L_max_resp_multi[:,2], facecolor='#0099ff', alpha=alph_shade)
    plt.plot(ww,AVG_S_noresp_multi,'k--',label = r'$\beta_1\bar{k}$',alpha = ln_shade)
    plt.plot(yy2,AVG_S_resp_multi,'k-',label = r'$\bar{w}\bar{k}$',alpha = ln_shade)
    plt.xlim(0,3)
    plt.xlabel('y',fontsize = 20)
    plt.ylabel(r'$\bar{s}$, $\lambda_{max}$',fontsize = 20)
    plt.legend(fontsize = 10.6,bbox_to_anchor = [0.48, 0.1])
    plt.tight_layout()
    sns.despine()
    plt.show()


def main():

    ins,dist,SAMPLES,outstring,plot,check = input_params()   

    if(check): 

        TOLLK = 1e-12
        yy = np.linspace(0,1,50)
        yy2 = np.linspace(0,3,50)
        M = 25

        adaptive_deg = AdaptiveDegradation(ins,dist,SAMPLES,outstring,plot,TOLLK)

        G  = adaptive_deg.initialize_network()
        ww = adaptive_deg.initialize_w_ax(G)
        kk = adaptive_deg.initialize_k_ax(G)
        pk0in = adaptive_deg.in_degree_distribution(G,kk)
        pw = adaptive_deg.weight_distribution(ww)
        pwk = adaptive_deg.weight_degree_distribution(pk0in,pw,kk)
        
        l_max_resp,l_max_noresp,AVG_S_resp,AVG_S_noresp = adaptive_deg.degradation_single_maxeig(M,yy,G,ww,pw,kk,pk0in)
        l_max_noresp_multi,l_max_resp_multi,AVG_S_noresp_multi,AVG_S_resp_multi = adaptive_deg.degradation_multi_maxeig(M,yy2,G,pwk,ww,pw,kk,pk0in,TOLLK)

        OUT = []
        OUT.append(AVG_S_resp)
        OUT.append(AVG_S_noresp)
        OUT.append(AVG_S_resp_multi)
        OUT.append(AVG_S_noresp_multi)
        OUT.append(l_max_resp)
        OUT.append(l_max_noresp)
        OUT.append(l_max_resp_multi)
        OUT.append(l_max_noresp_multi)
        OUT.append(yy)
        OUT.append(yy2)
        OUT.append(ww)
        
        OUT = np.array(OUT,dtype='object')
        np.save(outstring,OUT)

    if plot == 'plot':
        plotFIG5(outstring)

if __name__ == "__main__":
    main()