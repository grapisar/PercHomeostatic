import matplotlib.pyplot as plt
from sys import argv
import numpy as np
from deg_class import AdaptiveDegradation
from get_input_params import input_params

def main():

    ins,dist,SAMPLES,outstring,plot,check = input_params()   

    if(check): 
        
        TOLLK = 1e-12
        #yy = [0,0.25,0.50,0.75]
        yy = np.linspace(0,0.9,10)

        adaptive_deg = AdaptiveDegradation(ins,dist,SAMPLES,outstring,plot,TOLLK)

        G  = adaptive_deg.initialize_network()
        ww = adaptive_deg.initialize_w_ax(G)
        kk = adaptive_deg.initialize_k_ax(G)
        pk0in = adaptive_deg.in_degree_distribution(G,kk)
        pw = adaptive_deg.weight_distribution(ww)
        pwk = adaptive_deg.weight_degree_distribution(pk0in,pw,kk)

        PW_noad_T,WW_noad = adaptive_deg.degradation_no_adaptive(yy,G,ww,pw)
        PW_ad_T,WW_ad = adaptive_deg.degradation_adaptive_multi(yy,G,pwk,kk,pk0in,ww,pw,TOLLK)


        OUT = []
        OUT.append(PW_noad_T)
        OUT.append(WW_noad)
        OUT.append(PW_ad_T)
        OUT.append(WW_ad)
        OUT.append(yy)
        OUT.append(ww)

        OUT = np.array(OUT,dtype='object')
        np.save(outstring,OUT)       


        if plot == 'plot' :
            for i in range(0,len(PW_noad_T)):
                plt.figure()
                plt.hist(WW_noad[i],bins=int(np.sqrt(len(WW_noad[i]))),density=True)
                plt.title(r'$y = %.3f$ No Response' % yy[i])
                plt.xlabel(r'$x$')
                plt.ylabel(r'$w(x)$')
                plt.xlim(0,10)
                plt.plot(ww,PW_noad_T[i])
                plt.show()
                plt.close()

            for i in range(0,len(PW_ad_T)):
                plt.figure()
                plt.hist(WW_ad[i],bins=int(np.sqrt(len(WW_ad[i]))),density=True)
                plt.title(r'$y = %.3f$ Homeostatic Response' % yy[i])
                plt.xlabel(r'$x$')
                plt.ylabel(r'$w(x)$')
                plt.xlim(0,10)
                plt.plot(ww,PW_ad_T[i])
                plt.show()
                plt.close()

if __name__ == "__main__":
    main()