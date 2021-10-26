import numpy as np
from scipy.special import factorial, binom
from numpy.fft import fft,ifft

def integ(arr,df):
	return np.sum(arr)*df

def primitive(arr,ww,dw):
	out = []
	for t in ww:
		out.append(integ(arr*np.heaviside(t-ww,1),dw))
	return np.array(out)
	
def norm(arr,dw):
	if(integ(arr,dw) != 0):
		return arr/integ(arr,dw)
	if(integ(arr,dw) == 0):
		return arr

def scalex(arr,fac):
	split = arr[::int(fac)]
	add = np.zeros(len(arr)-len(split))
	return np.concatenate((split,add))

def gaussian(x, mu, sig):
    return np.exp(-np.power(x - mu, 2.) / (2 * np.power(sig, 2.)))
    
	
def conv_norm(a,b,n,df):		#convolution of a*b**n
	fa = fft(a)
	fb = (fft(b)*df)**n
	return np.real(ifft(fa*fb))

def mean_k(pk):
	s = 0
	for k in range(0,len(pk)):
		s += pk[k]*k
	return s
	
def mean_kex(pkex):
	s = 0
	for k in range(0,len(pkex)):
		s += pkex[k]*(k+1)
	return s	
	
def std_k(pk):
	ss = 0
	for k in range(0,len(pk)):
		ss += pk[k]*k*k
	return np.sqrt(ss - (mean_k(pk))**2 )
	
def mean_w(pw,ww,dw):
	return integ(pw*ww,dw)
	
def std_w(pw,ww,dw):
	return np.sqrt( integ(pw*ww*ww,dw) - (mean_w(pw,ww,dw))**2 )

def norm_wk(pwk,ww,kk,dw):
	s = 0
	for k in range(1,len(kk)):
		s += pwk[k-1]
	n = integ(s,dw)
	return pwk/n
	
def exp_wk(pwk,ww,kk,dw):
	s = 0
	for k in range(0,len(kk)):
		s += pwk[k-1]*k
	return integ(s*ww,dw)	

def PFKEX(wkpred,f,kk,ww,dw,kpred,TOLLK,multi,alpha):
	out = []
	for k in range(1,len(kpred)):
		#print(k)
		s = 0
		for j in range(k,len(kk)):
			a_k = norm(wkpred[j-1],dw)
			a = integ(a_k*np.heaviside(ww-f,1),dw)

			wpred_up = wkpred[j-1]*np.heaviside(ww-f,1)
			wpred_up_norm = norm(wpred_up,dw)
			wpredk_down = scalex(wkpred[j-1],k/alpha)*np.heaviside(f - (k*ww)/alpha,1)
			wpredk_down_norm = norm(wpredk_down,dw)

			conv_ = conv_norm(wpred_up,wpredk_down,(j-k),dw)
			conv = norm(conv_,dw)
			
			#PREVENTS DIVERGENCES DUE TO NUMERICAL INSTABILITIES
			if(np.sum(conv > 1e6) != 0):
				conv[conv > 1e6] = 0

			s += binom(j,k)*(a**k)*((1-a)**(j-k))*kpred[j]*conv			

		out.append(s*k)
	if(multi == 1):
		if(len(out) > 1):
			while True:
				a = np.array(out[-1])
				p = integ(a,dw)
				if(p > TOLLK):
					break
				else:
					del out[-1]
			
	return np.array(out)
	
	
def PFK_EX(pfk,kk):
	out = []
	for k in range(1,len(kk)):
		out.append(k*pfk[k])
	out = np.array(out)
	return out/np.sum(out)
	
def PFK_IN(kpred,wpred,f,ww,dw,kk):
	out = []
	a = integ(wpred*np.heaviside(ww-f,1),dw)
	for k in kk:
		sum = 0
		for j in range(k,len(kk)):
			sum += binom(j,k)*(a**k)*((1-a)**(j-k))*(kpred[j])

		out.append(sum)
	
	return(np.array(out))
	
def beta1_2(pw,ww,dw,wmax):
	eps = 1e-5
	b2 = []
	b1 = []
	w0 = integ(ww*pw,dw)
	for t in ww:
		x1 = integ(ww*pw*np.heaviside(t-ww,1),dw)
		x2 = integ(pw*np.heaviside(t-ww,1),dw)
		b2.append(x1/x2)
		if(np.abs((w0 - x1)/(1 - x2)) > eps):
			b1.append((w0 - x1)/(1 - x2))
		else:
			b1.append(t)		
	b1 = np.array(b1)
	b2 = np.array(b2)
			
			
	return b1, b2
	
def GEN(F,pkex):
	out = 0
	for k in range(1,len(pkex)+1):
		out += pkex[k-1] * F**k
	return out


		
