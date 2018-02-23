#!/home/coffee/anaconda3/bin/python
import numpy as np;
import argparse
import matplotlib.pyplot as plt;
import scipy
from scipy import special
import cmath as cm
from constants import EPS,PI,C_cmPps
from numpy.fft import fft as FFT
from numpy.fft import ifft as IFFT
from numpy.fft import fftfreq as FFTFREQS
import numpy.random as rng

parser = argparse.ArgumentParser(description='FFT <cos^2 theta > file');
parser.add_argument('-f', '--file', dest='file', type=str, default='', help='filename for the cossq data');
args = parser.parse_args();

def gauss(x,x0,w):
	return np.exp(-((x-x0)/w)**2);

def i2icm(x,twin):
	return x/(twin*C_cmPps);

def smooth(x,k):
	for i in range(len(x)):
		val = 0.;
		for j in range(k):
			val += x[(i+j)%len(x)];
		x[i] = val/k;
	return x;

def straighten(z):
	print('z[-1]=', z[-1], '\tz[0]=', z[0]);
	m = (z[-1] - z[0])/(len(z)-1);
	print('slope = ',m);
	z0=z[0];
	x = np.arange(len(z));
	z -= m*x+z0;
	return z;

def removemean(y):
	m = np.mean(y);
	y -= m;
	return y;

def rollends(y,k):
	rollvec=np.power(np.sin(np.linspace(0,PI/2.,k,dtype=float)),2);
	y[:len(rollvec)] *= rollvec;
	y[-len(rollvec):] *= np.flipud(rollvec);
	return y;
	
def zeroends(y,k):
	y=straighten(y);
	offset = y[0];
	y -= y[0];
	y=rollends(y,k);
	print('zeroends:\t',y[:10],'...',y[-10:]);
	return y;

def sinsqfilter(z,k):
	v=np.power(np.cos(np.linspace(0,PI/2.,k,dtype=float)),2);
	print(len(v),'\t<-- len(v) : k -->\t',k);
	z[:len(v)] *= v;
	z[-len(v):] *= np.flipud(v);
	z[len(v)+1:-len(v)-1] *= 0.;
	return z;
	

def normalize(z):
	c = float(np.real(np.dot(np.conj(z),z)));
	z /= np.sqrt(c);
	print('normalize:\t',z[:10]);
	return z;

def weinerfilter(vecFT):
	## so the program should be to weiner filter properly, using a fitting of the logarithm, it will be a linear fit
	lps = np.log(np.real(np.conj(vecFT)*vecFT) + EPS); # working in log space for y
	f=FFTFREQS(len(lps));
	msig1,bsig1=np.polyfit( f[len(f)//32:len(f)//8] , lps[len(f)//32:len(f)//8] , 1);
	msig2,bsig2=np.polyfit( f[len(f)//8:len(f)//4] , lps[len(f)//8:len(f)//4] ,1);
	mnoise,bnoise=np.polyfit( f[len(f)//4:len(f)//2] , lps[len(f)//4:len(f)//2] ,1);
	s = np.exp(bsig1) * np.exp(msig1*abs(f));
	s += np.exp(bsig2) * np.exp(msig2*abs(f));
	n = np.exp(bnoise) * np.exp(mnoise*abs(f));
	return s/(s+n);

### When constructing wavelets, the sum(psi(x)) = 0 and the normalization sum{ || psi(x) || } = 1

tvec = np.loadtxt(args.file,dtype=float,usecols=(0,));
cossq = np.loadtxt(args.file,dtype=float,usecols=(1,));
freqs = FFTFREQS(len(tvec));
#cossq = np.arange(len(tvec),dtype=float);
#print(cossq);
cossq = straighten(cossq);
cossq = zeroends(cossq,40);
cossq = removemean( cossq );
#cossq = removemean( cossq );
#cossq = zeroends(cossq,20);
#cossq = removemean( cossq );
#cossqFT = FFT(cossq); 
#wf = weinerfilter(cossqFT);
#cossqFT *= wf;
#cossq = IFFT(cossqFT);
#cossq = zeroends(cossq,40);
cossq = normalize(cossq);
print(np.dot(cossq,cossq), '<< dot : len >> ' , len(cossq));

cossqFT = FFT(cossq); 
ddcossqFT=np.zeros(cossqFT.shape,dtype=complex);
i=1;
ddcossqFT[i:len(ddcossqFT)//2] = cossqFT[:len(cossqFT)//2-i];
ddcossqFT[len(ddcossqFT)//2:-(i+1)] = cossqFT[len(cossqFT)//2+i:-1];
dcossqsq = np.real(IFFT(ddcossqFT));

cossqsq = np.power(cossq,2);
#cossqsq = zeroends(cossqsq,40);
# OK, here we FFT, then compute Weiner power spectrum, then multiply by (1j*freqs)^2 for second derivative, then multiply by Weiner then back transform to time.
cossqsq = straighten(cossqsq);
deriv = np.diff(cossqsq);
deriv = np.append(deriv,0.);
dderiv = np.diff(deriv);
dderiv = np.append(dderiv,0.);
dcossqsqFT = FFT(cossqsq);
dcossqsqFT *= -1.*(np.power(freqs,2));
##wf = weinerfilter(dcossqsqFT);
##dcossqsqFT *= wf;
#dcossqsq = np.real(IFFT(dcossqsqFT)); # careful here, the imaginary part should be inspected when doing this.
#dcossqsq = straighten(dcossqsq);
dcossqsq = straighten(dcossqsq);
dcossqsq = zeroends(dcossqsq,40);
dcossqsq = removemean(dcossqsq);

#eventually iterate over this to get closer to 0 mean
coeff = np.real(np.dot(dcossqsq,cossq));
print(np.dot(dcossqsq,cossq));
dcossqsq -= np.real(coeff * cossq);
dcossqsq = normalize(dcossqsq);


cos4sq = np.conj(dcossqsq)*dcossqsq;
#cos4sq = zeroends(cos4sq,40)
dcos4sqFT = FFT(cos4sq);
dcos4sqFT *= -1.*(np.power(freqs,2));
#wf = weinerfilter(dcos4sqFT);
#dcos4sqFT *= wf;
dcos4sq = np.real(IFFT(dcos4sqFT));
dcos4sq = normalize(dcos4sq);
coeff0 = np.real(np.dot(dcos4sq,cossq));
coeff1 = np.real(np.dot(dcos4sq,dcossqsq));
dcos4sq -= np.real(coeff0 * cossq);
dcos4sq -= np.real(coeff1 * cossqsq);



tvec.shape = (len(tvec),1);
cossq.shape = (len(cossq),1);
cossqsq.shape = (len(cossqsq),1);
dcossqsq.shape = (len(dcossqsq),1);
dcos4sq.shape = (len(dcos4sq),1);



#zFFT = FFT(np.roll(cossq,int(cossq.shape[0]/2))); # this moves a presumed central feature to the edges of the vector so the phase will have small slope
zFFT = FFT(cossq,axis=0); 
zFFT2 = FFT(dcossqsq,axis=0);
zFFT4 = FFT(dcos4sq,axis=0);
icmvec = FFTFREQS(len(tvec),tvec[1]-tvec[0]) / C_cmPps;
icmvec.shape = zFFT.shape;
outname = args.file + '.fft';
out = np.concatenate( (icmvec , abs(zFFT) , np.angle(zFFT), abs(zFFT2), np.angle(zFFT2), abs(zFFT4), np.angle(zFFT4)) , axis=1);
print( out.shape);
np.savetxt(outname,out,fmt='%.6e');
outname = args.file + '.real';
outreal = np.concatenate( (tvec , np.real(cossq) , np.real(cossqsq), np.real(dcossqsq), np.real(dcos4sq)) , axis=1);
np.savetxt(outname,outreal,fmt='%.6e');
outimag = np.concatenate( (tvec , np.imag(cossq) , np.imag(cossqsq), np.imag(dcossqsq), np.imag(dcos4sq)) , axis=1);
outname = args.file + '.imag';
np.savetxt(outname,outimag,fmt='%.6e');
