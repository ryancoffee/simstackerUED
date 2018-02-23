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
parser.add_argument('-f', '--file', dest='file', type=str, default='alignment.dat', help='filename for the cossq data (alignment.dat)');
parser.add_argument('-n', '--nbases', dest='nbases', type=int, default=6, help='number of basis functions (6)');
parser.add_argument('-w', '--rollwin', dest='rollwin', type=int, default=40, help='window size for rolloff ends (40)');
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

def straighten(r,k):
	print('r[-1,k]=', r[-1,k], '\tr[0,k]=', r[0,k]);
	m = (r[-1,k] - r[0,k])/(r.shape[0]-1);
	print('slope = ',m);
	r0=r[0,k];
	x = np.arange(r.shape[0]);
	x.shape=(r.shape[0],1);
	r[:,k] -= m*x[:,0]+r0;
	return r;

def removemean(r,k):
	m = np.mean(r[:,k]);
	r[:,k] -= m;
	return r;

def rollends(r,k,w):
	rollvec=np.power(np.sin(np.linspace(0,PI/2.,w,dtype=float)),2);
	rollvec.shape=(len(rollvec),1);
	r[:rollvec.shape[0],k] *= rollvec[:,0];
	r[-rollvec.shape[0]:,k] *= np.flipud(rollvec[:,0]);
	return r;
	
def zeroends(r,k,w):
	r=straighten(r,k);
	offset = r[0,k];
	r[:,k] -= r[0,k];
	r=rollends(r,k,w);
	print('zeroends:\t',r[:10,k],'...',r[-10:,k]);
	return r;

def sinsqfilter(r,k):
	v=np.power(np.cos(np.linspace(0,PI/2.,k,dtype=float)),2);
	print(r.shape,'\t<-- shape(r) : k -->\t',k);
	r[:len(v),k] *= v;
	r[-len(v):,k] *= np.flipud(v);
	r[len(v)+1:-len(v)-1,k] *= 0.;
	return r;
	

def normalize(r,k):
	c = float(np.real(np.sum(np.conj(r[:,k])*r[:,k])));
	r[:,k] /= c;
	print('normalize:\t',r[:10,k]);
	return r;

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


def orthonormalize(r,k,w):
	# first straighten, then zeroends then removemean, then remove the dots with previous (grahm), then normalize
	r = straighten(r,k);
	#r = zeroends(r,k,w);
	r = removemean(r,k);
	if k==0:
		return r;
	c = np.zeros((1,k),dtype=float);
	for i in range(k):
		c[0,i] = np.real(np.dot(r[:,k],r[:,i]));
	for i in range(k):
		r[:,k] -= np.real(c[0,i] * r[:,i]);
	r = normalize(r,k);
	return r;


### When constructing wavelets, the sum(psi(x)) = 0 and the normalization sum{ || psi(x) || } = 1

#########	Here we are actually constructing an orthonormal basis with ground state of the <cos**2 th> 	#########
#########	We do this by fourier transforming the mother, shifting up in frequency by one, then GS orthogonalizing	#######



tvec = np.loadtxt(args.file,dtype=float,usecols=(0,));
cossq = np.loadtxt(args.file,dtype=float,usecols=(1,));

rotor = np.zeros((len(cossq),args.nbases),dtype=complex);
rotorFT = np.zeros(rotor.shape,dtype=complex);


####	now fill accordingly	####
rotor[:,0] = cossq;
print(rotor[:10,:]);
rotor = orthonormalize(rotor,0,args.rollwin);
rotorFT[:,0] = FFT(rotor[:,0]); 
for k in range(1,rotor.shape[1]):
	# set new colFT to be prevcolFT shifted by one frequency unit #
	rotorFT[1:rotorFT.shape[0]//2,k] = rotorFT[:rotorFT.shape[0]//2-1,k-1];
	rotorFT[rotorFT.shape[0]//2:-2,k] = rotorFT[rotorFT.shape[0]//2+1:-1,k-1];
	# IFFT and orthonormalize #
	rotor[:,k] = np.real(IFFT(rotorFT[:,k]));
	rotor = orthonormalize(rotor,k,args.rollwin);
	rotorFT[:,k] = FFT(rotor[:,k]);


icmvec = FFTFREQS(len(tvec),tvec[1]-tvec[0]) / C_cmPps;
icmvec.shape = (rotor.shape[0],1);
tvec.shape = (rotor.shape[0],1);

outname = args.file + '.real';
np.savetxt(outname,np.concatenate((tvec,np.real(rotor)),axis=1),fmt='%.6e');
outname = args.file + '.imag';
np.savetxt(outname,np.concatenate((tvec,np.imag(rotor)),axis=1),fmt='%.6e');
outname = args.file + '.fft.abs';
np.savetxt(outname,np.concatenate((icmvec,np.abs(rotorFT)),axis=1),fmt='%.6e');
outname = args.file + '.fft.arg';
np.savetxt(outname,np.concatenate((icmvec,np.angle(rotorFT)),axis=1),fmt='%.6e');

############	Now try constructing the wavelet functions and FTs 	###########
############	Here we will use again <cos**2 th> as the mother	###########


wavelet = np.zeros((len(cossq),args.nbases),dtype=complex);
waveletFT = np.zeros(wavelet.shape,dtype=complex);
wavelet[:,0] = cossq;
wavelet = orthonormalize(wavelet,0,args.rollwin);
for j in range(1,args.nbases):
	wavelet[:wavelet.shape[0]//int(np.power(2,j)),j] = 1./np.sqrt(np.power(2,j))*wavelet[::int(np.power(2,j)),0];

waveletFT = FFT(wavelet,axis=1);
outname = args.file + '.wave';
np.savetxt(outname,np.concatenate((tvec,np.real(wavelet)),axis=1),fmt='%.6e');
outname = args.file + '.wave.abs';
np.savetxt(outname,np.concatenate((icmvec,abs(waveletFT)),axis=1),fmt='%.6e');
outname = args.file + '.wave.arg';
np.savetxt(outname,np.concatenate((icmvec,np.angle(waveletFT)),axis=1),fmt='%.6e');

##########	Use pywavelets to construct your own 	########
#[dec_lo dec_hi rec_lo rec_hi]
#dec_lo is the low frequency e.g. <cos**2 th> and the dec_hi is hte high frequency channel, or the derivative, (or the second deriv, or the freq shifted
