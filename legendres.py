#!/home/coffee/anaconda3/bin/python

import numpy as np
from scipy.special import lpmn
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as axes3d

njs = 4; 
nms = njs;

coeffs = np.zeros((nms,njs),dtype=float);
l,m=(0,0);
coeffs[m,l] = .5;
l,m=(1,0);
coeffs[m,l] = .15;
l,m=(1,1);
coeffs[m,l] = .1
l,m=(2,0);
coeffs[m,l] = .0225;
l,m=(2,1);
coeffs[m,l] = .025;
l,m=(2,2);
coeffs[m,l] = .2;
l,m=(3,0);
coeffs[m,l] = .03125;
l,m=(3,1);
coeffs[m,l] = .0325;
l,m=(3,2);
coeffs[m,l] = .035;
l,m=(3,3);
coeffs[m,l] = .3;

phi, theta = np.linspace(0, 2 * np.pi, 20), np.linspace(0, np.pi, 10)
PHI,THETA = np.meshgrid(phi,theta)
R = np.zeros(THETA.shape,dtype=float);
print(THETA.shape[0])
for i in range(THETA.shape[0]): # for ever phi sample point
	c = np.cos(THETA[i,0]);
	#print(THETA[i,0],c);
	Plm = lpmn(njs-1,njs-1,c)[0];
	r = 0.;
	for m in range(nms):
		for j in range(njs):
			r += coeffs[m,j]*Plm[m,j];
	R[i,:].fill(r);	
	if i == THETA.shape[0]//2:
		print('Plm',Plm);
print('coeffs',coeffs);
print(R);

X = R * np.sin(THETA) * np.cos(PHI)
Y = R * np.sin(THETA) * np.sin(PHI)
Z = R * np.cos(THETA)
fig = plt.figure()
ax = fig.add_subplot(1,1,1, projection='3d')
plot = ax.plot_surface(
    X, Y, Z, rstride=1, cstride=1, cmap=plt.get_cmap('jet'),
        linewidth=0, antialiased=False, alpha=0.5)

plt.show()
