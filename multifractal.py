#!/usr/bin/env python

import os
import sys
import math 
import numpy as np
import pywt


data = [1, 2, 3, 4, 5, 0, 3, 2, 4, 4, 6, 2, 1, 3, 5, 2, 7, 4, 5, 8, 9, 1, 2, 3, 4, 5, 6, 2, 4, 5 ,2, 3]
filt = [0.482962913145, 0.836516303738, 0.224143868042, -0.129409522551]

# **** my implementation of dwt **** #

def decimate(c,n,nn):
	i = n
	p = []
	while i<(n+nn-1):
		p.append(c[i])
		i=i+2
	return p

def periodic(x,y,l):
	p = []
	for i in range(-(x-1),0):	
		p.append(l[divmod(i,y)[1]])
	return p

n = len(filt)
limit = int(math.log(len(data), 2))   # calculating the maximum limit of dwt levels
c = data
dwtr = []      # list containing all levels dwt coefficients
 
g = [pow(-1,(y+1))*x for y,x in enumerate(filt)]
filt.reverse()
h = filt

for j in range(0,limit):
	nn = len(c)
	m = periodic(n,nn,c)
	c = m+c
	d = np.convolve(c,g)
	d = decimate(d,n,nn)
	c = np.convolve(c,h)
	c = decimate(c,n,nn)
	dwtr.append(d)
dwtr.append(c)

#--------------------------------------------------------------#

# **** dwt coefficients from pywavelet **** #

dwtc = pywt.wavedec(data, 'db2', level = limit)

#-------------------------------------------------------------#

sumy = []   # S function calculation - square of dwt coefficients used in summation 

for i in range(0, 5):
	sumy.append(sum(map(lambda x : x*x, dwtr[i])))
	sumy[i]/=len(dwtr[i])

#-------------------------------------------------------------#

# S function calculation of higher dimensions using leader wavelets 	

def leaders(l, gamma, j):
	p = []
	p.append(l[0])
	for i in range(1, len(l)-1) :
		p.append(pow(2, (gamma*j))*max(l[i-1], l[i], l[i+1]))
	p.append(l[len(l)-1])
	return p


q = 2  # represents the power to which the leader wavelets will be raised
gamma = 1  # need to be obtained from calculation haemodynamic function
sumdy = []  
for i in range(0, 5):
	modarr = [math.fabs(x) for x in dwtr[i]]
	sumdy.append(sum(map(lambda x : pow(x, q) , leaders(modarr , gamma, i))))
	sumdy[i]/=len(dwtr[i])

#-------------------------------------------------------------#

# zeta power to be calculated in both the above cases of S value



