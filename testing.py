#!/usr/bin/env python

import os
import sys
import math 
import numpy as np
import pywt

def leaders(l, gamma, j):
	p = []
	p.append(l[0])
	for i in range(1, len(l)-1) :
		p.append(pow(2, (gamma*j))*max(l[i-1], l[i], l[i+1]))
	p.append(l[len(l)-1])
	return p

def cumulants(data):
	cp1 = np.average(data)
	cp2 = np.average(np.square(data)) - np.square(cp1)
	return cp1, cp2

limit = 9
f1 = open("Testdata/Yale_10.1D","r")
f5 = open("Output/zeta210.out","w")

#--------------------------------------------------------------#

f1 = f1.readlines()
f1 = f1[1:]
	
for line in f1:
	data = [float(x) for x in line.split()]
# **** dwt coefficients from pywavelet **** #

	dwtc = pywt.wavedec(data, 'db2', level = limit)

#-------------------------------------------------------------#

	sumy = []   # S function calculation - square of dwt coefficients used in summation 

	for i in range(0, limit):
		sumy.append(sum(map(lambda x : x*x, dwtc[i])))
		sumy[i]/=len(dwtc[i])

#-------------------------------------------------------------#

# S function calculation of higher dimensions using leader wavelets 	

        j1 = 3
	j2 = 6
	q = 2  # represents the power to which the leader wavelets will be raised
	gamma = 1 
	sumdy = []  
	leaderwav = []
	for i in range(0, limit):
		modarr = [math.fabs(x) for x in dwtc[i]]
		leaderwav.append(leaders(modarr, gamma, i))
		sumdy.append(sum(map(lambda x : pow(x, q) , leaderwav[i])))
		sumdy[i]/=len(dwtc[i])

#-------------------------------------------------------------#

# weights calculation
	v = np.zeros(3)

	for i in range(0, 3):
		for k in range(j1, j2+1):
			v[i] += pow(k, i) * len(dwtc[k])

	w = np.zeros(j2+1)
	for i in range(j1, j2+1):
		w[i] = len(dwtc[i]) * (v[0]*i-v[1])/(v[0]*v[2] - v[1]*v[1])

#zeta power calculation 
	zeta = 0
	logsum = np.log2(sumdy)
	for i in range(j1, j2):
		zeta += w[i]*logsum[i]
	f5.write(str(zeta))

#-------------------------------------------------------------#

#c1 and c2 coefficients calculation
	

	c1 = 0
	c2 = 0

	for i in range(j1, j2+1):
		cp1,cp2 = cumulants(np.log(leaderwav[i]))
		c1 += w[i]*cp1
		c2 += w[i]*cp2

	c1 = np.log2(np.e) * c1 - gamma
	c2 = np.log2(np.e) * c2


	f5.write("\n")
#-------------------------------------------------------------#

#D(h) calculation and plotting
# D(h) = 1 - (h-c1)^2/c2

#-------------------------------------------------------------#





