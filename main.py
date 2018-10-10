#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 15:42:35 2018

@author: marbry
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

n = 5      # number atoms on each axis
a = 0.38    # nm, distance between atoms
N = n ** 3  # number of all atoms
m = 39.948  # mass
T_0 = 100     # Temperature
k = 0.00831 # Boltzman constant
save_to_file = False
print_3D_r = True
print_momentum_chart = False


L = 2.3 # nm
f = 10000 # nm
epsilon = 1
R = 0.38 # nm

###2.1
def e_kin_ax(T_0):
    log_x = np.log(np.random.random())
    log_y = np.log(np.random.random())
    log_z = np.log(np.random.random())
    e_arr = np.array([log_x, log_y, log_z])
    return -1/2*k*T_0*e_arr


def momentum_ax(T_0, m):
    temp = np.random.randint(2, size=3)
    for j in range(3):
		if temp[j] == 0:
			temp[j] = -1
    return np.multiply(np.sqrt(2*m*e_kin_ax(T_0)), temp)

###2.2
def force_P(epsilon, R, r_arr, index):
    r_i = r_arr[index]
    forces_P = []
    for ii in range(N-1):
        if ii != index:
            r_j = r_arr[ii]
            r_ij = np.sqrt(np.sum((r_i-r_j)**2))
            force = 12*epsilon*((R/r_ij)**12-(R/r_ij)**6)/(r_ij**2)
            force_ax = (r_i - r_j)*force
            forces_P.append(force_ax)
    return np.sum(forces_P, axis=0)

def force_S(r_i, L, f):
    r = np.sqrt(np.sum(r_i**2))
    if r < L:
        return r_i*0
    elif r >= L:
        return r_i*(f*(L - r))

def potential_S(r_i, L, f):
	r = np.sqrt(np.sum(r_i**2))
	if r < L:
		return r*0
	elif r >= L:
		return f/2*(r - L)**2

def potential_P(epsilon, R, r_arr, index, N):
	r_i = r_arr[index]
	v_P = []
	for ii in range(N):
		if ii != index:
			r_j = r_arr[ii]
			r_ij = np.sqrt(np.sum((r_i-r_j)**2))
			v = epsilon*((R/r_ij)**12-2*(R/r_ij)**6)
			v_P.append(v)
	return v_P

###2.3
#equation of move
def r_t(r, p, m, tau, t):
	r_tau = []
	r_t = r
	
###2.1
# vectors showing edges of cell
b_0 = np.array([a, 0, 0])
b_1 = np.array([a / 2, a * np.sqrt(3) / 2, 0])
b_2 = np.array([a / 2, a * np.sqrt(3) / 6, a * np.sqrt(2) / np.sqrt(3)])

# array of values dla i_0 i_1 i_2 (atoms' indexes)
i = range(n)
i = np.array(i)

# coordinates
r_arr = []
for j in range(n):
    for jj in range(n):
        for jjj in range(n):
            x = b_0[0] * (i[j] - (n - 1) / 2) + b_1[0] * (i[jj] - (n - 1) / 2) + b_2[0] * (i[jjj] - (n - 1) / 2)
            y = b_0[1] * (i[j] - (n - 1) / 2) + b_1[1] * (i[jj] - (n - 1) / 2) + b_2[1] * (i[jjj] - (n - 1) / 2)
            z = b_0[2] * (i[j] - (n - 1) / 2) + b_1[2] * (i[jj] - (n - 1) / 2) + b_2[2] * (i[jjj] - (n - 1) / 2)
            r_arr.append([x, y, z])
r_arr = np.array(r_arr)

# momentum
p_arr = []
for j in range(N):
    p_arr.append(momentum_ax(T_0=T_0, m=m))
p_arr = np.array(p_arr)

# normalize momentum
P = np.array([np.sum(p_arr[:, 0]), np.sum(p_arr[:, 1]), np.sum(p_arr[:, 2])])/N
p_arr[:] = p_arr[:] - P

###2.2
#forces
force_arr = []
force_S_arr = []
for j in range(N):
	f_S = force_S(r_i=r_arr[j], L=L, f=f)
	f_P = force_P(epsilon=epsilon, R=R, r_arr=r_arr, index=j)
	force_arr.append(f_P+f_S)
	force_S_arr.append(f_S)
force_arr = np.array(force_arr)
force_S_arr = np.array(force_S_arr)
preasure_walls = np.sum(force_S_arr)/4/np.pi/(L**2)

# potential (working good - time to improve)
v_arr = []
for j in range(N):
	v_P = potential_P(epsilon=epsilon, R=R, r_arr=r_arr, index=j, N=N)
	v_S = potential_S(r_i=r_arr[j], L=L, f=f)
	v_P_sum = 0
	for jj in range(j):
		v_P_sum = v_P_sum + v_P[jj]
	v_arr.append(v_P_sum+v_S)
#correct = -669.xxx

###2.3









### plots save to file thios kinde of stuff

if save_to_file:
    delimiter = " "
    f = open("data.txt","w+")
    for j in range(N):
        f.write(str(r_i[j][0])+delimiter+str(r_i[j][1])+delimiter+str(r_i[j][2]))
        f.write("\n")

    f.close()

if print_3D_r:
	fig = plt.figure()
	ax = Axes3D(fig)
	#ax.scatter3D(r_arr[:, 0], r_arr[:, 1], r_arr[:, 2])
	ax.scatter3D(force_arr[:,0], force_arr[:,1], force_arr[:,2])
	ax.set_ylim3d(ax.get_xlim3d())
	ax.set_zlim3d(ax.get_xlim3d())
	plt.show()

if print_momentum_chart:
	ay, ax = np.histogram(np.abs(p_arr[:,0]), 30)
	by, bx = np.histogram(np.abs(p_arr[:,1]), 30)
	cy, cx = np.histogram(np.abs(p_arr[:,2]), 30)
	plt.plot(.5*(ax[1:]+ax[:-1]), ay)
	plt.plot(.5*(bx[1:]+bx[:-1]), by)
	plt.plot(.5*(cx[1:]+cx[:-1]), cy)
	plt.show()

