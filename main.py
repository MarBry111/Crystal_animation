#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 15:42:35 2018

@author: marbry
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

n = 2      # number atoms on each axis
a = 1#0.38    # nm, distance between atoms
N = n ** 3  # number of all atoms
m = 1#39.948  # mass
T_0 = 1#100     # Temperature
k = 1#0.00831 # Boltzman constant
save_to_file = False


def e_kin_ax(T_0):
    lam_x = np.random.random()
    lam_y = np.random.random()
    lam_z = np.random.random()
    log_x = np.log(lam_x)
    log_y = np.log(lam_y)
    log_z = np.log(lam_z)
    e_arr = np.array([log_x, log_y, log_z])
    return -1/2*k*T_0*e_arr


def momentum_ax(T_0, m):
    temp = np.random.randint(2, size=3)
    for i in temp:
        if i == 0:
            i = -1
    return np.multiply(np.sqrt(2*m*e_kin_ax(T_0)), temp)


def force_P(epsilon, R, r_arr, index):
    r_i = r_arr[index]
    forces_P = []
    for ii in range(N-1):
        if ii != index:
            r_j = r_arr[ii]
            r_ij = np.sqrt(np.sum((r_i-r_j)**2))
            force = 12*epsilon*((R/r_ij)**12-(R/r_ij)**6)*(r_i-r_j)/r_ij**2
            force_ax = (r_i - r_j)*force
            forces_P.append(force_ax)
    return np.sum(forces_P, axis=0)


def force_S(r_i, L, f):
    r = np.sqrt(np.sum(r_i**2))
    if r < L:
        return r_i*0
    elif r >= L:
        return r_i*(f*(L - r))


# vectors showing edges of cell
b_0 = np.array([a, 0, 0])
b_1 = np.array([a / 2, a * np.sqrt(3) / 2, 0])
b_2 = np.array([a / 2, a * np.sqrt(3) / 2, a * np.sqrt(2) / np.sqrt(3)])

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

#forces
force_arr = []

L = 1
f = 1
epsilon = 1
R = 1

for j in range(N):
    force_arr.append(force_P(epsilon=epsilon, R=R, r_arr=r_arr, index=j)+force_S(r_i=r_arr[j], L=L, f=f))
    #print(force_P(epsilon=epsilon, R=R, r_arr=r_arr, index=j),force_S(r_i=r_arr[j], L=L, f=f))

force_arr = np.array(force_arr)
print(force_arr)

if save_to_file:
    delimiter = " "
    f = open("data.txt","w+")
    for j in range(N):
        f.write(str(r_i[j][0])+delimiter)
        f.write(str(r_i[j][1])+delimiter)
        f.write(str(r_i[j][2])+delimiter)
        f.write("\n")

    f.close()

#fig = plt.figure()
#ax = Axes3D(fig)
#ax.scatter3D(r_i[:, 0], r_i[:, 1], r_i[:, 2])
#ax.scatter3D(E_kin_i[:,0], E_kin_i[:,1], E_kin_i[:,2])
#ax.scatter3D(p_i[:,0], p_i[:,1], p_i[:,2])
#ax.set_ylim3d(ax.get_xlim3d())
#ax.set_zlim3d(ax.get_xlim3d())
'''
ay, ax = np.histogram(np.abs(p_arr[:,0]), 30)
by, bx = np.histogram(np.abs(p_arr[:,1]), 30)
cy, cx = np.histogram(np.abs(p_arr[:,2]), 30)
plt.plot(.5*(ax[1:]+ax[:-1]), ay)
plt.plot(.5*(bx[1:]+bx[:-1]), by)
plt.plot(.5*(cx[1:]+cx[:-1]), cy)
plt.show()
'''