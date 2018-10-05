# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 15:42:35 2018

@author: marbry
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

n = 20      # number atoms on each axis
a = 0.38    # nm, distance between atoms
N = n ** 3  # number of all atoms
m = 39.948  # mass
T_0 = 100     # Temperature
k = 0.00831 # Boltzman constant
save_to_file = False

def e_kin_ax(T_0):
    lam = np.random.random()
    return -1 / 2 * k * T_0 * np.log(lam)

'''
def momentum_ax(E_kin, m):
    temp = random.random()
    char = 0
    if temp > 0.5:
        char = 1
    else:
        char = -1
    return char * np.sqrt(2 * m * E_kin)
'''


def momentum_ax(T_0, m):
    temp = np.random.randint(0, 2)
    char = 0
    if temp == 1:
        char = 1
    else:
        char = -1
    return char * np.sqrt(2 * m * e_kin_ax(T_0))


# vectors showing edges of cell
b_0 = np.array([a, 0, 0])
b_1 = np.array([a / 2, a * np.sqrt(3) / 2, 0])
b_2 = np.array([a / 2, a * np.sqrt(3) / 2, a * np.sqrt(2) / np.sqrt(3)])

# array of values dla i_0 i_1 i_2 (atoms' indexes)
i = range(n)
i = np.array(i)

# coordinates
r_i = []
for j in range(n):
    for jj in range(n):
        for jjj in range(n):
            x = b_0[0] * (i[j] - (n - 1) / 2) + b_1[0] * (i[jj] - (n - 1) / 2) + b_2[0] * (i[jjj] - (n - 1) / 2)
            y = b_0[1] * (i[j] - (n - 1) / 2) + b_1[1] * (i[jj] - (n - 1) / 2) + b_2[1] * (i[jjj] - (n - 1) / 2)
            z = b_0[2] * (i[j] - (n - 1) / 2) + b_1[2] * (i[jj] - (n - 1) / 2) + b_2[2] * (i[jjj] - (n - 1) / 2)
            r_i.append([x, y, z])
r_i = np.array(r_i)

'''
# energie kinetyczne
E_kin_i = []
for j in range(n):
    for jj in range(n):
        for jjj in range(n):
            E_kin_x = e_kin_ax(T_0)
            E_kin_y = e_kin_ax(T_0)
            E_kin_z = e_kin_ax(T_0)
            E_kin_i.append([E_kin_x, E_kin_y, E_kin_z])
E_kin_i = np.array(E_kin_i)
'''

# momentum
p_i = []
for j in range(N):
    p_x = momentum_ax(T_0, m)
    p_y = momentum_ax(T_0, m)
    p_z = momentum_ax(T_0, m)
    p_i.append([p_x, p_y, p_z])
p_i = np.array(p_i)

# normalize momentum
P = np.array([np.sum(p_i[:, 0]), np.sum(p_i[:, 1]), np.sum(p_i[:, 2])])/N
p_i[:] = p_i[:] - P

#dopisac sily

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
ay, ax = np.histogram(np.abs(p_i[:,0]), 30)
by, bx = np.histogram(np.abs(p_i[:,1]), 30)
cy, cx = np.histogram(np.abs(p_i[:,2]), 30)
plt.plot(.5*(ax[1:]+ax[:-1]), ay)
plt.plot(.5*(bx[1:]+bx[:-1]), by)
plt.plot(.5*(cx[1:]+cx[:-1]), cy)
plt.show()
