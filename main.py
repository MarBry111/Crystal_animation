#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 15:42:35 2018

@author: marbry
"""

import numpy as np
import matplotlib.pyplot as plt
import itertools
from mpl_toolkits.mplot3d import Axes3D

n = 5      # number atoms on each axis
a = 0.38    # nm, distance between atoms
N = n ** 3  # number of all atoms
m = 39.948  # mass
T_0 = 100     # Temperature
k = 0.00831 # Boltzman constant
save_to_file = False
print_3D_r = False
print_momentum_chart = False
file_xyz = 'xyz.dat'
file_out = 'out.dat'

'''
# xyz.dat -> 
x y z
\n
\n 
x y z

# out.dat
t V E_k E_c T p
'''
L = 2.3 # nm
f = 10000 # nm
epsilon = 1
R = 0.38 # nm

tau = 0.01
s_d = 100000
s_0 = 1000
s_xyz = 10
s_out = 100


###2.1
def generate_r_arr(b0, b1, b2, n):
    i_arr = np.arange(n)
    r_arr = np.array(list(p for p in itertools.product(i_arr, repeat=3)), dtype=float)
    for j in range(3):
        r_arr[:, j] = (r_arr[:, 0] - (n - 1) / 2) * b0[j] + (r_arr[:, 1] - (n - 1) / 2) * b1[j] + (r_arr[:, 2] - (n - 1) / 2) * b2[j]
    return r_arr


def generate_e_kin_arr(temp_0, N):
    rand_arr = np.random.rand(N, 3)
    e_arr = np.log(rand_arr)
    e_arr = -1/2*k*temp_0*e_arr
    return e_arr


def generate_momentum_arr(e_kin_arr, m, N):
    temp = np.random.randint(2, size=(N, 3))
    temp[temp == 0] = -1
    p_arr = np.sqrt(2*m*e_kin_arr)
    #normalize momentum
    P = np.array([np.sum(p_arr[:, 0]), np.sum(p_arr[:, 1]), np.sum(p_arr[:, 2])]) / N
    p_arr = p_arr - P
    return p_arr

###2.2
###TO DOOO
def force_P(epsilon, R, r_arr, index):
    r_i = r_arr[index]
    forces_P = []
    for ii in range(N):
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
        r_out = r_i
        try:
            r_out = r_i*(f*(L - r))
        except TypeError:
            print(r_i,f,L,r)
        return r_out



def generate_force_arr(r_arr, R, L, f, epsilon, N):
    f_arr = []
    f_S_arr = []
    for jj in range(N):
        f_S = force_S(r_i=r_arr[jj], L=L, f=f)
        f_P = force_P(epsilon=epsilon, R=R, r_arr=r_arr, index=jj)
        f_arr.append(f_P+f_S)
        f_S_arr.append(f_S)
    f_arr = np.array(f_arr)
    f_S_arr = np.array(f_S_arr)
    return f_arr, f_S_arr

###DONE
def potential_S(r_arr, L, f):
    r = np.sqrt(np.sum(r_arr**2,axis=1))
    v_s = f/2*(r - L)**2
    v_s[v_s<L] =0
    return v_s

###TO DO
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

def generate_v(r_arr, R, L, f, epsilon, N):
    v_arr = []
    for j in range(N):
        v_P = potential_P(epsilon=epsilon, R=R, r_arr=r_arr, index=j, N=j)
        v_S = potential_S(r_i=r_arr[j], L=L, f=f)
        v_P_sum = 0
        for jj in range(j):
            v_P_sum = v_P_sum + v_P[jj]
        v_arr.append(v_P_sum+v_S)
    v_arr = np.array(v_arr)
    return v_arr


###2.3
#p (t i 1/2 tau)
def generate_p_arr_tau(p_arr, f_arr, tau):
    p_arr_tau = p_arr + f_arr*tau/2
    return p_arr_tau


def generate_r_arr_t(r_arr, p_arr_tau, m, tau):
    r_arr_t = r_arr + p_arr_tau*tau/m
    return r_arr_t


def generate_p_arr_t(p_arr_tau, f_arr_t, tau):
    p_arr_t = p_arr_tau + f_arr_t*tau/2
    return p_arr_t


def save_to_file_xyz(file_xyz, r_arr):
    delimiter = " "
    with open(file_xyz,"a") as file_1:
        for j in range(N):
            file_1.write(str(r_arr[j][0])+delimiter+str(r_arr[j][1])+delimiter+str(r_arr[j][2])+"\n")
        file_1.write("\n \n")


def save_to_file_out(file_out, t, v, E_k, E_c, T, p):
    delimiter = " "
    with open(file_out,"a") as file_1:
        file_1.write(str(t)+delimiter+str(v)+delimiter+str(E_k)+delimiter+str(E_c)+delimiter+str(T)+delimiter+str(E_c)+"\n \n")




###2.1
with open(file_xyz,"w") as file_1:
    file_1.write("")

with open(file_out,"w") as file_2:
    file_2.write("")

# vectors showing edges of cell
b_0 = np.array([a, 0, 0])
b_1 = np.array([a / 2, a * np.sqrt(3) / 2, 0])
b_2 = np.array([a / 2, a * np.sqrt(3) / 6, a * np.sqrt(2) / np.sqrt(3)])

# array of values dla i_0 i_1 i_2 (atoms' indexes)
i = range(n)
i = np.array(i)

# coordinates
r_arr = generate_r_arr(b0=b_0, b1=b_1, b2=b_2, n=n)

#energy and momentum
e_kin_arr = generate_e_kin_arr(temp_0=T_0, N=N)
p_arr = generate_momentum_arr(e_kin_arr=e_kin_arr,m=m, N=N)

###2.2
#forces
f_arr, f_S_arr = generate_force_arr(r_arr=r_arr, R=R, L=L, f=f, epsilon=epsilon, N=N)

v_arr = generate_v(r_arr=r_arr, R=R, L=L, f=f, epsilon=epsilon, N=N)

'''
for j in range(N):
    f_S = force_S(r_i=r_arr[j], L=L, f=f)
    f_P = force_P(epsilon=epsilon, R=R, r_arr=r_arr, index=j)
    f_arr.append(f_P+f_S)
    f_S_arr.append(f_S)
    
    v_P = potential_P(epsilon=epsilon, R=R, r_arr=r_arr, index=j, N=N)
    v_S = potential_S(r_i=r_arr[j], L=L, f=f)
    v_P_sum = 0
    for jj in range(j):
        v_P_sum = v_P_sum + v_P[jj]
    v_arr.append(v_P_sum+v_S)
'''
#correct = -669.xxx	
#print(np.sum(v_arr))

#f_arr = np.array(f_arr)
#f_S_arr = np.array(f_S_arr)

preasure_walls = np.sum(f_S_arr)/4/np.pi/(L**2)

###2.3
for j in range(s_d+s_0):
    p_arr_tau = generate_p_arr_tau(p_arr=p_arr, f_arr=f_arr, tau=tau)

    r_arr = generate_r_arr_t(r_arr=r_arr, p_arr_tau=p_arr_tau, m=m, tau=tau)
    '''
    f_arr_t = []
    for jj in range(N):
        f_S = force_S(r_i=r_arr[jj], L=L, f=f)
        f_P = force_P(epsilon=epsilon, R=R, r_arr=r_arr, index=jj)
        f_arr_t.append(f_P+f_S)
    f_arr = np.array(f_arr_t)
    '''
    f_arr, f_S_arr = generate_force_arr(r_arr=r_arr, R=R, L=L, f=f, epsilon=epsilon, N=N)

    p_arr = generate_p_arr_t(p_arr_tau=p_arr_tau, f_arr_t=f_arr, tau=tau)
    if(j%s_xyz == 0):
        save_to_file_xyz(file_xyz=file_xyz, r_arr=r_arr)
    if(j%s_out == 0):
        save_to_file_out(file_out=file_out, t=j, v=np.sum(v_arr), E_k=1, E_c=1, T=T_0, p=1)





### plots save to file thios kinde of stuff

if print_3D_r:
    fig = plt.figure()
    ax = Axes3D(fig)
    #ax.scatter3D(r_arr[:, 0], r_arr[:, 1], r_arr[:, 2])
    ax.scatter3D(force_arr[:,0], force_arr[:,1], force_arr[:,2])
    ax.set_ylim3d(ax.get_xlim3d())
    ax.set_zlim3d(ax.get_xlim3d())
    plt.show()

if print_momentum_chart:
    ay, ax = np.histogram(p_arr[:,0], 30)
    by, bx = np.histogram(p_arr[:,1], 30)
    cy, cx = np.histogram(p_arr[:,2], 30)
    plt.plot(.5*(ax[1:]+ax[:-1]), ay)
    plt.plot(.5*(bx[1:]+bx[:-1]), by)
    plt.plot(.5*(cx[1:]+cx[:-1]), cy)
    plt.show()

