#!/usr/bin/python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 15:42:35 2018

@author: marbry
"""

from functions_crystal_simulation import *
import time 

n = 5      # number atoms on each axis
a = 0.38    # nm, distance between atoms
m = 39.948  # mass
T_0 = 100     # Temperature
k = 0.00831 # Boltzman constant
file_xyz = 'xyz.dat'
file_out = 'out.dat'
temp = []
temp.append(T_0)
energy_k = 0
energy_c = 0

L = 2.3 # nm
f = 10000 # nm
epsilon = 1
R = 0.38 # nm

tau = 0.01
s_d = 10000
s_0 = 1000
s_xyz = 10
s_out = 100

###2.1
with open(file_xyz, "w") as file_1:
    file_1.write("")

with open(file_out, "w") as file_2:
    file_2.write("")

N = n ** 3  # number of all atoms (3)

# vectors showing edges of cell (4)
b_0 = np.array([a, 0, 0])
b_1 = np.array([a / 2, a * np.sqrt(3) / 2, 0])
b_2 = np.array([a / 2, a * np.sqrt(3) / 6, a * np.sqrt(2) / np.sqrt(3)])

# array of values dla i_0 i_1 i_2 (atoms' indexes)
i = range(n)
i = np.array(i)

# coordinates (5)
r_arr = generate_r_arr(b0=b_0, b1=b_1, b2=b_2, n=n)

#energy (6)
e_kin_arr = generate_e_kin_arr(temp_0=T_0, N=N, k=k)
#momentum (7) (8)
p_arr = generate_momentum_arr(e_kin_arr=e_kin_arr,m=m, N=N)

#print_momentum_chart(p_arr=p_arr)
###2.2
#forces
#f_arr, f_S_arr = generate_force_arr(r_arr=r_arr, R=R, L=L, f=f, epsilon=epsilon, N=N)

f_S_arr, v_S_arr = force_potential_S(r_arr=r_arr, L=L, f=f, N=N)
#v_arr = generate_v(r_arr=r_arr, R=R, L=L, f=f, epsilon=epsilon, N=N)

f_P_arr, v_P_arr = force_potential_P(epsilon=epsilon, R=R, r_arr=r_arr, N=N)

v_arr = v_S_arr+v_P_arr
f_arr = f_S_arr+f_P_arr

#correct = -669.xxx
#print(np.sum(v_arr))

#f_arr = np.array(f_arr)
#f_S_arr = np.array(f_S_arr)

preasure_walls = np.sum(f_S_arr)/4/np.pi/(L**2)
#print(preasure_walls)

start = time.time()
###2.3
for j in range(s_d+s_0):
    p_arr_tau = generate_p_arr_tau(p_arr=p_arr, f_arr=f_arr, tau=tau)

    r_arr = generate_r_arr_t(r_arr=r_arr, p_arr_tau=p_arr_tau, m=m, tau=tau)

    #f_arr, f_S_arr = generate_force_arr(r_arr=r_arr, R=R, L=L, f=f, epsilon=epsilon, N=N)
    #v_arr = generate_v(r_arr=r_arr, R=R, L=L, f=f, epsilon=epsilon, N=N)
    
    f_S_arr, v_S_arr = force_potential_S(r_arr=r_arr, L=L, f=f, N=N)
    f_P_arr, v_P_arr = force_potential_P(epsilon=epsilon, R=R, r_arr=r_arr, N=N)
    
    v_arr = v_S_arr+v_P_arr
    f_arr = f_S_arr+f_P_arr
    
    p_arr = generate_p_arr_t(p_arr_tau=p_arr_tau, f_arr_t=f_arr, tau=tau)
    
    preasure_walls = np.sum(f_S_arr)/4/np.pi/(L**2)
    temp.append(generate_temp(p_arr=p_arr, N=N, k=k, m=m))
    energy_c, energy_k = generate_H_Ek(p_arr, m, v_arr)

    if(j==2000):
        end = time.time()
        print("2k pekło: ", end - start)
    if(j%s_xyz == 0):
        save_to_file_xyz(file_xyz=file_xyz, r_arr=r_arr, N=N)
    if(j%s_out == 0):
        save_to_file_out(file_out=file_out, t=j, v=np.sum(v_arr), E_k=energy_k, E_c=energy_c, T=temp[-1], p=preasure_walls)

