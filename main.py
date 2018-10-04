# -*- coding: utf-8 -*-
"""
Created on Wed Oct  3 15:42:35 2018

@author: marbry
"""

import numpy as np
import random
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


a = 0.38 #nm, odleglosci miedzy atomami
n = 5 #ilosc atomow na krawedzi
N = n**3
m = 39.948 #masa
T_0 = 100

def E_kin (T_0):
    k = 0.00831
    lam = random.random()
    return -1/2*k*T_0*np.log(lam)
    
def momentum (E_kin, m):
    temp = random.random()
    char = 0
    if temp> 0.5: char = 1
    else: char = -1
    return char*np.sqrt(2*m*E_kin)
    

#wektory wyznaczajace krawedzi komorek
b_0 = np.array([a, 0, 0])
b_1 = np.array([a/2, a*np.sqrt(3)/2, 0])
b_2 = np.array([a/2, a*np.sqrt(3)/2, a*np.sqrt(2)/np.sqrt(3)])

#tablica z wartosciami dla i_0 i_1 i_2  (indeksy atomow)
i = range(n)
i = np.array(i)

#wspolrzedne
r_i0 = []
for j in range(n):
    for jj in range(n):
        for jjj in range(n):
            x = b_0[0]*(i[j]-(n-1)/2)+b_1[0]*(i[jj]-(n-1)/2)+b_2[0]*(i[jjj]-(n-1)/2)
            y = b_0[1]*(i[j]-(n-1)/2)+b_1[1]*(i[jj]-(n-1)/2)+b_2[1]*(i[jjj]-(n-1)/2)
            z = b_0[2]*(i[j]-(n-1)/2)+b_1[2]*(i[jj]-(n-1)/2)+b_2[2]*(i[jjj]-(n-1)/2)         
            r_i0.append([x,y,z])
r_i0 = np.array(r_i0)

#energie kinetyczne
E_kin_i = []
for j in range(n):
    for jj in range(n):
        for jjj in range(n):
            E_kin_x = E_kin(T_0)
            E_kin_y = E_kin(T_0)
            E_kin_z = E_kin(T_0)        
            E_kin_i.append([E_kin_x,E_kin_y,E_kin_z])
E_kin_i = np.array(E_kin_i)

#tu normalizacja pedu

#pedy
p_i0 = []
for j in range(N):
    p_x = momentum(E_kin_i[j][0], m)
    p_y = momentum(E_kin_i[j][1], m)
    p_z = momentum(E_kin_i[j][2], m)        
    p_i0.append([p_x,p_y,p_z])
p_i0 = np.array(p_i0)


P = np.array([np.sum(p_i0[:,0]),np.sum(p_i0[:,1]),np.sum(p_i0[:,2])])

for j in range(N):
    p_i0[j] = p_i0[j] - P

'''
delimiter = " "
f = open("data.txt","w+") 
for j in range(N):
    f.write(str(r_i0[j][0])+delimiter) 
    f.write(str(r_i0[j][1])+delimiter)
    f.write(str(r_i0[j][2])+delimiter)
    f.write("\n")
    
f.close()
'''

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
#ax.scatter3D(r_i0[:,0], r_i0[:,1], r_i0[:,2])
ax.scatter3D(E_kin_i[:,0], E_kin_i[:,1], E_kin_i[:,2])
#ax.scatter3D(p_i0[:,0], p_i0[:,1], p_i0[:,2])
plt.show()
