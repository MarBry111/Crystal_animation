import numpy as np
import matplotlib.pyplot as plt
import itertools
from mpl_toolkits.mplot3d import Axes3D

'''
# xyz.dat -> 
x y z
\n
\n 
x y z

# out.dat
t V E_k E_c T p
'''

###2.1
def generate_r_arr(b0, b1, b2, n):
    i_arr = np.arange(n)
    r_arr = np.array(list(p for p in itertools.product(i_arr, repeat=3)), dtype=float)
    for j in range(3):
        r_arr[:, j] = (r_arr[:, 0] - (n - 1) / 2) * b0[j] + (r_arr[:, 1] - (n - 1) / 2) * b1[j] + (r_arr[:, 2] - (n - 1) / 2) * b2[j]
    return r_arr


def generate_e_kin_arr(temp_0, N, k):
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
## TO DO
'''
def force_P(epsilon, R, r_arr, index, N):
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
'''

def force_potential_P(epsilon, R, r_arr, N):
	f_P = np.empty(shape=(N,3))
	v_P = np.zeros(shape=(N,1))
	for j in range(N):
		#ri - rj (array)
		r_i = r_arr[j].reshape(1,3) 
		# sprawdzic dokadnie wiersz a nie czy zawierta
		r_j = r_arr[:,(r_arr[:]!=r_i[:])]
		print(r_j.shape)
		r_j = r_j.reshape(N-1, 3)
		ri_j = (r_i-(r_arr!=r_i))
		#skalar r_ij
		r_ij = np.sqrt(np.sum((ri_j**2), axis=1))
		r_ij = r_ij.reshape(N,1)
		f_P[j] = np.sum((12*epsilon*((R/r_ij)**12-(R/r_ij)**6)/(r_ij**2)*ri_j), axis=0)
		r_ij = r_ij.reshape(1,N)
		v = epsilon*((R/r_ij)**12-2*(R/r_ij)**6)
		v_P[j] = np.sum(epsilon*((R/r_ij[:j])**12-2*(R/r_ij[:j])**6))
	return f_P, v_P
		


def force_potential_S(r_arr, L, f):
    r = np.sqrt(np.sum(r_arr**2, axis=1))
    f_s = r_arr*f*(L-r[:, None])
    f_s[r < L] = 0
    v_s = f/2*(r - L)**2
    v_s[r < L] = 0
    return f_s, v_s


## TO DO
'''
def generate_force_arr(r_arr, R, L, f, epsilon, N):
    f_arr = []
    f_S_arr = force_S(r_arr=r_arr, L=L, f=f)
    for jj in range(N):
        f_P = force_P(epsilon=epsilon, R=R, r_arr=r_arr, index=jj, N=N)
        f_arr.append(f_P)
    f_arr = np.array(f_arr)
    return f_arr, f_S_arr
'''

def potential_S(r_arr, L, f):
    r = np.sqrt(np.sum(r_arr**2, axis=1))
    v_s = f/2*(r - L)**2
    v_s[r < L] = 0
    return v_s


## TO DO
'''
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


## TO DO
def generate_v(r_arr, R, L, f, epsilon, N):
    v_arr = []
    v_S = potential_S(r_arr=r_arr, L=L, f=f)
    for j in range(N):
        v_P = potential_P(epsilon=epsilon, R=R, r_arr=r_arr, index=j, N=j)
        v_P_sum = 0
        for jj in range(j):
            v_P_sum = v_P_sum + v_P[jj]
        v_arr.append(v_P_sum+v_S[j])
    v_arr = np.array(v_arr)
    return v_arr
'''

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


def generate_temp(p_arr, N, k, m):
	temp = 2/3/N/k*np.sum(p_arr**2)/2/m
	return temp


def generate_H_Ek(p_arr, m, v_arr):
	e_k = np.sum(p_arr**2)/2/m
	e_c = e_k + np.sum(v_arr)
	return e_c, e_k



def save_to_file_xyz(file_xyz, r_arr, N):
    delimiter = " "
    with open(file_xyz,"a") as file_1:
        for j in range(N):
            file_1.write(str(r_arr[j][0])+delimiter+str(r_arr[j][1])+delimiter+str(r_arr[j][2])+"\n")
        file_1.write("\n \n")


def save_to_file_out(file_out, t, v, E_k, E_c, T, p):
    delimiter = " "
    with open(file_out,"a") as file_1:
        file_1.write(str(t)+delimiter+str(v)+delimiter+str(E_k)+delimiter+str(E_c)+delimiter+str(T)+delimiter+str(E_c)+"\n \n")


def print_momentum_chart(p_arr):
    ay, ax = np.histogram(p_arr[:,0], 30)
    by, bx = np.histogram(p_arr[:,1], 30)
    cy, cx = np.histogram(p_arr[:,2], 30)
    plt.plot(.5*(ax[1:]+ax[:-1]), ay)
    plt.plot(.5*(bx[1:]+bx[:-1]), by)
    plt.plot(.5*(cx[1:]+cx[:-1]), cy)
    plt.show()


def print_3D(arr):
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter3D(arr[:, 0], arr[:, 1], arr[:, 2])
    ax.set_ylim3d(ax.get_xlim3d())
    ax.set_zlim3d(ax.get_xlim3d())
    plt.show()
