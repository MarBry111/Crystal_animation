import numpy as np
import matplotlib.pyplot as plt
import itertools
from mpl_toolkits.mplot3d import Axes3D

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
###TO DOOO
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


def force_S(r_arr, L, f):
    r = np.sqrt(np.sum(r_arr**2, axis=1))
    f_s = r_arr*f*(L-r[:, None])
    f_s[f_s < L] = 0
    return f_s
    '''
    if r < L:
        return r_i*0
    elif r >= L:
        r_out = r_i
        try:
            r_out = r_i*(f*(L - r))
        except TypeError:
            print(r_i,f,L,r)
        return r_out
    '''


def generate_force_arr(r_arr, R, L, f, epsilon, N):
    f_arr = []
    f_S_arr = force_S(r_arr=r_arr, L=L, f=f)
    for jj in range(N):
        f_P = force_P(epsilon=epsilon, R=R, r_arr=r_arr, index=jj, N=N)
        f_arr.append(f_P)
    f_arr = np.array(f_arr)
    return f_arr, f_S_arr


###DONE
def potential_S(r_arr, L, f):
    r = np.sqrt(np.sum(r_arr**2, axis=1))
    v_s = f/2*(r - L)**2
    v_s[r < L] = 0
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
    v_S = potential_S(r_arr=r_arr, L=L, f=f)
    for j in range(N):
        v_P = potential_P(epsilon=epsilon, R=R, r_arr=r_arr, index=j, N=j)
        v_P_sum = 0
        for jj in range(j):
            v_P_sum = v_P_sum + v_P[jj]
        v_arr.append(v_P_sum+v_S[j])
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


def print_3D(force_arr):
    fig = plt.figure()
    ax = Axes3D(fig)
    #ax.scatter3D(r_arr[:, 0], r_arr[:, 1], r_arr[:, 2])
    ax.scatter3D(force_arr[:,0], force_arr[:,1], force_arr[:,2])
    ax.set_ylim3d(ax.get_xlim3d())
    ax.set_zlim3d(ax.get_xlim3d())
    plt.show()
