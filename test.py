import itertools
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

n = 2      # number atoms on each axis
a = 0.38    # nm, distance between atoms
N = n ** 3  # number of all atoms
m = 39.948  # mass


def generate_r_arr(b0, b1, b2, n):
    i_arr = np.arange(n)
    r_arr = np.array(list(p for p in itertools.product(i_arr, repeat=3)), dtype=float)
    for j in range(3):
        r_arr[:, j] = (r_arr[:, 0]-(n-1)/2)*b0[j] + (r_arr[:, 1]-(n-1)/2)*b1[j] + (r_arr[:, 2]-(n-1)/2)*b2[j]
    return r_arr


def print_3D_r(r_arr):
    fig = plt.figure()
    ax = Axes3D(fig)
    ax.scatter3D(r_arr[:, 0], r_arr[:, 1], r_arr[:, 2])
    ax.set_ylim3d(ax.get_xlim3d())
    ax.set_zlim3d(ax.get_xlim3d())
    plt.show()


b_0 = np.array([a, 0, 0])
b_1 = np.array([a / 2, a * np.sqrt(3) / 2, 0])
b_2 = np.array([a / 2, a * np.sqrt(3) / 6, a * np.sqrt(2) / np.sqrt(3)])
r_arr1 = generate_r_arr(b0=b_0, b1=b_1, b2=b_2, n=n)

p = np.zeros(shape=(3, 3))
p = p-np.array([1,1,1])
print(p)