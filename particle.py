import numpy as np
import random
from main import k, n


# building class here, not finished
class Particle:
    m = 0
    T = 0
    r = np.array([])
    energy_kin = np.array([])
    p = np.array([])

    @staticmethod
    def energy_kin_ax(temp):
        lam = random.random()
        return -1 / 2 * k * temp * np.log(lam)

    @staticmethod
    def momentum_ax(energy_kin, m):
        temp = random.random()
        char = 0
        if temp > 0.5:
            char = 1
        else:
            char = -1
        return char * np.sqrt(2 * m * energy_kin)

    def __init__(self, mass, temp0, b0, b1, b2, iter1, iter2, iter3):
        # temperature
        self.m = mass
        self.T = temp0
        # x, y, z coordinates
        x = b0[0] * (iter1 - (n - 1) / 2) + b1[0] * (iter2 - (n - 1) / 2) + b2[0] * (iter3 - (n - 1) / 2)
        y = b0[1] * (iter1 - (n - 1) / 2) + b1[1] * (iter2 - (n - 1) / 2) + b2[1] * (iter3 - (n - 1) / 2)
        z = b0[2] * (iter1 - (n - 1) / 2) + b1[2] * (iter2 - (n - 1) / 2) + b2[2] * (iter3 - (n - 1) / 2)
        self.r = np.append(self.r, [x, y, z])
        # Energy in each axis
        e_kin_x = Particle.energy_kin_ax(self.T)
        e_kin_y = Particle.energy_kin_ax(self.T)
        e_kin_z = Particle.energy_kin_ax(self.T)
        self.energy_kin = np.append(self.energy_kin, [e_kin_x, e_kin_y, e_kin_z])
        # momentum in each axis
        p_x = Particle.momentum_ax(e_kin_x, self.m)
        p_y = Particle.momentum_ax(e_kin_y, self.m)
        p_z = Particle.momentum_ax(e_kin_z, self.m)
        np.append(self.p, [p_x, p_y, p_z])