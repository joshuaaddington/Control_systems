import numpy as np
import massParam as P
import scipy.integrate as scint

class MassDynamics:
    def __init__(self. alpha = 0.0)
    #Initial conditions
    self.state = np.array([[P.z0],
                           [P.zdot0]])
    # Mass of mass
    self.m = P.m
    # Spring constant
    self.k = P.k
    # Damping coefficient
    self.b = P.b
    # Sample rate
    self.Ts = P.Ts
    # Max force input
    self.force_limit = P.F_max


    def update(self, u):
        u = saturate(u, self.force_limit)



def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u
