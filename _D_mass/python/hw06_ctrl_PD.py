import numpy as np
import massParam as P
from hw06_mass_Kp_Kd import *


class ctrlPD:
    def __init__(self):
        #  tuning parameters
        #tr = 0.8 # part (a)
        tr = 0.37 # tuned for faster rise time before saturation.
        zeta = 0.707

        # desired natural frequency
        wn = 0.5 * np.pi / (tr * np.sqrt(1 - zeta**2))
        alpha1 = 2.0 * zeta * wn
        alpha0 = wn**2

        # compute PD gains
        self.kd = P.m * (alpha1 - (P.b/P.m))
        self.kp = P.m * (alpha0 - (P.k/P.m))
        print('kp: ', self.kp)
        print('kd: ', self.kd)

    def update(self, x_ref, state):
        x = state[0][0]
        x_dot = state[1][0]

        # compute feedback linearizing torque tau_fl
        F_fl = 0

        # compute the linearized torque using PD
        F_tilde = self.kp * (x_ref - x) \
                    - self.kd * x_dot
        
        # compute total torque
        F = F_fl + F_tilde
        F = saturate(F, P.F_max)
        return F


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u