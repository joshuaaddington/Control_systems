#%%
import sympy as sp
import hummingbirdParam as P
from sympy.physics.vector import dynamicsymbols
from IPython.display import display, Math
from sympy.physics.vector import vlatex
from sympy import sin, cos, tan
import numpy as np

# longitudinal controller class
class ctrl_longitudinal:
    def __init__(self):
        # Define the Transfer Function
        theta = dynamicsymbols('theta') # we have to tell SymPy that these are functions of time
        theta_dot = sp.diff(theta)
        s = sp.symbols('s') # Laplace symbol
        
        self.F_ctrl = sp.symbols('F_ctrl')
        F_fl = ((P.m1*P.ell1 + P.m2*P.ell2)*9.810*cos(theta))/P.ellT

        self.b0 = (P.ellT/(P.m1*P.ell1**2 + P.m2*P.ell2**2 + P.J1y + P.J2y))
        theta_ddot = self.b0*self.F_ctrl

        x = sp.Matrix([theta, theta_dot])
        x_dot = sp.Matrix([theta_dot, theta_ddot])

        self.F = self.F_ctrl + F_fl

        self.open_loop = self.b0 * self.F_ctrl / s^2

        #  tuning parameters
        #tr = 0.8 # part (a)
        tr = 1 # tuned for faster rise time before saturation.
        zeta = 0.707

        # desired natural frequency
        wn = 0.5 * np.pi / (tr * np.sqrt(1 - zeta**2))
        alpha1 = 2.0 * zeta * wn
        alpha0 = wn**2

        # compute PD gains
        self.kp = alpha0 / self.b0
        self.kd = alpha1 / self.b0

        def update(self, theta_d, state):
            theta = state[0][0]
            thetadot = state[1][0]

            # compute the linearized torque using PD
            F_ctrl = self.kp * (theta_d - theta) - self.kd * thetadot
            
            # compute total torque
            F = F_ctrl + F_fl
            tau = saturate(F, P.force_max)
            return tau

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u



ctrl_test = ctrl_longitudinal()


# %%
