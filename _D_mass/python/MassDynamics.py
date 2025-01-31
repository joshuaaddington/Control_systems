import numpy as np
import massParam as P
import scipy.integrate as scint
from massStateEOM import eom

class massDynamics:
    def __init__(self, alpha = 0.0):
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

    def f (self, state, tau):
        # Return xdot = f(x,u), the system state update equations
        # re-label states for readability
        z = state[0,0]
        zd = state[1,0]
        zdd = (1/self.m) * (tau - self.b*zd) - (self.k*z)
        xdot = np.array([[zd], [zdd]])
        return xdot
        
    def rk4_step(self, u):
        # Integrate ODE using Runge-Kutta RK4 algorithm
        F1 = self.f(self.state, u)
        F2 = self.f(self.state + self.Ts / 2 * F1, u)
        F3 = self.f(self.state + self.Ts / 2 * F2, u)
        F4 = self.f(self.state + self.Ts * F3, u)
        self.state = self.state + self.Ts / 6 * (F1 + 2*F2 + 2*F3 + F4)

    def update(self, u):
        u = saturate(u, self.force_limit)
        self.rk4_step(u)
        y = np.array(self.state[0,0])
        return y



def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u
