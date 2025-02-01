import numpy as np
import blockbeamParam as P
import scipy.integrate as scint
from blockbeamStateEOM import eom

class blockbeamDynamics:
    def __init__(self, alpha = 0.0):
        #Initial conditions
        self.state = np.array([[P.z0],
                            [P.zdot0],
                            [P.theta0],
                            [P.thetadot0]])
        # mass of blockbeam
        self.m1 = P.m1
        self.m2 = P.m2
        # Sample rate
        self.Ts = P.Ts
        # Max force input
        self.torque_limit = P.F_max

    def f (self, state, tau):
        # Return xdot = f(x,u), the system state update equations
        # re-label states for readability
        z = state[0,0]
        zd = state[1,0]
        theta = state[2,0]
        thetad = state[3,0]
        zdd = (z*thetad**2 - 9.8*np.sin(theta))
        thetadd = ((3 * tau * np.cos(theta) - 58.8 * self.m1 * z * np.cos(theta) - 12.0 * self.m1 * z * thetad * zd - 14.7 * self.m2 * np.cos(theta)) / (6.0 * self.m1 * z**2 + 0.5 * self.m2))
        xdot = np.array([[zd], [zdd], [thetad], [thetadd]])
        return xdot
        
    def rk4_step(self, u):
        # Integrate ODE using Runge-Kutta RK4 algorithm
        F1 = self.f(self.state, u)
        F2 = self.f(self.state + self.Ts / 2 * F1, u)
        F3 = self.f(self.state + self.Ts / 2 * F2, u)
        F4 = self.f(self.state + self.Ts * F3, u)
        self.state = self.state + self.Ts / 6 * (F1 + 2*F2 + 2*F3 + F4)

    def update(self, u):
        u = saturate(u, self.torque_limit)
        self.rk4_step(u)
        y = np.array(self.state[0,0])
        return y



def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u
