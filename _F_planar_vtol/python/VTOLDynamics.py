import numpy as np
import VTOLParam as P
import scipy.integrate as scint
from VTOLStateEOM import eom

class VTOLDynamics:
    def __init__(self, alpha = 0.0):
        # Physical parameters of the VTOL
        self.mc = P.mc  # mass of the central body
        self.mr = P.mr  # mass of each rotor
        self.Jc = P.Jc  # moment of inertia of central body
        self.d = P.d    # distance to each rotor
        self.mu = P.mu  # damping coefficient
        self.g = P.g    # gravitational acceleration
        self.F_wind = P.F_wind  # wind disturbance force

    def f (self, state, tau):
        # Return xdot = f(x,u), the system state update equations
        # re-label states for readability
        z = state[0,0]
        zd = state[1,0]
        h = state[2,0]
        hd = state[3,0]
        theta = state[4,0]
        thetad = state[5,0]
        zdd = 
        hdd = 
        thetadd = 
        xdot = np.array([[zd], [zdd], [hd],[hdd], [thetad], [thetadd]])
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
