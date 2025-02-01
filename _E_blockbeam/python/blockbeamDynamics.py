#%%
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
        self.force_limit = P.F_max

    def f (self, state, tau):
        # Return xdot = f(x,u), the system state update equations
        # re-label states for readability
        z = state[0][0]
        zd = state[1][0]
        theta = state[2][0]
        thetad = state[3][0]
        zdd = eom(state, tau, self.m1, self.m2)[0][0]
        thetadd = eom(state, tau, self.m1, self.m2)[1][0]
        print("zd = ", zd, "type is: ", type(zd), "shape is: ", zd.shape)
        print("zdd = ", zdd, "type is: ", type(zdd), "shape is: ", zdd.shape)
        print("thetad = ", thetad, "type is: ", type(thetad), "shape is: ", thetad.shape)
        print("thetadd = ", thetadd, "type is: ", type(thetadd), "shape is: ", thetadd.shape)
        xdot = np.array([[zd], [zdd], [thetad], [thetadd]])
        print(xdot)
        return xdot
        
    def rk4_step(self, u):
        # Integrate ODE using Runge-Kutta RK4 algorithm
        F1 = self.f(self.state, u)
        F2 = self.f(self.state + self.Ts / 2 * F1, u)
        F3 = self.f(self.state + self.Ts / 2 * F2, u)
        F4 = self.f(self.state + self.Ts * F3, u)
        self.state = self.state + self.Ts / 6 * (F1 + 2*F2 + 2*F3 + F4)

    def update(self, u):
        u = np.array(saturate(u, self.force_limit)).reshape(1,1)
        print("u = ", u, "type is: ", type(u), "shape is: ", u.shape)
        self.rk4_step(u)
        y = np.array(self.state[0,0])
        return y



def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u

#%% 
blockbeam = blockbeamDynamics()
blockbeam.state = np.array([[0.25], [0.1], [np.pi/180], [3]])
u = 1
blockbeam.update(u)
# %%
