import numpy as np 
import rodMassParam as P
import dill 
dill.settings['recurse'] = True

class rodMassDynamics:
    def __init__(self, alpha=0.0):
        # Initial state conditions
        self.state = np.array([
            [0],  # Theta initial orientation
            [0],  # Thetadot initial velocity
        ])
        # simulation time step
        self.Ts = P.Ts
        # Mass of the pendulum, kg
        self.m = P.m * (1.+alpha*(2.*np.random.rand()-1.))
        # Mass of the cart, kg
        self.ell = P.ell * (1.+alpha*(2.*np.random.rand()-1.))
        # Damping coefficient, Ns
        self.b = P.b * (1.+alpha*(2.*np.random.rand()-1.))
        # gravity constant is well known, don't change.
        self.g = P.g
        self.force_limit = P.tau_max
        
    def update(self, u):
        # This is the external method that takes the input u at time
        # t and returns the output y at time t.
        # saturate the input force
        u = saturate(u, self.force_limit)
        self.rk4_step(u)  # propagate the state by one time sample
        y = self.h()  # return the corresponding output
        return y

    def f(self, state, u):
        # Return xdot = f(x,u)

        # make sure u is a 1D array - this is bit hacky, but we made u a 1D array 
        # when we generated the eom function using "lambdify"
        u = np.array([u]) 

        theta = state[0][0]
        thetadot = state[1][0]
        denom = P.ell**2 * P.m
        thetaddot = u/(denom) - (P.b * thetadot)/(denom) - (P.g * np.cos(theta))/(P.ell) - (P.k1 * theta) / (denom) - (P.k2 * theta**3)/(denom)

        xdot = np.array([[thetadot], [thetaddot[0]]])

        return xdot

    def h(self):
        # return y = h(x) - this is a model of what we can measure, and the output we are
        # trying to control. In this case, we can measure the position of the cart, and the angle
        # of the pendulum. 

        theta = self.state[0][0]
        y = np.array([[theta]])
        return y

    def rk4_step(self, u):
        # Integrate ODE using Runge-Kutta RK4 algorithm
        F1 = self.f(self.state, u)
        F2 = self.f(self.state + self.Ts / 2 * F1, u)
        F3 = self.f(self.state + self.Ts / 2 * F2, u)
        F4 = self.f(self.state + self.Ts * F3, u)
        self.state = self.state + self.Ts / 6 * (F1 + 2*F2 + 2*F3 + F4)

        
def saturate(u, limit):
    if abs(u) > limit:
        u = limit*np.sign(u)
    return u
