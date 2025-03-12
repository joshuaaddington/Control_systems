import numpy as np
import hummingbirdParam as P
import dill

dill.settings['recurse'] = True

class HummingbirdDynamics:
    def __init__(self, alpha=0.0):
        # Initial state conditions
        self.state = np.array([
            [P.phi0],  # roll angle
            [P.theta0],  # pitch angle
            [P.psi0],  # yaw angle
            [P.phidot0],  # roll rate
            [P.thetadot0],  # pitch rate
            [P.psidot0],  # yaw rate
        ])

        # Vary the actual physical parameters
        self.ell1 = P.ell1 * (1. + alpha * (2. * np.random.rand() - 1.))
        self.ell2 = P.ell2 * (1. + alpha * (2. * np.random.rand() - 1.))
        self.ell3x = P.ell3x * (1. + alpha * (2. * np.random.rand() - 1.))
        self.ell3y = P.ell3y * (1. + alpha * (2. * np.random.rand() - 1.))
        self.ell3z = P.ell3z * (1. + alpha * (2. * np.random.rand() - 1.))
        self.ellT = P.ellT * (1. + alpha * (2. * np.random.rand() - 1.))
        self.d = P.d * (1. + alpha * (2. * np.random.rand() - 1.))
        self.m1 = P.m1 * (1. + alpha * (2. * np.random.rand() - 1.))
        self.m2 = P.m2 * (1. + alpha * (2. * np.random.rand() - 1.))
        self.m3 = P.m3 * (1. + alpha * (2. * np.random.rand() - 1.))
        self.J1x = P.J1x * (1. + alpha * (2. * np.random.rand() - 1.))
        self.J1y = P.J1y * (1. + alpha * (2. * np.random.rand() - 1.))
        self.J1z = P.J1z * (1. + alpha * (2. * np.random.rand() - 1.))
        self.J2x = P.J2x * (1. + alpha * (2. * np.random.rand() - 1.))
        self.J2y = P.J2y * (1. + alpha * (2. * np.random.rand() - 1.))
        self.J2z = P.J2z * (1. + alpha * (2. * np.random.rand() - 1.))
        self.J3x = P.J3x * (1. + alpha * (2. * np.random.rand() - 1.))
        self.J3y = P.J3y * (1. + alpha * (2. * np.random.rand() - 1.))
        self.J3z = P.J3z * (1. + alpha * (2. * np.random.rand() - 1.))
        self.km = P.km * (1. + alpha * (2. * np.random.rand() - 1.))

        self.param_vals = [
            self.m1, self.m2, self.m3,
            self.J1x, self.J1y, self.J1z,
            self.J2x, self.J2y, self.J2z,
            self.J3x, self.J3y, self.J3z,
            self.ell1, self.ell2, self.ell3x,
            self.ell3y, self.ell3z, self.ellT, self.d
        ]

        self.M = dill.load(open("hb_M_func.pkl", "rb"))
        self.C = dill.load(open("hb_C_func.pkl", "rb"))
        self.dP_dq = dill.load(open("hb_dP_dq_func.pkl", "rb"))
        self.tau = dill.load(open("hb_tau_func.pkl", "rb"))
        self.B = P.beta * np.eye(3)  # friction-based terms

    def update(self, u):
        # Saturate the input force
        u = saturate(u, P.torque_max)
        self.rk4_step(u)  # propagate the state by one time sample
        y = self.h()  # return the corresponding output
        return y

    def f(self, state, u):
        f_l = u[0, 0] * self.km
        f_r = u[1, 0] * self.km

        # Calculating the terms in the equations of motion.
        M = self.M(state, self.param_vals)  # mass matrix
        C = self.C(state, self.param_vals)  # C matrix
        dP_dq = self.dP_dq(state, self.param_vals)  # gravity-based terms
        tau = self.tau(state, [f_l, f_r], self.param_vals)  # generalized force vector

        # Calculate the second derivative of q
        qddot = np.linalg.inv(M) @ (tau - self.B @ state[3:6] - C - dP_dq)

        # Extract first derivative (velocity-based) terms from the state
        phidot = state[3][0]
        thetadot = state[4][0]
        psidot = state[5][0]

        # Define the second derivatives from qddot
        phiddot = qddot[0][0]
        thetaddot = qddot[1][0]
        psiddot = qddot[2][0]

        # Build xdot and return
        xdot = np.array([
            [phidot],
            [thetadot],
            [psidot],
            [phiddot],
            [thetaddot],
            [psiddot]
        ])
        return xdot

    def h(self):
        # Return y = h(x)
        phi = self.state[0][0]
        theta = self.state[1][0]
        psi = self.state[2][0]
        y = np.array([[phi], [theta], [psi]])
        return y

    def rk4_step(self, u):
        # Integrate ODE using Runge-Kutta RK4 algorithm
        F1 = self.f(self.state, u)
        F2 = self.f(self.state + P.Ts / 2 * F1, u)
        F3 = self.f(self.state + P.Ts / 2 * F2, u)
        F4 = self.f(self.state + P.Ts * F3, u)
        self.state += P.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)


def saturate(u, limit):
    for i in range(u.shape[0]):
        if abs(u[i][0]) > limit:
            u[i][0] = limit * np.sign(u[i][0])
    return u