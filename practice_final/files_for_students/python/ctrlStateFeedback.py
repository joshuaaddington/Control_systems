import numpy as np
import control as cnt
import rodMassParam as P

class ctrlStateFeedback:
    # dirty derivatives to estimate thetadot
    def __init__(self):
        #  tuning parameters
        tr = 0.4
        zeta = 0.9
        integrator_pole = -10

        # State Space Equations
        # xdot = A*x + B*u
        # y = C*x

        A = np.array([[0.0, 1.0],
                      [-.02 / (P.ell**2 * P.m), -P.b/(P.ell**2*P.m)]])
        B = np.array([[0.0],
                      [1/(P.ell**2 * P.m)]])        
        C = np.array([[1.0, 0.0]])

        A1 = np.vstack((np.hstack((A, np.zeros((np.size(A,1),1)))), 
                        np.hstack((-C, np.array([[0.0]]))) ))
        B1 = np.vstack( (B, 0.0) )

        # gain calculation
        wn = 37.07  # natural frequency
        des_char_poly = np.convolve([1, 2 * zeta * wn, wn**2], [1, -integrator_pole])
        des_poles = np.roots(des_char_poly)
        #des_poles = [-1.0, -1.0]

        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A1, B1)) != A1.shape[0]:
            print("The system is not controllable")
        else:
            K1 = cnt.place(A1, B1, des_poles)
            self.K = K1[0][0:2]
            self.ki = K1[0][2]
        print('K: ', self.K)
        print('ki: ', self.ki)
        print(des_poles)
        self.integrator = 0.0  # integrator
        self.error_d1 = 0.0  # error signal delayed by 1 sample

    def update(self, theta_r, x):
        theta = x[0][0]
        # integrate error
        error = theta_r - theta
        self.integrator = self.integrator \
                          + (P.Ts / 2.0) * (error + self.error_d1)
        self.error_d1 = error
        # compute feedback linearizing torque tau_fl
        tau_fl = 0
        # Compute the state feedback controller
        tau_tilde = -self.K @ x - self.ki * self.integrator
        # compute total torque
        tau = saturate(tau_fl + tau_tilde[0], P.tau_max)
        return tau


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u


