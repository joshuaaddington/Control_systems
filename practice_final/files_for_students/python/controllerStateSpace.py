import numpy as np
import control as cnt
import rodMassParam as P

class controllerStateSpace:
    # dirty derivatives to estimate thetadot
    def __init__(self):
        A = np.array([[0.0, 1.0],
                      [-.02 / (P.ell**2 * P.m), -P.b/(P.ell**2*P.m)]])
        B = np.array([[0.0],
                      [1/(P.ell**2 * P.m)]])        
        C = np.array([[1.0, 0.0]])

        # gain calculation
        zeta = .9
        wn = 37.07  # natural frequency
        des_char_poly = [1, 2 * zeta * wn, wn**2]
        des_poles = np.roots(des_char_poly)
        #des_poles = [-1.0, -1.0]

        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A, B)) != A.shape[0]:
            print("The system is not controllable")
        else:
            self.K = cnt.place(A, B, des_poles)
            self.kr = -1.0 / (C @ np.linalg.inv(A - B @ self.K) @ B)
        print('K: ', self.K)
        print('kr: ', self.kr)
        print(des_poles)

    def update(self, theta_r, x):
        theta = x[0][0]
        # compute feedback linearizing torque tau_fl

        # Compute the state feedback controller
        tau_tilde = -self.K @ x + self.kr * theta_r

        # compute total torque
        tau = self.saturate(tau_tilde[0][0], P.tau_max)

        return tau

    def saturate(u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u)
        return u

