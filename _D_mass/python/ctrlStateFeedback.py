#%%
import numpy as np
import control as cnt
import massParam as P

class ctrlStateFeedback:
    # dirty derivatives to estimate thetadot
    def __init__(self):
        #  tuning parameters
        tr = 2.0
        zeta = 0.707

        # State Space Equations
        # xdot = A*x + B*u
        # y = C*x

        # gain calculation
        wn = 2.2 / tr  # natural frequency
        des_char_poly = [1, 2 * zeta * wn, wn**2]
        des_poles = np.roots(des_char_poly)

        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(P.A, P.B)) != 2:
            print("The system is not controllable")
        else:
            self.K = (cnt.place(P.A, P.B, des_poles))
            self.kr = -1.0 / (P.C @ np.linalg.inv(P.A - P.B @ self.K) @ P.B)
        print('K: ', self.K)
        print('kr: ', self.kr)
        print(des_poles)

    def update(self, z_r, x):
        # compute feedback linearizing torque tau_fl
        force_fl = 0

        # Compute the state feedback controller
        force_tilde = -self.K @ x + self.kr * z_r

        # compute total torque
        force = saturate(force_fl + force_tilde[0][0], P.F_max)

        return force


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u


