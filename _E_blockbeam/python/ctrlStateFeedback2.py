#%%
import numpy as np
import control as cnt
import blockbeamParam as P

class ctrlStateFeedback:
    # dirty derivatives to estimate thetadot
    def __init__(self):
        #  tuning parameters
        tr_z = 1.2
        zeta_z = 0.95
        M = 10
        tr_theta = tr_z/M
        zeta_th  = 0.95

        # State Space Equations
        # xdot = A*x + B*u
        # y = C*x

        # gain calculation
        wn_th = 2.2 / tr_theta  # natural frequency for angle
        wn_z = 2.2 / tr_z  # natural frequency for position
        des_char_poly = np.convolve(
            [1, 2 * zeta_z * wn_z, wn_z**2],
            [1, 2 * zeta_th * wn_th, wn_th**2])
        self.des_poles = np.roots(des_char_poly)

        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(P.A, P.B)) != 4:
            print("The system is not controllable")
        else:
            self.K = (cnt.place(P.A, P.B, self.des_poles))
            self.kr = -1.0 / (P.C @ np.linalg.inv(P.A - P.B @ self.K) @ P.B)
            print("Check denominator: ", P.C @ np.linalg.inv(P.A - P.B @ self.K) @ P.B)
        print('K: ', self.K)
        print('kr: ', self.kr)
        print(self.des_poles)

    def update(self, z_r, x):
        # Compute the state feedback controller
        F_tilde = -self.K @ (x - P.xe) + self.kr * (z_r - P.ze)

        # compute total torque
        F = saturate(P.Fe + F_tilde[0][0], P.F_max)
        return F


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u


ctrlTest = ctrlStateFeedback()
#%%