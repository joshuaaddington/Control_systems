import numpy as np
import control as cnt
import massParam as P

class ctrlDisturbanceObserver :
    def __init__(self):
        #--------------------------------------------------
        # State Feedback Control Design
        #--------------------------------------------------
        #  tuning parameters
        tr = 0.8
        zeta = 0.707
        integrator_pole = [-2]
        tr_obs = tr/10           # natural frequency for observer
        zeta_obs = 0.707       # damping ratio for observer
        dist_obsv_pole = [-1]  # pole for disturbance observer

        # State Space Equations
        # xdot = A*x + B*u
        # y = C*x
        A = np.array([[0.0, 1.0],
                      [-P.k / P.m, -P.b / P.m]])
        B = np.array([[0.0],
                     [1.0 / P.m]])
        C = np.array([[1.0, 0.0]])
        # form augmented system
        Cr = C
        A1 = np.zeros((3, 3))
        A1[0:2, 0:2] = A
        A1[2, 0:2] = -Cr

        B1 = np.zeros((3, 1))
        B1[0:2, :] = B

        # gain calculation
        wn = 2.2 / tr  # natural frequency
        #wn = 0.5*np.pi/(tr*np.sqrt(1-zeta**2)) # natural frequency
        des_char_poly = np.convolve([1, 2 * zeta * wn, wn**2],
                                    np.poly(integrator_pole))
        des_poles = np.roots(des_char_poly)
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A1, B1)) != A1.shape[0]:
            print("The system is not controllable")
        else:
            K1 = cnt.place(A1, B1, des_poles)
            self.K = K1[0][0:2]
            self.ki = K1[0][2]
        # observer design
        # Augmented Matrices
        self.A2 = np.zeros((A.shape[0]+1, A.shape[1]+B.shape[1]))
        self.A2[0:A.shape[0], 0:A.shape[1]] = A
        self.A2[0:B.shape[0], A.shape[1]] = B.flatten()
        self.B2 = np.zeros((B.shape[0]+B.shape[1], B.shape[1]))
        self.B2[0:B.shape[0], :] = B
        self.C2 = np.concatenate((C, np.zeros((1, 1))), axis=1)


        wn_obs = 2.2 / tr_obs  # natural frequency for observer
        des_char_est = np.array([1., 2.*zeta*wn_obs, wn_obs**2.])
        des_obsv_char_poly = np.convolve([1, 2 * zeta_obs * wn_obs, wn_obs**2],
                                         np.poly(dist_obsv_pole))
        des_obsv_poles = np.roots(des_obsv_char_poly)
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(self.A2.T, self.C2.T)) != 3:
            print("The system is not observerable")
        else:
            self.L2 = cnt.acker(self.A2.T, self.C2.T, des_obsv_poles).T
        print('K: ', self.K)
        print('ki ', self.ki)
        print('L^T: ', self.L2.T)
        #--------------------------------------------------
        # variables to implement integrator
        self.integrator = 0.0  # integrator
        self.error_d1 = 0.0  # error signal delayed by 1 sample
        self.obsv_state = np.array([
            [0.0],  # theta_hat_0
            [0.0],  # thetadot_hat_0
            [0.0],  # estimate of the disturbance
        ])
        self.force_d1 = 0.0  # control torque, delayed 1 sample
        self.L = self.L2
        self.A = self.A2
        self.B = B1
        self.C = self.C2

    def update(self, z_r, y_m):
        # update the observer and extract theta_hat
        x_hat, d_hat = self.update_observer(y_m)
        z_hat = x_hat[0][0]
        # integrate error
        error = z_r - z_hat
        self.integrator = self.integrator + (P.Ts / 2.0) * (error + self.error_d1)
        self.error_d1 = error
        # compute feedback linearizing torque tau_fl
        force_unsat = -self.K @ x_hat - self.ki * self.integrator - d_hat


        # compute total torque
        force = saturate(force_unsat, P.F_max)
        if self.ki != 0.0:
            self.integrator = self.integrator + P.Ts/self.ki*(force-force_unsat)
        self.force_prev = force
        self.force_d1 = force
        return force, x_hat, d_hat

    def update_observer(self, y_m):
        # update the observer using RK4 integration
        F1 = self.observer_f(self.obsv_state, y_m)
        F2 = self.observer_f(self.obsv_state + P.Ts / 2 * F1, y_m)
        F3 = self.observer_f(self.obsv_state + P.Ts / 2 * F2, y_m)
        F4 = self.observer_f(self.obsv_state + P.Ts * F3, y_m)
        self.obsv_state = self.obsv_state + P.Ts / 6 *  (F1 + 2*F2 + 2*F3 + F4)
        x_hat = self.obsv_state[0:2]
        d_hat = self.obsv_state[2,0]
        return x_hat, d_hat

    def observer_f(self, obsv_state, y_m):
        # compute feedback linearizing torque tau_fl
        observer_state_dot = self.A2 @ obsv_state \
                   + self.B2 * self.force_d1 \
                   + self.L2 @ (y_m - self.C2 @ obsv_state)
        return observer_state_dot


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u


