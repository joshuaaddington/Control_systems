import numpy as np
import control as cnt
import massParam as P

class controllerStateSpace:
    # dirty derivatives to estimate zdot
    def __init__(self):
        # LQR Control Matrices
        Q = np.array([[100.0, 0.0, 0.0],   # z
                      [0.0, 0.1, 0.0],   # z_dot
                      [0.0, 0.0, 1000.0]])  # e_z_int
        R = np.array([[5.0]])  # [F]

        # State Space Equations
        A = np.array([[0, 1], [-P.k1/P.m, -P.b/P.m]])
        B = np.array([[0.0],
                      [1/(P.m)]])        
        C = np.array([[1.0, 0.0]])

        Cr = C
        # form augmented system
        A1 = np.vstack((np.hstack((A, np.zeros((np.size(A,1),1)))), 
                        np.hstack((-Cr, np.array([[0.0]]))) ))
        B1 = np.vstack( (B, 0.0) )

        # gain calculation
        # zeta = .95
        # wn = 3.178  # natural frequency
        des_poles = np.array([-2 + 2j, -2 - 2j, -5])

        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A, B)) != A.shape[0]:
            print("The system is not controllable")
        else:
            # Gains Matrix for normal State Space Control
            # K1 = cnt.place(A1, B1, des_poles)
            # self.K = K1[0][0:2]
            # self.ki = K1[0][2]
            # Gains Matrix for LQR Control
            K1_lqr, _, _ = cnt.lqr(A1, B1, Q, R)
            self.K = K1_lqr[0][0:2].reshape(1,2)
            self.ki = K1_lqr[0][2]
        # Augmented Matrices for observer Design
        A2 = np.concatenate((
                            np.concatenate((A, B), axis=1),
                            np.zeros((1, 3))),
                        axis=0)
        B2 = np.concatenate((B, np.zeros((1, 1))), axis=0)
        C2 = np.concatenate((C, np.zeros((1, 1))), axis=1)
        des_obsv_poles = np.array([-10 + 2j, -10 - 2j, -2])
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A2.T, C2.T)) != 3:
            print("The system is not observerable")
        else:
            L2 = cnt.place(A2.T, C2.T, des_obsv_poles).T
        
        print('K: ', self.K)
        print('ki ', self.ki)
        print('L^T: ', L2.T)
        print("obsv_poles:", des_obsv_poles)
        #--------------------------------------------------
        # variables to implement integrator
        self.obsv_state = np.array([
            [0.0],  # theta_hat_0
            [0.0],  # thetadot_hat_0
            [0.0],  # estimate of the disturbance
        ])
        self.F_d1 = 0
        self.L = L2
        self.A = A2
        self.B = B1
        self.C = C2
        self.integrator = 0.0  # integrator
        self.error_d1 = 0.0  # error signal delayed by 1 sample

    def update(self, z_r, y_m):
        x_hat, d_hat = self.update_observer(y_m)
        z_hat = x_hat[0][0]
        # compute feedback linearizing torque F_fl
        error = z_r - z_hat
        self.integrator = self.integrator \
                          + (P.Ts / 2.0) * (error + self.error_d1)
        self.error_d1 = error
        # Compute the state feedback controller
        F_tilde = -self.K @ x_hat - self.ki * self.integrator - d_hat
        F = F_tilde + P.F_e

        # compute total torque
        F = self.saturate(F[0], P.F_max)
        self.F_d1 = F
        return F, x_hat, d_hat
    
    def update_observer(self, y_m):
        # update the observer using RK4 integration
        F1 = self.observer_f(self.obsv_state, y_m)
        F2 = self.observer_f(self.obsv_state + P.Ts / 2 * F1, y_m)
        F3 = self.observer_f(self.obsv_state + P.Ts / 2 * F2, y_m)
        F4 = self.observer_f(self.obsv_state + P.Ts * F3, y_m)
        self.obsv_state = self.obsv_state + P.Ts / 6 *  (F1 + 2*F2 + 2*F3 + F4)
        x_hat = self.obsv_state[0:2]
        d_hat = self.obsv_state[2][0]
        return x_hat, d_hat
    
    def observer_f(self, x_hat, y_m):
        # compute feedback linearizing torque tau_fl
        z_hat = x_hat[0][0]
        F_e = P.F_e
        xhat_dot = self.A @ x_hat + self.B * (self.F_d1 - F_e) + self.L * (y_m - self.C @ x_hat)
        return xhat_dot

    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit * np.sign(u)
        return u

