import numpy as np
import control as cnt
import rodMassParam as P

class ctrlDisturbanceObserver:
    # dirty derivatives to estimate thetadot
    def __init__(self):
        #  tuning parameters
        tr = 0.4
        zeta = 0.9
        integrator_pole = -10
        zeta_obs = 0.707       # damping ratio for observer
        dist_obsv_pole = -15  # pole for disturbance observer

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

        A2 = np.concatenate((
                            np.concatenate((A, B), axis=1),
                            np.zeros((1, 3))),
                        axis=0)
        B2 = np.concatenate((B, np.zeros((1, 1))), axis=0)
        C2 = np.concatenate((C, np.zeros((1, 1))), axis=1)
        des_obsv_poles = np.hstack(((des_poles * 5) [0:2], dist_obsv_pole))
        print("des_obsv_poles:", des_obsv_poles)
        if np.linalg.matrix_rank(cnt.ctrb(A2.T, C2.T)) != 3:
            print("The system is not observerable")
        else:
            L2 = cnt.place(A2.T, C2.T, des_obsv_poles).T

        print('K: ', self.K)
        print('ki: ', self.ki)
        print(des_poles)
        print('L^T: ', L2.T)
        self.integrator = 0.0  # integrator
        self.error_d1 = 0.0  # error signal delayed by 1 sample
        self.obsv_state = np.array([
            [0.0],  # theta_hat_0
            [0.0],  # thetadot_hat_0
            [0.0],  # estimate of the disturbance
        ])
        self.L = L2
        self.A = A2
        self.B = B1
        self.C = C2
        self.tau_d1 = 0.0

    def update(self, theta_r, y_m):
        x_hat, d_hat = self.update_observer(y_m)
        theta_hat = x_hat[0][0]
        # integrate error
        error = theta_r - theta_hat
        self.integrator = self.integrator + (P.Ts / 2.0) * (error + self.error_d1)
        self.error_d1 = error
        # compute feedback linearizing torque tau_fl
        tau_fl = 0
        # Compute the state feedback controller
        tau_tilde = -self.K @ x_hat - self.ki * self.integrator
        # compute total torque
        tau = saturate(tau_fl + tau_tilde[0], P.tau_max)
        self.tau_d1 = tau
        return tau, x_hat, d_hat
    
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
        theta_hat = x_hat[0][0]
        tau_fl = 0
        # xhatdot = A*xhat + B*(u-ue) + L(y-C*xhat)
        xhat_dot = self.A @ x_hat\
                   + self.B * (self.tau_d1 - tau_fl)\
                   + self.L * (y_m - self.C @ x_hat)
        return xhat_dot

def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u


