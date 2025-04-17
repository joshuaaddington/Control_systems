import numpy as np
import control as cnt
import VTOLParam as P

class ctrlLQRVTOL:
    def __init__(self):
        #--------------------------------------------------
        # State Feedback Control Design
        #--------------------------------------------------
        #  tuning parameters
        tr_h    = 5.0 
        zeta_h  = 0.95
        wn_h    = 0.5*np.pi/(tr_h*np.sqrt(1-zeta_h**2))
        tr_z    = 5.0 
        zeta_z  = 0.95
        wn_z    = 0.5*np.pi/(tr_z*np.sqrt(1-zeta_z**2))
        tr_th   = tr_z/10.0 
        zeta_th = 0.95
        wn_th   = 0.5*np.pi/(tr_th*np.sqrt(1-zeta_th**2))
        integrator_h = 1.0
        integrator_z = 1.0
        wn_h_obs    = 10.0*wn_h
        wn_z_obs    = 10.0*wn_z
        wn_th_obs   = 5.0*wn_th
        dist_obsv_pole_lon = 10.0
        dist_obsv_pole_lat = 10.0

        Q_lon = np.array([[15.0, 0.0, 0.0], # h
                          [0.0, 20.0, 0.0], # h_dot
                          [0.0, 0.0, 10.0]]) # e_h_int
        R_lon = np.array([[0.5]]) # [F]

        Q_lat = np.array([[2, 0, 0, 0, 0], # z
                          [0, 1, 0, 0, 0], # z_dot
                          [0, 0, 15, 0, 0], # theta
                          [0, 0, 0, 10, 0], # theta_dot
                          [0, 0, 0, 0, 4]]) # e_z_int
        R_lat = np.array([[1.0]]) # [tau]

        # State Space Equations
        # xdot = A*x + B*u
        # y = C*x
        self.Fe = (P.mc + 2.0 * P.mr) * P.g  # equilibrium force 
        A_lon = np.array([[0.0, 1.0],
                        [0.0, 0.0]])
        B_lon = np.array([[0.0],
                        [1.0 / (P.mc + 2.0 * P.mr)]])
        C_lon = np.array([[1.0, 0.0]])

        A_lat = np.array([[0.0, 0.0, 1.0, 0.0],
                        [0.0, 0.0, 0.0, 1.0],
                        [0.0, -self.Fe / (P.mc + 2.0 * P.mr), \
                         -(P.mu / (P.mc + 2.0 * P.mr)), 0.0],
                        [0.0, 0.0, 0.0, 0.0]])
        B_lat = np.array([[0.0],
                        [0.0],
                        [0.0],
                        [1.0 / (P.Jc + 2 * P.mr * P.d ** 2)]])
        C_lat = np.array([[1.0, 0.0, 0.0, 0.0],
                          [0.0, 1.0, 0.0, 0.0]])
        Cr_lon = C_lon
        A1_lon = np.zeros((A_lon.shape[0]+1, A_lon.shape[0]+1))        
        A1_lon[0:2, 0:2] = A_lon
        A1_lon[2, 0:2] = -Cr_lon
        B1_lon = np.vstack((B_lon, np.zeros((1,1))))
        Cr_lat = C_lat[0,:].reshape(1,4)
        A1_lat = np.vstack((
                np.hstack((A_lat, np.zeros((4,1)))),
                np.hstack((-Cr_lat, np.zeros((1,1))))))
        B1_lat = np.vstack((B_lat, np.zeros((1,1))))

        # gain calculation
        # des_char_poly_lon = np.convolve([1.0, 2.0 * zeta_h * wn_h, wn_h ** 2],
        #                                 [1, integrator_h])
        # des_poles_lon = np.roots(des_char_poly_lon)
        
        # des_char_poly_lat = np.convolve(
        #     np.convolve([1.0, 2.0 * zeta_z * wn_z, wn_z ** 2],
        #                 [1.0, 2.0 * zeta_th * wn_th, wn_th ** 2]),
        #     [1, integrator_z])
        # des_poles_lat = np.roots(des_char_poly_lat)
        # Compute the gains if the system is controllable
        if np.linalg.matrix_rank(cnt.ctrb(A1_lon, B1_lon)) != A1_lon.shape[0]:
            print("The longitudinal system is not controllable")
        else:
            K1_lon,_,_ = cnt.lqr(A1_lon, B1_lon,Q_lon, R_lon)
            self.K_lon = K1_lon[0][0:2]
            self.ki_lon = K1_lon[0][2]

        if np.linalg.matrix_rank(cnt.ctrb(A1_lat, B1_lat)) != A1_lat.shape[0]:
            print("The lateral system is not controllable")
        else:
            K1_lat,_,_ = cnt.lqr(A1_lat, B1_lat, Q_lat, R_lat)
            self.K_lat = K1_lat[0][0:4]
            self.ki_lat = K1_lat[0][4]
        # observer design
        # Augmented Matrices
        self.A2_lon = np.concatenate((
                np.concatenate((A_lon, B_lon), axis=1),
                np.zeros((1, 3))),
                axis=0)
        self.B2_lon = np.vstack((B_lon, np.zeros((1, 1))))
        self.C2_lon = np.concatenate((C_lon, np.zeros((1, 1))), axis=1)
        self.A2_lat = np.concatenate((
                np.concatenate((A_lat, B_lat), axis=1),
                np.zeros((1, 5))),
                axis=0)
        self.B2_lat = np.vstack((B_lat, np.zeros((1, 1))))
        self.C2_lat = np.concatenate((C_lat, np.zeros((2, 1))), axis=1)


        des_obs_char_poly_lon = np.convolve([1.0, 2.0*zeta_h*wn_h_obs, wn_h_obs**2],
                                            [1, dist_obsv_pole_lon])
        des_obs_poles_lon = np.roots(des_obs_char_poly_lon)

        des_obs_char_poly_lat = np.convolve(
                                    np.convolve([1.0, 2.0*zeta_z*wn_z_obs, wn_z_obs**2],
                                                [1.0, 2.0*zeta_th*wn_th_obs, wn_th_obs**2]),
                                    [1, dist_obsv_pole_lat])
        des_obs_poles_lat = np.roots(des_obs_char_poly_lat)

        if np.linalg.matrix_rank(cnt.ctrb(self.A2_lon.T, self.C2_lon.T)) != self.A2_lon.shape[0]:
            print("The longitudinal system is not observable")
        else:
            self.L2_lon = cnt.place(self.A2_lon.T, self.C2_lon.T, des_obs_poles_lon).T

        if np.linalg.matrix_rank(cnt.ctrb(self.A2_lat.T, self.C2_lat.T)) != self.A2_lat.shape[0]:
            print("The lateral system is not observable")
        else:
            self.L2_lat = cnt.place(self.A2_lat.T, self.C2_lat.T, des_obs_poles_lat).T
        print('K_lon: ', self.K_lon)
        print('ki_lon: ', self.ki_lon)
        print('L_lon^T: ', self.L2_lon.T)
        print('K_lat: ', self.K_lat)
        print('ki_lat: ', self.ki_lat)
        print('L_lat^T: ', self.L2_lat.T)
        #--------------------------------------------------
        # variables to implement integrator
        self.observer_state_lon = np.array([[0.0],
                                            [0.0], 
                                            [0.0]])
        self.observer_state_lat = np.array([[0.0], 
                                            [0.0], 
                                            [0.0], 
                                            [0.0], 
                                            [0.0]])

        self.integrator_z = 0.0      # integrator on position z
        self.error_z_prev = 0.0        # error signal delayed by 1 sample
        self.integrator_h = 0.0      # integrator on altitude h
        self.error_h_prev = 0.0        # error signal delayed by 1 sample
        self.F_prev = 0.0  # Force signal delayed by 1 sample
        self.tau_prev = 0.0  # torque signal delayed by 1 sample
        self.F_limit = P.F_max * 2.0
        self.tau_limit = P.F_max * P.d * 2.0
        self.Ts = P.Ts

    def update(self, r, y):
        # assign reference inputs and measured outputs
        z_r = r[0][0]
        h_r = r[1][0]
        y_lat = np.array([[y[0][0]],[y[2][0]]])
        y_lon = y[1][0]

        # update the observers
        xhat_lat, dhat_lat = self.update_lat_observer(y_lat)
        xhat_lon, dhat_lon = self.update_lon_observer(y_lon)
        z_hat = xhat_lat[0][0]
        h_hat = xhat_lon[0][0]
        theta_hat = xhat_lat[1][0]

        # integrate error for h and z
        error_z = z_r - z_hat
        self.integrator_z += (P.Ts/2.0)*(error_z + self.error_z_prev)
        self.error_z_prev = error_z
        error_h = h_r - h_hat
        self.integrator_h += (P.Ts/2.0)*(error_h + self.error_h_prev)
        self.error_h_prev = error_h

        # Compute the state feedback controllers: 
        #   including the disturbance estimate in the controller 
        #   helps improve performance
        F_tilde = -self.K_lon @ xhat_lon - self.ki_lon * self.integrator_h
        F = self.Fe + F_tilde[0] - dhat_lon
        F_sat = saturate(F[0], self.F_limit)
        self.integratorAntiWindup(F_sat, F, self.ki_lon, self.integrator_h)

        tau = -self.K_lat @ xhat_lat - self.ki_lat*self.integrator_z - dhat_lat
        tau_sat = saturate(tau[0], self.tau_limit)
        self.integratorAntiWindup(tau_sat, tau, self.ki_lat, self.integrator_z)

        u = np.array([[F_sat], [tau_sat]])
        self.F_prev = F_sat
        self.tau_prev = tau_sat
        return u, xhat_lat, xhat_lon, dhat_lat, dhat_lon

    def integratorAntiWindup(self, u_sat, u_unsat, ki, integrator):
        if ki != 0.0:
            integrator = integrator + P.Ts/ki*(u_sat-u_unsat)

    def update_lat_observer(self, y_m):
        # update the observer using RK4 integration
        F1 = self.observer_f_lat(self.observer_state_lat, y_m)
        F2 = self.observer_f_lat(self.observer_state_lat + self.Ts / 2 * F1, y_m)
        F3 = self.observer_f_lat(self.observer_state_lat + self.Ts / 2 * F2, y_m)
        F4 = self.observer_f_lat(self.observer_state_lat + self.Ts * F3, y_m)
        self.observer_state_lat += self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
        xhat_lat = self.observer_state_lat[0:-1]
        dhat_lat = self.observer_state_lat[-1]
        return xhat_lat, dhat_lat

    def observer_f_lat(self, observer_state_lat, y_m):
        # observer_state_dot = 
        #           A2*(observer_state-equil_state) + B2*(u-ue) 
        #               + L2(y-C2*observer_state)
        observer_state_lat_dot = self.A2_lat @ observer_state_lat \
                + self.B2_lat * self.tau_prev \
                + self.L2_lat @ (y_m - self.C2_lat @ observer_state_lat)
        return observer_state_lat_dot

    def update_lon_observer(self, y_m):
        # update the observer using RK4 integration
        F1 = self.observer_f_lon(self.observer_state_lon, y_m)
        F2 = self.observer_f_lon(self.observer_state_lon + self.Ts / 2 * F1, y_m)
        F3 = self.observer_f_lon(self.observer_state_lon + self.Ts / 2 * F2, y_m)
        F4 = self.observer_f_lon(self.observer_state_lon + self.Ts * F3, y_m)
        self.observer_state_lon += self.Ts / 6 * (F1 + 2 * F2 + 2 * F3 + F4)
        xhat_lon = self.observer_state_lon[0:-1]
        dhat_lon = self.observer_state_lon[-1]
        return xhat_lon, dhat_lon

    def observer_f_lon(self, observer_state_lon, y_m):
        # observer_state_dot = 
        #           A2*(observer_state-equil_state) + B2*(u-ue) 
        #               + L2(y-C2*observer_state)
        xhat_dot = self.A2_lon @ observer_state_lon \
                   + self.B2_lon * (self.F_prev - self.Fe) \
                   + self.L2_lon @ (y_m - self.C2_lon @ observer_state_lon)
        return xhat_dot


def saturate(u, limit):
    if abs(u) > limit:
        u = limit * np.sign(u)
    return u


