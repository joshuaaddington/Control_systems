import numpy as np
import hummingbirdParam as P


class ctrlPD:
    def __init__(self):
        # tuning parameters
        tr_pitch = 0.6   # rise time for pitch
        zeta_pitch = .707 # damping ratio for pitch
        tr_yaw = 2 # rise time for yaw
        zeta_yaw = .707  # damping ratio for yaw
        M = 10  # bandwidth separation
        tr_roll = tr_yaw / M  # rise time for roll
        zeta_roll = 0.707 # damping ratio for roll

        # gain calculation
        wn_pitch = 2.2 / tr_pitch  # natural frequency for pitch
        wn_yaw = 2.2 / tr_yaw  # natural frequency for yaw
        wn_roll = 2.2 / tr_roll  # natural frequency for roll
        self.kp_pitch = wn_pitch**2 / P.b_theta  
        self.kd_pitch = 2 * zeta_pitch * wn_pitch / P.b_theta  
        self.kp_roll = P.J1x * wn_roll**2
        self.kd_roll = P.J1x * 2 * zeta_roll * wn_roll
        self.kp_yaw = wn_yaw**2 / P.b_psi  
        self.kd_yaw = 2 * zeta_yaw * wn_yaw / P.b_psi  

        # print gains to terminal
        print('kp_pitch: ', self.kp_pitch)
        print('kd_pitch: ', self.kd_pitch) 
        print('kp_roll: ', self.kp_roll)
        print('kd_roll: ', self.kd_roll) 
        print('kp_yaw: ', self.kp_yaw)
        print('kd_yaw: ', self.kd_yaw) 
        
        self.Ts = P.Ts     # sample rate of the controller
        self.sigma = 0.05  # cutoff freq for dirty derivative

        # delayed variables
        self.phi_prev = 0.
        self.phi_dot = 0.
        self.theta_prev = 0.
        self.theta_dot = 0.
        self.psi_prev = 0.
        self.psi_dot = 0.
        self.error_theta_prev = 0.  # pitch error delayed by 1
        self.error_psi_prev = 0.  # yaw error delayed by 1

    def update(self, r, y):
        theta_ref = r[0][0]
        psi_ref = r[1][0]
        phi = y[0][0]
        theta = y[1][0]
        psi = y[2][0]
        force_equilibrium = P.g * (P.m1 * P.ell1 + P.m2 * P.ell2) \
                            * np.cos(theta) / P.ellT
        # compute errors
        error_theta = theta_ref - theta
        error_psi = psi_ref - psi

        # update differentiators
        self.phi_dot = (2 * self.sigma - P.Ts) / (2 * self.sigma + P.Ts) * self.phi_dot \
                       + (2.0 / (2.0*self.sigma + P.Ts)) * (phi - self.phi_prev)
        self.theta_dot = (2 * self.sigma - P.Ts) / (2 * self.sigma + P.Ts) * self.theta_dot \
                       + (2.0 / (2.0*self.sigma + P.Ts)) * (theta - self.theta_prev)
        self.psi_dot = (2 * self.sigma - P.Ts) / (2 * self.sigma + P.Ts) * self.psi_dot \
                       + (2.0 / (2.0*self.sigma + P.Ts)) * (psi - self.psi_prev)

        # pitch control - if you choose to use the F_e only instead, that is acceptable too
        force_unsat = force_equilibrium \
                      + self.kp_pitch * error_theta \
                      - self.kd_pitch * self.theta_dot
        force = saturate(force_unsat, -P.force_max, P.force_max)

        # outer loop yaw control
        phi_ref_unsat = self.kp_yaw * error_psi \
                        - self.kd_yaw * self.psi_dot
        phi_ref = saturate(phi_ref_unsat, -np.pi/4, np.pi/4)

        # inner loop pitch control
        error_phi = phi_ref - phi
        torque_unsat = self.kp_roll * error_phi - self.kd_roll * self.phi_dot
        torque = saturate(torque_unsat, -P.torque_max, P.torque_max)

        # convert force and torque to pwm signals
        pwm = np.array([[force + torque / P.d],               # u_left
                      [force - torque / P.d]]) / (2 * P.km)   # r_right          
        pwm = saturate(pwm, 0, 1)
 
        # update all delayed variables
        self.phi_prev = phi
        self.theta_prev = theta
        self.psi_prev = psi
        self.error_theta_prev = error_theta
        self.error_psi_prev = error_psi

        # return pwm plus reference signals
        return pwm, np.array([[phi_ref], [theta_ref], [psi_ref]])

def saturate(u, low_limit, up_limit):
    if isinstance(u, float) is True:
        if u > up_limit:
            u = up_limit
        if u < low_limit:
            u = low_limit
    else:
        for i in range(0, u.shape[0]):
            if u[i][0] > up_limit:
                u[i][0] = up_limit
            if u[i][0] < low_limit:
                u[i][0] = low_limit
    return u




