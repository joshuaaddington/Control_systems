import numpy as np
import hummingbirdParam as P

class ctrlLatPD:
    def __init__(self):
        #  tuning parameters for roll
        tr_roll = .2
        zeta_roll = 0.9

        # desired natural frequency
        wn_roll = 0.5 * np.pi / (tr_roll * np.sqrt(1 - zeta_roll**2))
        alpha1_roll = 2.0 * zeta_roll * wn_roll
        alpha0_roll = wn_roll**2

        # PD gains for roll
        self.kp_roll = alpha0_roll * P.J1x
        self.kd_roll = alpha1_roll * P.J1x
        print("kp:", self.kp_roll, "kd:", self.kd_roll)

        self.phi_d1 = 0
        self.phi_dot = 0
        self.error_phi_d1 = 0

        """ YAW PARAMETERS"""
        # tuning parameters for yaw
        M = 10 # Separation factor
        tr_yaw = tr_roll * M
        zeta_yaw = 0.9

        wn_yaw = 0.5 * np.pi / (tr_yaw * np.sqrt(1 - zeta_yaw**2)) # Calculates the natural frequency of the desired system
        alpha1_yaw = 2 * zeta_yaw * wn_yaw
        alpha0_yaw = wn_yaw**2

        # PD gains for yaw
        b_psi =  (((P.ellT * P.m1*P.ell1 + P.m2*P.ell2)*9.810)/P.ellT) / (P.J1z + P.J2z + P.J3z + P.ell1**2*P.m1 + P.ell2**2*P.m2 + P.ell3x**2*P.m3 + P.ell3y**2*P.m3) # This will have to be changed with the real value
        self.kp_yaw = alpha0_yaw / b_psi
        self.kd_yaw = alpha1_yaw / b_psi

        self.psi_d1 = 0
        self.psi_dot = 0
        self.error_psi_d1 = 0


    def update(self, r, y):
        psi = y[1][0]
        phi = y[2][0]
        psidot = P.beta * self.psi_dot + (1 - P.beta) * ((psi - self.psi_d1) / P.Ts)
        phidot = P.beta * self.phi_dot + (1 - P.beta) * ((phi - self.phi_d1) / P.Ts)
        psi_ref = r[0][0]

        # the reference angle for theta comes from the
        # outer loop PD control
        phi_r = self.kp_yaw * (psi_ref - psi) - self.kd_yaw * psidot

        # Caluclate torque to achieve desired roll angle
        torque = self.kp_roll * (phi_r - phi) - self.kd_roll * phidot
        torque = saturate(torque, -P.torque_max, P.torque_max)
        force = ((P.m1*P.ell1 + P.m2*P.ell2)*9.810)/P.ellT
        # convert force and torque to pwm signals
        print("torque: ", torque)
        print("force: ", force)
        pwm = P.mixing @ np.array([[force], [torque]]) / P.km
        pwm = saturate(pwm, 0, 1)
        print("pwm: ", pwm)

        self.phi_d1 = phi
        self.psi_d1 = psi
        self.phi_dot = phidot
        self.psi_dot = psidot

        return pwm, np.array([[phi_r],[0], [psi_ref]])


def saturate(u, low_limit, up_limit):
    if isinstance(u, float) is True:
        u = np.max((np.min((u, up_limit)), low_limit))
    else:
        for i in range(0, u.shape[0]):
            u[i][0] = np.max((np.min((u[i][0], up_limit)), low_limit))
    return u