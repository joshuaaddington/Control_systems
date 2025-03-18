import numpy as np
import hummingbirdParam as P

class ctrlLonPD:
    def __init__(self):
        F_fl = ((P.m1*P.ell1 + P.m2*P.ell2)*9.810)/P.ellT
        self.b0 = (P.ellT/(P.m1*P.ell1**2 + P.m2*P.ell2**2 + P.J1y + P.J2y))
        
        #  tuning parameters
        #tr = 0.8 # part (a)
        tr = 1 # tuned for faster rise time before saturation.
        zeta = 0.707

        # desired natural frequency
        wn = 0.5 * np.pi / (tr * np.sqrt(1 - zeta**2))
        alpha1 = 2.0 * zeta * wn
        alpha0 = wn**2

        # compute PD gains
        self.kp = alpha0 / self.b0
        self.kd = alpha1 / self.b0
        print("kp:", self.kp, "kd:", self.kd)

        self.theta_d1 = 0
        self.theta_dot = 0
        self.error_theta_d1 = 0


    def update(self, r, y):
        theta_ref = r[0][0]
        theta = y[1][0]
        thetadot = P.beta * self.theta_dot + (1 - P.beta) * ((theta - self.theta_d1) / P.Ts)
        
        force_equilibrium = (P.m1*P.ell1 + P.m2*P.ell2)*P.g*np.cos(theta)/P.ellT
        force_control = self.kp * (theta_ref - theta) - self.kd * thetadot
        force = force_equilibrium + force_control
        force = saturate(force, 0, P.force_max)
        torque = 0
        # convert force and torque to pwm signals
        pwm = P.mixing @ np.array([[force], [torque]]) / P.km
        pwm = saturate(pwm, 0, 1) 

        self.theta_d1 = theta

        return pwm, np.array([[0], [theta_ref], [0]])


def saturate(u, low_limit, up_limit):
    if isinstance(u, float) is True:
        u = np.max((np.min((u, up_limit)), low_limit))
    else:
        for i in range(0, u.shape[0]):
            u[i][0] = np.max((np.min((u[i][0], up_limit)), low_limit))
    return u