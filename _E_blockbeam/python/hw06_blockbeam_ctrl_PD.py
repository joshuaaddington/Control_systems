import numpy as np
import blockbeamParam as P

class ctrlPD:
    def __init__(self):
        ####################################################
        #       PD Control: Time Design Strategy
        ####################################################
        # tuning parameters
        tr_th = .2          # Rise time for inner loop (theta)
        zeta_th = 0.707       # inner loop Damping Coefficient

        # saturation limits
        F_max = 5             		  # Max Force, N
        error_max = 1        		  # Max step size,m
        theta_max = 30.0 * np.pi / 180.0  # Max theta, rads

        #---------------------------------------------------
        #                    Inner Loop
        #---------------------------------------------------
        # parameters of the open loop transfer function
        b0_th = 3*P.ell / (P.ell**2*P.m2 + 3*P.m1*P.ell/2)
        a0_th = 0

        # coefficients for desired inner loop
        wn_th = wn = 0.5 * np.pi / (tr_th * np.sqrt(1 - zeta_th**2))     # Natural frequency

        # compute gains
        self.kp_th = (wn_th**2 + a0_th) / b0_th
        self.kd_th = (2.0 * zeta_th * wn_th) / b0_th
        DC_gain = b0_th * self.kp_th / (b0_th * self.kp_th + a0_th)

        #---------------------------------------------------
        #                    Outer Loop
        #---------------------------------------------------
        # coefficients for desired outer loop
        M = 10.0           # Time scale separation 
        zeta_z = 0.707     # outer loop Damping Coefficient
        tr_z = M * tr_th   # desired rise time, s
        wn_z = wn = 0.5 * np.pi / (tr_z * np.sqrt(1 - zeta_z**2))  # desired natural frequency

        # compute gains
        self.kp_z = -wn_z**2 / P.g
        self.kd_z = -2.0 * zeta_z * wn_z / P.g

        # print control gains to terminal        
        print('DC_gain', DC_gain)
        print('kp_th: ', self.kp_th)
        print('kd_th: ', self.kd_th)
        print('kp_z: ', self.kp_z)
        print('kd_z: ', self.kd_z)

    def update(self, z_r, state):
        z = state[0][0]
        theta = state[1][0]
        zdot = state[2][0]
        thetadot = state[3][0]

        # the reference angle for theta comes from the
        # outer loop PD control
        theta_r = self.kp_z * (z_r - z) - self.kd_z * zdot

        # the force applied to the cart comes from the
        # inner loop PD control
        F = self.kp_th * (theta_r - theta) - self.kd_th * thetadot
        F = saturate(F)

        return F
    
    def equilibrium(self, state):
        # compute the equilibrium force
        z = state[0][0]
        
        F_eq = P.m1*P.g*z/P.ell + P.m2*P.g/2.0

        return F_eq
def saturate(u):
    if abs(u) > P.F_max:
        u = np.sign(u) * P.F_max
    return u






