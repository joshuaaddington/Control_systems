import numpy as np
import massParam as P

class controllerPID:

    def __init__(self):
        self.limit = P.F_max
        #  tuning parameters
        self.ki = 3  # integrator gain


        # compute PD gains
        self.kp = 5
        self.kd = 2.92

        # print gains
        print('kp: ', self.kp)
        print('ki: ', self.ki)
        print('kd: ', self.kd)       
        
        # dirty derivative gain
        self.sigma = 0.05  

        #----------------------------------------------------------
        # variables for integrator and differentiator
        self.z_dot = 0  # estimated derivative of z
        self.z_dot_prev = 0
        self.z_prev = 0  # z delayed by one sample
        self.error_dot = 0.0  # estimated derivative of error
        self.error_prev = 0.0  # Error delayed by one sample
        self.integrator = 0.0  # integrator

    def update(self, z_r, y):
        z = y[0][0]

        # compute feedback linearized torque F_fl
        # F_fl = P.F_eq

        # compute the linearized torque using PID
        # Compute the current error
        error = z_r - z

        # differentiate z
        self.z_dot = (2.0*self.sigma - P.Ts) / (2.0*self.sigma + P.Ts) * self.z_dot_prev \
            + (2.0 / (2.0*self.sigma + P.Ts)) * ((z - self.z_prev))
        

        
        # Anti-windup scheme: only integrate z when z_dot is small
        if abs(self.z_dot < 0.1):
            self.integrator = self.integrator \
                + (P.Ts / 2) * (error + self.error_prev)
            
        F_e = P.F_e
        # PID control
        F_unsat = self.kp * error \
            + self.ki * self.integrator \
                - self.kd * self.z_dot + F_e
        
        # compute total torque
        F = self.saturate(F_unsat, P.F_max)

        # update delayed variables
        self.error_prev = error
        self.z_prev = z
        self.z_dot_prev = self.z_dot 
        return F

    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit*np.sign(u)
        return u







