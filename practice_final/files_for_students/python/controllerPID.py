import numpy as np
import rodMassParam as P

class controllerPID:

    def __init__(self):
        self.limit = P.tau_max
        #  tuning parameters
        self.ki = 1.5  # integrator gain


        # compute PD gains
        self.kp = 8.57
        self.kd = .317

        # print gains
        print('kp: ', self.kp)
        print('ki: ', self.ki)
        print('kd: ', self.kd)       
        
        # dirty derivative gain
        self.sigma = 0.05  

        #----------------------------------------------------------
        # variables for integrator and differentiator
        self.theta_dot = 0  # estimated derivative of theta
        self.theta_dot_prev = 0
        self.theta_prev = 0  # theta delayed by one sample
        self.error_dot = 0.0  # estimated derivative of error
        self.error_prev = 0.0  # Error delayed by one sample
        self.integrator = 0.0  # integrator

    def update(self, theta_r, y):
        theta = y[0][0]

        # compute feedback linearized torque tau_fl
        # tau_fl = P.tau_eq

        # compute the linearized torque using PID
        # Compute the current error
        error = theta_r - theta

        # differentiate theta
        self.theta_dot = (2.0*self.sigma - P.Ts) / (2.0*self.sigma + P.Ts) * self.theta_dot_prev \
            + (2.0 / (2.0*self.sigma + P.Ts)) * ((theta - self.theta_prev))
        

        
        # Anti-windup scheme: only integrate theta when theta_dot is small
        if abs(self.theta_dot < 0.08):
            self.integrator = self.integrator \
                + (P.Ts / 2) * (error + self.error_prev)
            
        # PID control
        tau_unsat = self.kp * error \
            + self.ki * self.integrator \
                - self.kd * self.theta_dot
        
        # compute total torque
        tau = self.saturate(tau_unsat, P.tau_max)

        # update delayed variables
        self.error_prev = error
        self.theta_prev = theta
        self.theta_dot_prev = self.theta_dot 
        return tau

    def saturate(self, u, limit):
        if abs(u) > limit:
            u = limit*np.sign(u)
        return u







