#%%
import sympy as sp
from h6_State_Space_Models import *
from h5_transfer_functions import *

# longitudinal controller class
class ctrl_longitudinal:
    def __init__(self):
        # Turn theta_over_F into a sympy fraction and pull off b0 and a1
        self.theta_over_F = theta_over_F
        self.theta_over_F_frac = sp.fraction(self.theta_over_F)
        self.theta_over_F_num = self.theta_over_F_frac[0]
        self.theta_over_F_den = self.theta_over_F_frac[1]
        self.b0 = self.theta_over_F_num.coeff(F_ctrl, 1)
        self.a1 = self.theta_over_F_den.coeff(s, 1)

        F = F_ctrl + F_fl

ctrl_test = ctrl_longitudinal()


# %%
