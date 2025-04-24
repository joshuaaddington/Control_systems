import numpy as np
import massParam as P

class controllerStateSpace:
    def __init__(self):
        pass

    def update(self, z_r, y):
        pass
        #return F, x_hat, d_hat

    def update_observer(self, y):
        pass
        #return observer_state

    def observer_f(self, x_hat, y):
        pass
        #return observer_state_dot

    def saturate(self,u):
        if abs(u) > self.limit:
            u = self.limit*np.sign(u)
        return u