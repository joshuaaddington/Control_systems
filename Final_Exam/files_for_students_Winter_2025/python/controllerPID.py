import numpy as np
import massParam as P

class controllerPID:

    def __init__(self):
        pass

    def update(self, z_r, y):
        pass
        #return F

    def saturate(self, u):
        if abs(u) > self.limit:
            u = self.limit*np.sign(u)
        return u







