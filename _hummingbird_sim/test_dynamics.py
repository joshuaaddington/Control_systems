# use this file to test your vectors and matrices from the hummingbirdDynamics.py file
import numpy as np
from hummingbirdDynamics import HummingbirdDynamics as dynamics
import hummingbirdParam as P
import pickle as pkl

data = pkl.load(open("./test_matrices.pkl", "rb"))

# states are defined in the following order: [phi, theta, psi, phi_dot, theta_dot, psi_dot]
states = [np.array([0, 0, 0, 0, 0, 0]).reshape(6, 1),
          np.array([0.1, 0.1, -0.1, 0.1, -0.1, 0.1]).reshape(6, 1),
          np.array([0.05, -0.3, 0.22, 0.2, 0.1, -0.1]).reshape(6, 1)]

# inputs are defined in the following order: [f_l, f_r]
inputs = [np.array([[0], [0]]),
          np.array([[0.1], [0.1]]),
          np.array([[0.05], [-0.05]])]

hb_dynamics = dynamics(alpha=0.0)
np.set_printoptions(formatter={'float': lambda x: f"{x:14.3g}"})

for i in range(len(states)):
    M_test = hb_dynamics.M(states[i], hb_dynamics.param_vals)
    C_test = hb_dynamics.C(states[i], hb_dynamics.param_vals)
    dP_dq_test = hb_dynamics.dP_dq(states[i], hb_dynamics.param_vals)
    tau_test = hb_dynamics.tau(states[i], inputs[i], hb_dynamics.param_vals)

    print("Test case index", i)
    print("Your 'M' matrix is:\n", M_test)    
    print("Key 'M' matrix is:\n", data['M'][i], "\n\n")    

    print("Your 'C' matrix is:\n", C_test)
    print("Key 'C' matrix is:\n", data['C'][i], "\n\n")    

    print("Your 'dP_dq' matrix is:\n", dP_dq_test)
    print("Key 'dP_dq' matrix is:\n", data['dP_dq'][i], "\n\n")    

    print("Your 'tau' matrix is:\n", tau_test)
    print("Key 'tau' matrix is:\n", data['tau'][i], "\n\n")    

    input("\nPress Enter to continue...")
