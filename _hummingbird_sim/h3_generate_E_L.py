#%%
from h2_generate_KE import * # import everything from the previous file, including our calculations for "M"
from sympy import sin, cos # import the sine and cosine functions

# symbols we need to define the right hand side of the equations of motion (tau - B@q_dot) and
# for potential energy (gravity)
# TODO define symbols needed for the RHS of the equations of motion (g, beta, d, f_l, f_r, ell_T, km)

g, beta, d, f_l, f_r, ell_T, km = sp.symbols('g beta d f_l f_r ell_T km')


# TODO calculate the kinetic energy using the definition with the mass matrix "M" 
# and "q_dot"
K = 1.0/2.0 * q_dot.T @ M @ q_dot # calculate the kinetic energy

# then make sure to grab just the scalar part of the result
K = K[0,0]


# TODO now calculat the potential energy "P"
P = m1*g*ell_1*sp.sin(theta) + m2*g*ell_2*sp.sin(theta) + m3*g*ell_3z
# this method would also work. 
# g_vec = sp.Matrix([[0], [0], [-g]]) # This would be negative because the z-axis is pointing down
#                                     # in figure 2-3 of the hummingbird lab manual. 
# P = m1*g_vec.T@p1_in_w + m2*g_vec.T@p2_in_w + m3*g_vec.T@p3_in_w
# P = P[0,0]

# TODO now calculate and define tau, C, and dP/dq for the hummingbird based on the definitions in the lab manual
tau = sp.Matrix([[d*(f_l-f_r)], 
                 [ell_T*(f_l+f_r)*cos(phi)], 
                 [ell_T*(f_l+f_r)*cos(theta)*sin(phi)-d*(f_l-f_r)*sin(theta)]])

C = M.diff(t)@q_dot - 0.5*sp.Matrix([[q_dot.T@M.diff(q[0])], [q_dot.T@M.diff(q[1])], [q_dot.T@M.diff(q[2])]])@q_dot

dP_dq = P.diff(q)

#%% TODO run this code to verify that your calculations match the lab manual
# display(Math(vlatex(M)))
# display(Math(vlatex(C)))
# display(Math(vlatex(dP_dq)))
# display(Math(vlatex(tau)))

# if yours does not match, you can use some combinations of sp.trigexpand(sp.simiplify()) to help match the lab manual, but this is not required
# and if you can't easily see the comparison, it may be worth moving on to using the "test_dynamics.py" file to perform numerical comparisons. 


#%%
# TODO You do not need to modify anything below this line, but it is worth looking at the code to understand what is happening
import numpy as np
import hummingbirdParam as P

params_sub = [(g, P.g)]

M = M.subs(params_sub)
C = C.subs(params_sub)
dP_dq = dP_dq.subs(params_sub)
tau = tau.subs(params_sub)

# TODO now you can either generate functions to calculate M, C, tau, friction, and dP_dq using lambdify,
# or you can hand-code the results you have found above into the hummingbirdDynamics.py file.

state = np.array([phi, 
        theta, 
        psi, 
        q_dot[0], 
        q_dot[1], 
        q_dot[2]]).reshape(6,1)

params = [m1, m2, m3, \
          J1[0,0], J1[1,1], J1[2,2], \
          J2[0,0], J2[1,1], J2[2,2], \
          J3[0,0], J3[1,1], J3[2,2], \
          ell_1, ell_2, ell_3x, ell_3y, ell_3z, ell_T, \
          d]

M_num = sp.lambdify([state, params], np.array(M), modules='numpy')
C_num = sp.lambdify([state, params], np.array(C), modules='numpy')
dP_dq_num = sp.lambdify([state, params], np.array(dP_dq), modules='numpy')
tau_num = sp.lambdify([state, [f_l, f_r], params], np.array(tau), modules='numpy')

import dill
dill.settings['recurse'] = True
dill.dump(M_num, open("hb_M_func.pkl", "wb"))
dill.dump(C_num, open("hb_C_func.pkl", "wb"))
dill.dump(dP_dq_num, open("hb_dP_dq_func.pkl", "wb"))
dill.dump(tau_num, open("hb_tau_func.pkl", "wb"))


#%%
# TODO, if you want to check your functions, you can do so as follows:
import hummingbirdParam as P
param_vals = [P.m1, P.m2, P.m3,
              P.J1x, P.J1y, P.J1z, P.J2x, P.J2y, P.J2z, P.J3x, P.J3y, P.J3z,
              P.ell1, P.ell2, P.ell3x, P.ell3y, P.ell3z, P.ellT, P.d]

x = np.array([[0., 0., 0., 0., 0., 0.]]).T
u = np.array([[0.24468/2, 0.22468/2]]).T
M_val = M_num(x, param_vals)
C_val = C_num(x, param_vals)
dP_dq_val = dP_dq_num(x, param_vals)
tau_val = tau_num(x, u, param_vals)
B_val = 0.001*np.eye(3)

# %%
