#%%
# these are libraries and functions that will help us with the calculations for h2
import sympy as sp
from sympy.physics.vector import dynamicsymbols
from sympy.physics.vector.printing import vlatex
from IPython.display import Math, display


#%%
# these functions are defined in our helper functions file in the public repository
from EL_helper_functions import rotx, roty, rotz, calc_omega, find_coeffs

# Defining necessary symbols and variables to use in the calculations
t, ell_1, ell_2, ell_3x, ell_3y, ell_3z, m1, m2, m3= sp.symbols('t, ell_1, ell_2, ell_3x, ell_3y, ell_3z, m1, m2, m3')

theta = dynamicsymbols('theta') # we have to tell SymPy that these are functions of time
psi = dynamicsymbols('psi')
phi = dynamicsymbols('phi')

q = sp.Matrix([[phi], [theta], [psi]])
q_dot = q.diff(t)

# TODO define the position of each mass in body frame, then rotate
# it into the world or inertial frame. 
p1_in_b = sp.Matrix([ell_1, 0, 0])
p1_in_w = rotz(psi)@roty(theta)@rotx(phi)@p1_in_b

p2_in_2 = sp.Matrix([ell_2, 0, 0])
p2_in_w = rotz(psi)@roty(theta)@p2_in_2

p3_in_1 = sp.Matrix([ell_3x, ell_3y, ell_3z])
p3_in_w = rotz(psi)@p3_in_1


#%% 
# TODO, take the time derivative of the position vectors to get the linear velocity,
# then use the "find_coeffs" function to calculate the "V_i" matrices
v1_in_w = p1_in_w.diff(t)
V1 = find_coeffs(v1_in_w, q_dot)

v2_in_w = p2_in_w.diff(t)
V2 = find_coeffs(v2_in_w, q_dot)

v3_in_w = p3_in_w.diff(t)
V3 = find_coeffs(v3_in_w, q_dot)


# now calculate the rotation matrices for each rigid body
# TODO use the "rotx", "roty", and "rotz" functions to calculate the rotation matrices
# for each rigid body
R1 =  rotz(psi)@roty(theta)@rotx(phi)
R2 =  rotz(psi)@roty(theta)
R3 =  rotz(psi)

# we can use the rotation matrices to calculate the angular velocity of each rigid body
# TODO use the "calc_omega" function to calculate the angular velocity of each rigid body 
omega_1 = calc_omega(R1)
omega_2 = calc_omega(R2)
omega_3 = calc_omega(R3)

# TODO use the "find_coeffs" function to calculate the "W_i" matrices
W1 = find_coeffs(omega_1, q_dot)
W2 = find_coeffs(omega_2, q_dot)
W3 = find_coeffs(omega_3, q_dot)


#%%
# TODO define the inertia tensors for each rigid body
J1x = sp.symbols('J1x')
J1y = sp.symbols('J1y')
J1z = sp.symbols('J1z')
J2x = sp.symbols('J2x')
J2y = sp.symbols('J2y')
J2z = sp.symbols('J2z')
J3x = sp.symbols('J3x')
J3y = sp.symbols('J3y')
J3z = sp.symbols('J3z')

J1 = sp.diag(J1x, J1y, J1z)
J2 = sp.diag(J2x, J2y, J2z)
J3 = sp.diag(J3x, J3y, J3z)

# TODO calculate M using the masses and the V, W, R, and J matrices
M = sp.zeros(3,3)
M += m1*V1.T*V1 + W1.T*R1*J1*R1.T*W1
M += m2*V2.T*V2 + W2.T*R2*J2*R2.T*W2
M += m3*V3.T*V3 + W3.T*R3*J3*R3.T*W3
    


# simplifying functions and displaying the result
M = sp.trigsimp(M)
M = sp.simplify(M)
M = sp.expand_trig(M)
display(Math(vlatex(M)))
# %%
