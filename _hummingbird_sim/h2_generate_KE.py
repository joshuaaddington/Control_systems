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
J1 = sp.diag(0.000189, 0.001953, 0.001894)
J1x = 0.000189
J1y = 0.001953
J1z = 0.001894
J2 = sp.diag(0.000231, 0.003274, 0.003416)
J2x = 0.000231
J2y = 0.003274
J2z = 0.003416
J3 = sp.diag(0.0002222, 0.0001956, 0.000027)
J3x = 0.0002222
J3y = 0.0001956
J3z = 0.000027


# TODO calculate M using the masses and the V, W, R, and J matrices
M = sp.zeros(3,3)
M22 = -J1y*sp.sin(phi)**2 + J1y + J1z*sp.sin(phi)**2 + J2y + m1*ell_1**2 + m2*ell_2**2
M23 = (J1y - J1z)*sp.sin(phi)*sp.cos(phi)*sp.cos(theta)
M33 = -J1x*sp.cos(theta)**2 + J1x - J1y*sp.cos(phi)**2*sp.cos(theta)**2 + J1y*sp.cos(theta)**2 + J1z*sp.cos(phi)**2*sp.cos(theta)**2 - J2x*sp.cos(theta)**2 \
    + J2x + J2z*sp.cos(theta)**2 + J3z + ell_1**2*m1*sp.cos(theta)**2 + ell_2**2*m2*sp.cos(theta)**2 + ell_3x**2*m3 + ell_3y**2*m3

M = sp.Matrix([
    [J1x,                  0,   -J1x*sp.sin(theta)],
    [0,                  M22,   M23               ],
    [-J1x*sp.sin(theta), M23,   M33               ]
])


# simplifying functions and displaying the result
M = sp.simplify(M)
M = sp.trigsimp(M)
M = sp.expand_trig(M)
display(Math(vlatex(M)))
# %%
