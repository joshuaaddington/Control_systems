#%%
import sympy as sp
from sympy.physics.vector import dynamicsymbols
from sympy.physics.vector.printing import vlatex
from IPython.display import Math, display

from sympy import sin, cos, diff, Matrix, diag, symbols, Function, pretty_print, simplify, init_printing, latex

###################################################################################
# find Kinetic Energy (with vectors/matrices) for case study E blockbeam
###################################################################################
# importing these functions directly to make life a little easier and code a little
from sympy import sin, cos, Matrix, diag, symbols, simplify

#defining mathematical variables (called symbols in sp) and time varying functions
t, m1, m2, ell, g = symbols('t, m1, m2, ell, g')
z = dynamicsymbols('z')
theta = dynamicsymbols('theta')

#defining generalized coords and derivatives
q = Matrix([[z], [theta]])
qdot = q.diff(t)

#defining the kinetic energy
p1 = Matrix([[z*cos(theta)], [z*sin(theta)], [0]])
p2 = Matrix([[ell/2*cos(theta)], [ell/2*sin(theta)], [0]])
v1 = p1.diff(t)
v2 = p2.diff(t)

# by inspection, we can find the angular velocity of the rod
omega = Matrix([[0], [0], [theta.diff(t)]])

# if we are uncomfortable with the above, we can use the rotation matrix describing
# the orientation of the beam relative to the inertial frame to find the angular ve

R = Matrix([[cos(theta), -sin(theta), 0], [sin(theta), cos(theta), 0], [0, 0, 1]])
# omega = Omeg_mat.vee()
# next we define the inertia tensor for the beam, modeled as a thin rod
J = diag(0, m2*ell**2/12.0, m2*ell**2/12.0)
Ke = simplify(0.5*m1*v1.T @ v1 + 0.5*m2*v2.T @ v2 + 0.5*omega.T @ R @ J @ R.T @ omega)

# just grabbing the scalar inside this matrix so that we can do L = K-P, since P is
Ke = Ke[0,0]

Pe = m1*g*p1[1] + m2*g*p2[1]

L = sp.simplify(Ke - Pe)


#%%
# Solution for Euler-Lagrange equations, but this does not include right-hand side (like friction and tau)
EL = simplify( diff(diff(L, qdot), t) - diff(L, q) )

#%% 
############################################################
### Including friction and generalized forces, then solving for highest order derivatives
############################################################

# these are just convenience variables
zd = z.diff(t)
zdd = zd.diff(t)
thetad = theta.diff(t)
thetadd = thetad.diff(t)

# defining symbols for external force and friction
F = symbols('F')

# defining the right-hand side of the equation and combining it with E-L part
RHS = Matrix([[0],[F*ell*cos(theta)]])
full_eom = EL - RHS

# finding and assigning zdd and thetadd
# if our eom were more complicated, we could rearrange, solve for the mass matrix, and invert it to move it to the other side and find qdd and thetadd
result = simplify(sp.solve(full_eom, (zdd, thetadd)))

#TODO - add an example of finding the same thing, but not using sp.solve


# result is a Python dictionary, we get to the entries we are interested in
# by using the name of the variable that we were solving for
zdd_eom = result[zdd]  # EOM for zdd, as a function of states and inputs
thetadd_eom = result[thetadd]  # EOM for thetadd, as a function of states and inputs
display(Math(vlatex(zdd_eom)))
display(Math(vlatex(thetadd_eom)))

#%% [markdown]
# OK, now we can get the state variable form of the equations of motion.

#%%
import blockbeamParam as P

# If we want to define fixed parameters that are not states or inputs we can here
# But I'm going to keep them as inputs to the equation for now just so that I can change stuff on the fly
params = [(g, P.g), (ell, P.ell)]

# substituting parameters into the equations of motion
# zdd_eom = zdd_eom.subs(params)
# thetadd_eom = thetadd_eom.subs(params)
# print("after subsituting parameters")
# display(Math(vlatex(zdd_eom)))
# display(Math(vlatex(thetadd_eom)))

# now defining the state variables that will be passed into f(x,u) 
state = sp.Matrix([z, theta, zd, thetad])
ctrl_input = sp.Matrix([F])

# defining the function that will be called to get the derivatives of the states
state_dot = sp.Matrix([zd, thetad, zdd_eom, thetadd_eom])
display(Math(vlatex(Matrix(state_dot))))


#%%
import numpy as np

# converting the function to a callable function that uses numpy to evaluate and 
# return a list of state derivatives
eom = sp.lambdify([state, ctrl_input, m1, m2], np.array(state_dot), 'numpy')

# calling the function as a test to see if it works:
cur_state = [0, 0, 0, 0]
cur_input = [1]
# print("x_dot = ", eom(cur_state, cur_input, P.m1, P.m2))

# %%
