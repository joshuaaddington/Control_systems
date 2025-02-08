#%%
import sympy as sp
from sympy.physics.vector import dynamicsymbols
from sympy.physics.vector.printing import vlatex
from IPython.display import Math, display

from sympy import sin, cos, diff, Matrix, diag, symbols, Function, pretty_print, simplify, init_printing, latex

###################################################################################
# find Kinetic Energy (with vectors/matrices) for case study E VTOL
###################################################################################
# importing these functions directly to make life a little easier and code a little
from sympy import sin, cos, Matrix, diag, symbols, simplify

#defining mathematical variables (called symbols in sp) and time varying functions
t, mc, mr, Jc, g, d, mu = symbols('t, mc, mr, Jc, g, d, mu')
z = dynamicsymbols('z')
h = dynamicsymbols('h')
theta = dynamicsymbols('theta')

#defining generalized coords and derivatives
q = Matrix([[z], [h], [theta]])
qdot = q.diff(t)

#defining the kinetic energy

# by inspection, we can find the angular velocity of the rod

# if we are uncomfortable with the above, we can use the rotation matrix describing
# the orientation of the beam relative to the inertial frame to find the angular ve
# omega = Omeg_mat.vee()
# next we define the inertia tensor for the beam, modeled as a thin rod
# just grabbing the scalar inside this matrix so that we can do L = K-P, since P is


#%%
# Solution for Euler-Lagrange equations, but this does not include right-hand side (like friction and tau)

#%% 
############################################################
### Including friction and generalized forces, then solving for highest order derivatives
############################################################

# these are just convenience variables
zd = z.diff(t)
zdd = zd.diff(t)
hd = h.diff(t)
hdd = hd.diff(t)
thetad = theta.diff(t)
thetadd = thetad.diff(t)

# defining symbols for external force and friction
Fr, Fl = symbols('Fr, Fl')

# defining the right-hand side of the equation and combining it with E-L part

full_eom = sp.Matrix([
    [(-(Fr + Fl)*sin(theta) - mu*zd)/(mc + 2*mr)],
    [(Fr + Fl)*cos(theta)/(mc + 2*mr) -g],
    [d*(Fr - Fl)/(Jc + 2*mr*d**2)]
])

# finding and assigning zdd and thetadd
# if our eom were more complicated, we could rearrange, solve for the mass matrix, and invert it to move it to the other side and find qdd and thetadd
zdd_eom = full_eom[0]
hdd_eom = full_eom[1]
thetadd_eom = full_eom[2]


# result is a Python dictionary, we get to the entries we are interested in
# by using the name of the variable that we were solving for
# display(Math(vlatex(zdd_eom)))
# display(Math(vlatex(hdd_eom)))
# display(Math(vlatex(thetadd_eom)))

#%% [markdown]
# OK, now we can get the state variable form of the equations of motion.

#%%
import VTOLParam as P

# If we want to define fixed parameters that are not states or inputs we can here
# But I'm going to keep them as inputs to the equation for now just so that I can change stuff on the fly
# params = [(g, P.g), (ell, P.ell)]

# # substituting parameters into the equations of motion
# zdd_eom = zdd_eom.subs(params)
# thetadd_eom = thetadd_eom.subs(params)
# print("after subsituting parameters")
# display(Math(vlatex(zdd_eom)))
# display(Math(vlatex(thetadd_eom)))

# now defining the state variables that will be passed into f(x,u) 
state = sp.Matrix([z, h, theta, zd, hd, thetad])
ctrl_input = sp.Matrix([Fr, Fl])

# defining the function that will be called to get the derivatives of the states
state_dot = sp.Matrix([zd, hd, thetad, zdd_eom, hdd_eom, thetadd_eom])
# display(Math(vlatex(Matrix(state_dot))))


#%%
# import numpy as np

# # converting the function to a callable function that uses numpy to evaluate and 
# # return a list of state derivatives
# eom = sp.lambdify([state, ctrl_input, m1, m2], np.array(state_dot), 'numpy')

# # calling the function as a test to see if it works:
# cur_state = [0, 0, 0, 0]
# cur_input = [1]
# # print("x_dot = ", eom(cur_state, cur_input, P.m1, P.m2))

# # %%
