#%%
import sympy as sp
from sympy.physics.vector import dynamicsymbols
from sympy.physics.vector.printing import vlatex
from IPython.display import Math, display

from sympy import sin, cos, diff, Matrix, diag, symbols, Function, pretty_print, simplify, init_printing, latex

# defining mathematical variables (called symbols in sp) and time varying functions like z and theta
t, m1, k = sp.symbols('t, m1, k')

# defining generalized coord
z = dynamicsymbols('z')
zd = z.diff(t)
zdd = zd.diff(t)
q = sp.Matrix([[z]])
qdot = q.diff(t)

# position of mass in inertial frame
p = sp.Matrix([[z], [0], [0]])
v = p.diff(t)

# calculate the kinetic energy
Ke = 0.5*m1*v.T @ v
Ke = Ke[0,0]

Pe = .5*k*z**2

L = sp.simplify(Ke - Pe)


#%%
# Solution for Euler-Lagrange equations, but this does not include right-hand side (like friction and tau)
EL = simplify( diff(diff(L, qdot), t) - diff(L, q) )

# display(Math(vlatex(EL)))



#%% 
############################################################
### Including friction and generalized forces, then solving for highest order derivatives
############################################################

# these are just convenience variables
zd = z.diff(t)
zdd = zd.diff(t)

# defining symbols for external force and friction
F, b = symbols('F, b')

# defining the right-hand side of the equation and combining it with E-L part
RHS = Matrix([[F - b*zd]])
full_eom = EL - RHS

# finding and assigning zdd and thetadd
# if our eom were more complicated, we could rearrange, solve for the mass matrix, and invert it to move it to the other side and find qdd and thetadd
result = simplify(sp.solve(full_eom, (zdd)))

#TODO - add an example of finding the same thing, but not using sp.solve


# result is a Python dictionary, we get to the entries we are interested in
# by using the name of the variable that we were solving for
zdd_eom = result[zdd]  # EOM for zdd, as a function of states and inputs

# display(Math(vlatex(zdd_eom)))


#%% [markdown]
# OK, now we can get the state variable form of the equations of motion.

#%%
import massParam as P

# If we want to define fixed parameters that are not states or inputs we can here
# But I'm going to keep them as inputs to the equation for now just so that I can change stuff on the fly
params = []

# substituting parameters into the equations of motion
zdd_eom = zdd_eom.subs(params)

# now defining the state variables that will be passed into f(x,u) 
state = sp.Matrix([[z], [zd]])
ctrl_input = sp.Matrix([[F]])

# defining the function that will be called to get the derivatives of the states
state_dot = sp.Matrix([[zd], [zdd_eom]])

display(Math(vlatex(state_dot)))

#%%
import numpy as np

# converting the function to a callable function that uses numpy to evaluate and 
# return a list of state derivatives
eom = sp.lambdify([state, ctrl_input, m1, k, b], np.array(state_dot), 'numpy')

# calling the function as a test to see if it works:
cur_state = [0, 0]
cur_input = [1]
print("x_dot = ", eom(cur_state, cur_input, P.m, P.k, P.b))

# %%
