#%%
import sympy as sp
from sympy.physics.vector import dynamicsymbols
from sympy.physics.vector.printing import vlatex
from IPython.display import Math, display
from blockbeamStateEOM import *

from sympy import sin, cos, diff, Matrix, diag, symbols, Function, pretty_print, simplify, init_printing, latex

#%%
# display(Math(vlatex(state_dot)))
# display(Math(vlatex(state)))

# Linearize the system
Theta_e = sp.Matrix([[0]])
# Ze = sp.Matrix([[z]])
Fe = sp.Matrix([[g*(m2*ell + 2*m1*z)/2*ell]])
A = state_dot.jacobian(state)
B = state_dot.jacobian(ctrl_input)

# When trying to figure this out, substituting z for 0 gets me really close. There's something funny going on with z and the way it's being handled
A_linear = sp.simplify(A.subs([(zd, 0), (z, 0), (F, Fe[0,0]), (thetad, 0), (theta, Theta_e[0,0])]))
B_linear = sp.simplify(B.subs([(zd, 0), (z, 0), (F, Fe[0,0]), (thetad, 0), (theta, Theta_e[0,0])]))

# print("A_linear = ")
# display(Math(vlatex(A_linear)))
# print("B_linear = ")
# display(Math(vlatex(B_linear)))

# Feedback linearization


# %%
