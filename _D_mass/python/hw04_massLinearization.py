#%%
import sympy as sp
from sympy.physics.vector import dynamicsymbols
from sympy.physics.vector.printing import vlatex
from IPython.display import Math, display
from massStateEOM import *

from sympy import sin, cos, diff, Matrix, diag, symbols, Function, pretty_print, simplify, init_printing, latex

#%%
# display(Math(vlatex(state_dot)))

# Linearize the system
Fe = sp.Matrix([[0]])
Ze = sp.Matrix([[0]])
A = state_dot.jacobian(state)
B = state_dot.jacobian(ctrl_input)

A_linear = A.subs([(zd, 0), (z, Ze[0,0]), (F, Fe[0,0])])
B_linear = B.subs([(zd, 0), (z, Ze[0,0]), (F, Fe[0,0])])

print("A_linear = ")
display(Math(vlatex(A_linear)))
print("B_linear = ")
display(Math(vlatex(B_linear)))

# Feedback linearization


# %%
