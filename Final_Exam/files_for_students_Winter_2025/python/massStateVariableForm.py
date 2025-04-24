#%%
import sympy as sp
from sympy.physics.vector import dynamicsymbols
from sympy.physics.vector.printing import vpprint, vlatex
from IPython.display import Math, display

from sympy import sin, cos, diff, Matrix, diag, symbols, Function, pretty_print, simplify, init_printing, latex

#%%

F, t, m1, c, g, k1, k2, b = symbols('F, t, m1, c, g, k1, k2, b')

z = dynamicsymbols('z')
zdot = z.diff(t)
zddot = zdot.diff(t)

q = Matrix([[z]])
qdot = q.diff(t)

p1 = Matrix([[0], [z], [0]])

v1 = diff(p1, t)

K = simplify(0.5*m1*v1.T@v1)
K = K[0,0]
display(Math(vlatex(K)))
#%%
# Find PE

P_spring = 0.5 * k1 * z**2 + 0.25 * k2 * z**4
P = m1 * g * z + P_spring

L = simplify(K-P)

#%%
EL_mass = simplify(diff(diff(L, qdot), t) - diff(L, q) )
#%%

RHS = Matrix([[F - c*sp.sign(zdot)]])
full_eom = EL_mass - RHS

result = simplify(sp.solve(full_eom, (zddot)))

denom = m1
zddot_eom = (F - b*zdot - k1*z - k2*z**3 + ((1/sp.sqrt(2))*m1*g))/denom

display(Math(vlatex(zddot_eom)))
#%%