#%%
import sympy as sp
from sympy.physics.vector import dynamicsymbols
from sympy.physics.vector.printing import vlatex
from IPython.display import Math, display

# defining mathematical variables (called symbols in sp) and time varying functions like z and theta
t, m1, k = sp.symbols('t, m1, k')

# defining generalized coord
z = dynamicsymbols('z')
zd = z.diff(t)
zdd = zd.diff(t)
q = sp.Matrix([[z]])
qd = q.diff(t)

# position of mass in inertial frame
p = sp.Matrix([[z], [0], [0]])
v = p.diff(t)

# calculate the kinetic energy
Ke = 0.5*m1*v.T @ v
Ke = Ke[0,0]

Pe = .5*m1*k*z**2

L = sp.simplify(Ke - Pe)

MassEomLeft = sp.simplify(sp.diff(sp.diff(L, qd), t) - sp.diff(L, q) )

# Now include generalized forces and friction
F, b = sp.symbols('F, b')

MassEomRight = sp.Matrix([[F - b*zd]])
fullEom = MassEomLeft - MassEomRight

result = sp.simplify(sp.solve(fullEom, qd))
zdd_eom = result[qd[0]]

display(Math(vlatex(zdd_eom)))

#%%
# Now we can import the parameters and initial conditions from the parameter file
import massParam as P

params = [(m1, P.m), (k, P.k), (b, P.b)]

zdd_eom = zdd_eom.subs(params)
display(Math(vlatex(zdd_eom)))

state = [z], [zd]
ctrl = [F]

state_dot = sp.Matrix([[zd], [zdd_eom]])

#%%
# Convert to numpy function
import numpy as np

eom = sp.lambdify([state, ctrl], np.array(state_dot), "numpy")

cur_state = np.array([[0], [0]])
cur_input = np.array([[1]])
output = eom(cur_state, cur_input)

print("x_dot = ", output)
# %%
