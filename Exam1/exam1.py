#%%
# Problem 6
import sympy as sp
from sympy.physics.vector import dynamicsymbols
from sympy import Matrix, diff, sin, cos, simplify, diag

theta = dynamicsymbols("theta")
y = dynamicsymbols("y")
t = sp.symbols("t")
J, mb, r, = sp.symbols("J, mb, r")

q = sp.Matrix([[theta], [y]])
qdot = q.diff(t)

thetadot = theta.diff(t)
ydot = y.diff(t)

Lagr = (0.5*J*thetadot**2 + 0.5*mb*ydot**2 + r*thetadot*ydot - y*cos(theta))

LHS = simplify(diff(diff(Lagr, qdot),t) - diff(Lagr, q))

#%%
# Problem 10
V = sp.symbols("V")
i, z = dynamicsymbols("i, z")
L, R, m, k, g = sp.symbols("L, R, m, k, g")

idot = i.diff(t)
zdot = z.diff(t)
zddot = zdot.diff(t)

idot = (V - R*i)/L
zddot = ((k*(i/z)**2) - m*g)/m

u = sp.Matrix([[V]])

X = sp.Matrix([[i],[zdot],[z]])
X_dot = sp.Matrix([[idot],[zddot],[zdot]])

A = X_dot.jacobian(X)
b = X_dot.jacobian(u)

# %%
# Problem 12
s = sp.symbols("s")

A = sp.Matrix([[-10, 0, 0], [0, 0, 1], [0, -.5, -1.1]])
B = sp.Matrix([[10],[0],[.5]])

C = sp.Matrix([[1, 0, 0], [0, 1, 0]])
D = sp.Matrix([[0],[0]])

ans = simplify(C@(s*sp.eye(3)-A).inv()@B+D)
Y_over_theta = ans[1]/ans[0]
# %%
# Problem 13
A = sp.Matrix([[.4, .2],[1, 0]])
B = sp.Matrix([[.2],[0]])

C = sp.Matrix([[0, 1]])
D = sp.Matrix([[0]])

ans = simplify(C@(s*sp.eye(2)-A).inv()@B+D)
# %%
