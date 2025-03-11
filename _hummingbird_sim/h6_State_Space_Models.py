#%%
import sympy as sp
from IPython.display import display, Math
from sympy.physics.vector import vlatex
from h4_linearize_dyn import *
from h5_transfer_functions import *

# Longitudinal State Space Model
x_long = sp.Matrix([q[1], q_dot[1]])
x_long_dot = sp.Matrix([q_ddot[1], Thetadd_EOM])
u_long = sp.Matrix([F_ctrl])

A_long = x_long_dot.jacobian(x_long)
B_long = x_long_dot.jacobian(u_long)
C_long = sp.Matrix([[1, 0]])
D_long = sp.Matrix([[0]])

# display('Longitudinal State Space Model:')
# display("A:")
# display(Math(vlatex(A_long)))
# display("X:")
# display(Math(vlatex(x_long)))
# display("B:")
# display(Math(vlatex(B_long)))
# display("U:")
# display(Math(vlatex(u_long)))

#%%
# Lateral State Space Model
x_lat = sp.Matrix([q[0], q[2], q_dot[0], q_dot[2]])
x_dot_lat = sp.Matrix([q_dot[0], q_dot[2], phidd_eom, psidd_eom])
u_lat = sp.Matrix([Tau])

A_lat= sp.simplify(x_dot.jacobian(x).subs(equilibrium_subs).subs([(ell_1, P.ell1), (ell_2, P.ell2), (ell_T, P.ellT), (m1, P.m1), (ell_3x, P.ell3x), (ell_3y, P.ell3y), (ell_3z, P.ell3z), (m2, P.m2), (m3, P.m3), (g, P.g), (J1x, P.J1x), (J1y, P.J1y), (J1z, P.J1z), (J2x, P.J2x), (J2y, P.J2y), (J2z, P.J2z), (J3x, P.J3x), (J3y, P.J3y), (J3z, P.J3z)]))
B_lat = sp.simplify(x_dot.jacobian(u).subs(equilibrium_subs).subs([(ell_1, P.ell1), (ell_2, P.ell2), (ell_T, P.ellT), (m1, P.m1), (ell_3x, P.ell3x), (ell_3y, P.ell3y), (ell_3z, P.ell3z), (m2, P.m2), (m3, P.m3), (g, P.g), (J1x, P.J1x), (J1y, P.J1y), (J1z, P.J1z), (J2x, P.J2x), (J2y, P.J2y), (J2z, P.J2z), (J3x, P.J3x), (J3y, P.J3y), (J3z, P.J3z)]))
C_lat = sp.Matrix([[1, 0]])
D_lat = sp.Matrix([[0]])

# display('Lateral State Space Model:')
# display("A:")
# display(Math(vlatex(A_lat)))
# display("X:")
# display(Math(vlatex(x_lat)))
# display("B:")
# display(Math(vlatex(B_lat)))
# display("U:")
# display(Math(vlatex(u_lat)))

# %%
