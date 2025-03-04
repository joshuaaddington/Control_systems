#%%
import sympy as sp
from sympy import laplace_transform
from h4_linearize_dyn import *

#%%
# Laplace symbol
s = sp.symbols('s')
Thetadd_EOM = Thetadd_EOM.subs(theta, 0)
Theta_over_F = laplace_transform(Thetadd_EOM, t, s)
Theta_over_F = sp.simplify(Theta_over_F[0])
display(Math(vlatex(Theta_over_F)))
# %%
# Lateral Transfer Functions
A_lin_phi = A_num[[0, 2], :][:, [0, 2]]
B_lin_phi = B_num[[0, 2], :]
C = sp.Matrix([[1, 0]])
D = sp.Matrix([[0]])
phi_over_tau = sp.simplify(C@(s*sp.eye(2)-A_lin_phi).inv() @B_lin_phi+D)[0]
# %%
A_lin_psi = A_num[[1, 3], :][:, [1, 3]]
B_lin_psi = A_num[3,0]
C = sp.Matrix([[1,0]])
D = sp.Matrix([[0]])
psi_over_phi = sp.simplify(C@(s*sp.eye(2)-A_lin_psi).inv() *B_lin_psi)[0]
# %%
