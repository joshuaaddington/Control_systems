#%%
from hw5_blockbeam_transfer_functions import *
import blockbeamParam as P
import sympy as sp
from numpy import pi, sqrt

# Rearrange Transfer Function to isolate s**2 term
Z_over_Theta = transfer_func_partC[0]/transfer_func_partC[1]

params = [(m1, P.m1), (m2, P.m2), (ell, P.ell), (g, P.g), (ze, P.ze)]
Z_over_Theta = Z_over_Theta.subs(params)

poles = sp.solve(Z_over_Theta.as_numer_denom()[1], s)

print("Poles of the system are:")
display(Math(vlatex(poles)))
# %%
# Isolate squared term in the denominator (this is with everything substituted already)
Z_over_F_num = Z_over_Theta.as_numer_denom()[0]
Z_over_F_den = Z_over_Theta.as_numer_denom()[1]
squared_coeff = Z_over_F_den.coeff(s, 2)
Z_over_F_num = Z_over_F_num/squared_coeff
Z_over_F_den = Z_over_F_den/squared_coeff

a0 = Z_over_F_den.coeff(s, 0)
a1 = Z_over_F_den.coeff(s, 1)
b0 = Z_over_F_num.coeff(s, 0)

# Find Closed Loop Transfer Function based on given Kp and Kd
Kp, Kd = symbols('Kp Kd')
Closed_loop = 0.688* s**2 + (a1 + b0*Kd)*s + (a0 + b0*Kp)
# Find the roots of the closed loop transfer function
closed_root1, closed_root2 = sp.solve(Closed_loop, s)
closed_root1 = sp.simplify(closed_root1)
closed_root2 = sp.simplify(closed_root2)
print("Roots of the closed loop transfer function are:")
display(Math(vlatex(closed_root1)))
display(Math(vlatex(closed_root2)))
#%%

# Calculate Kp and Kd based on desired roots
des_root1 = -1
des_root2 = -1.5
pole1 = sp.Eq((-a1 + b0*Kd)/2 + sp.sqrt(((a1 + b0*Kd)/2)**2 - (a0 + b0*Kp)), des_root1)
pole2 = sp.Eq((-a1 + b0*Kd)/2 - sp.sqrt(((a1 + b0*Kd)/2)**2 - (a0 + b0*Kp)), des_root2)
Kp_from_root, Kd_from_root = sp.solve((pole1, pole2), (Kp, Kd))[0]

print("Kp and Kd for desired roots are:")
display(Math(vlatex(Kp_from_root)))
display(Math(vlatex(Kd_from_root)))

# %%
# Calculate Kp and Kd based on desired rise time and damping ratio
tr = 10 # tuned for faster rise time before saturation.
zeta = 0.707

# desired natural frequency
wn = 0.5 * pi / (tr * sqrt(1 - zeta**2))
alpha1 = 2.0 * zeta * wn
alpha0 = wn**2
b0 = P.g

# compute PD gains
kd = alpha1 / b0
kp = alpha0 / b0
print("Kp and Kd for desired rise time and damping ratio are:")
display(Math(vlatex(kp)))
display(Math(vlatex(kd)))

# %%