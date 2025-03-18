#%%
from h3_generate_E_L import * 
from IPython.display import display, Math
from sympy.physics.vector import vlatex

#%%
# TODO - calculate the full equations of motion for the hummingbird, you can do this
# as a single equation = 0 (i.e. LHS - RHS), or you can treat each part individually and
# do the calculations for LHS and RHS separately. 
q_ddot = q_dot.diff(t)
B = beta*sp.eye(3)  # damping matrix
LHS = M*q_ddot + C + dP_dq
RHS = tau - B@q_dot

#%%[markdown]
# Longitudinal Dynamics are derived in this section: 

#%%
# step 0 - just print the equations for theta_ddot

# display(Math(vlatex(sp.expand_trig(LHS[1]))))
# display(Math(vlatex(sp.expand_trig(RHS[1]))))


#%%
# step 1 - substitute zero in for all the suggested variables
long_zero_subs = [(phi, 0), (q_dot[0], 0), (q_dot[2], 0), (q_ddot[2], 0), (beta, 0)]
LHS_subbed_long = LHS[1].subs(long_zero_subs)
RHS_subbed_long = RHS[1].subs(long_zero_subs)


# display('Longitudinal Dynamics (step 1):')
# display(Math(vlatex(LHS_subbed_long)))
# display(Math(vlatex(RHS_subbed_long)))

#%%

# step 2 - substitute F and Tau for f_l, f_r, and d*(f_l-f_r)
F = sp.symbols('F')

long_input_subs = [(f_l + f_r, F)]
LHS_subbed_long = LHS_subbed_long.subs(long_input_subs)
RHS_subbed_long = RHS_subbed_long.subs(long_input_subs)

FULL_SUBBED_LONG = LHS_subbed_long - RHS_subbed_long

# display('Longitudinal Dynamics (step 2):')
# display(Math(vlatex(FULL_SUBBED_LONG)))

#%%
# step 3 - find the feedback linearization force F_fl
F_ctrl = sp.symbols('F_ctrl')
F_fl = sp.symbols('F_fl')
F_fl = ((m1*ell_1 + m2*ell_2)*g)/ell_T
FULL_SUBBED_LONG = sp.simplify(FULL_SUBBED_LONG.subs(F, F_fl + F_ctrl))
FULL_SUBBED_LONG = FULL_SUBBED_LONG.subs(g, 9.810)
# display('Feedback Linearization Force (step 3):')
# display(Math(vlatex(F_fl)))

#%%
# step 4 - find the linearized equation of motion for theta_ddot
display('Linearized Longitudinal Dynamics (step 4):')
Thetadd_EOM = sp.simplify(sp.solve(FULL_SUBBED_LONG, q_ddot[1])[0])
Thetadd_EOM = Thetadd_EOM.subs([(J1x, P.J1x), (J1y, P.J1y), (J1z, P.J1z), (J2x, P.J2x), (J2y, P.J2y), (J2z, P.J2z), (J3x, P.J3x), (J3y, P.J3y), (J3z, P.J3z), (ell_1, P.ell1), (ell_2, P.ell2), (ell_T, P.ellT), (m1, P.m1), (m2, P.m2), (m3, P.m3), (theta, 0)])
# display(Math(vlatex(Thetadd_EOM)))

#%% [markdown]
# **Lateral Dynamics are derived in this section:**

#%%
# step 1 - substitute zero in for all the suggested variables
Tau = sp.symbols('Tau')
F = sp.symbols('F')
lat_subs = [(q[1], 0), (q_dot[1], 0), (q_ddot[1], 0), (f_l + f_r, F), (d*(f_l - f_r), Tau), (beta, 0)]
LHS_lat = LHS[:1, :]  # First row
LHS_lat = LHS_lat.col_join(LHS[2:, :])  # Join with the third row onward
RHS_lat = RHS[:1, :]  # First row
RHS_lat = RHS_lat.col_join(RHS[2:, :])  # Join with the third row onward

LHS_zeroes_subbed_lat = LHS_lat.subs(lat_subs)
RHS_zeroes_subbed_lat = RHS_lat.subs(lat_subs)
FULL_zeroes_subbed_lat = LHS_zeroes_subbed_lat - RHS_zeroes_subbed_lat

# display('Lateral Dynamics (step 1):')
# display(Math(vlatex(sp.simplify(FULL_zeroes_subbed_lat))))

#%%
# step 2 - find the equilibrium force F_e
F_e_lat = sp.symbols('F_e_lat')

#%%
# step 3 - find the equilibrium variables Tau_e and phi_e
Tau_e_lat = sp.symbols('Tau_e')
Phi_e_lat = sp.symbols('phi_e')
Tau_e_lat = 0
Phi_e_lat = 0

equilibrium_subs = [(F, F_e_lat), (Tau, Tau_e_lat), (phi, Phi_e_lat), (q_dot[0], 0), (q_dot[2], 0), (q_ddot[0], 0), ]

#%%
# step 4 - find the linearized equation of motion for phi_ddot and psi_ddot
phidd_eom = sp.solve(FULL_zeroes_subbed_lat, q_ddot[0])[q_ddot[0]]
psidd_eom = sp.solve(FULL_zeroes_subbed_lat, q_ddot[2])[q_ddot[2]]

# display(Math(vlatex(sp.simplify(phidd_eom))))
# display(Math(vlatex(sp.simplify(psidd_eom))))

x = sp.Matrix([q[0], q[2], q_dot[0], q_dot[2]])
x_dot = sp.Matrix([q_dot[0], q_dot[2], phidd_eom, psidd_eom])
u = sp.Matrix([Tau])

A = sp.simplify(x_dot.jacobian(x).subs(equilibrium_subs))
B = sp.simplify(x_dot.jacobian(u).subs(equilibrium_subs))

display('Linearized Lateral Dynamics (step 4):')
display(Math(vlatex(sp.simplify(A))))
display(Math(vlatex(sp.simplify(B))))
#%%
# Substitute numerical values for the parameters
A_num = sp.simplify(A.subs([(ell_1, P.ell1), (ell_2, P.ell2), (ell_T, P.ellT), (m1, P.m1), (ell_3x, P.ell3x), (ell_3y, P.ell3y), (ell_3z, P.ell3z), (m2, P.m2), (m3, P.m3), (g, P.g), (J1x, P.J1x), (J1y, P.J1y), (J1z, P.J1z), (J2x, P.J2x), (J2y, P.J2y), (J2z, P.J2z), (J3x, P.J3x), (J3y, P.J3y), (J3z, P.J3z)]))
B_num = sp.simplify(B.subs([(ell_1, P.ell1), (ell_2, P.ell2), (ell_T, P.ellT), (m1, P.m1), (ell_3x, P.ell3x), (ell_3y, P.ell3y), (ell_3z, P.ell3z), (m2, P.m2), (m3, P.m3), (g, P.g), (J1x, P.J1x), (J1y, P.J1y), (J1z, P.J1z), (J2x, P.J2x), (J2y, P.J2y), (J2z, P.J2z), (J3x, P.J3x), (J3y, P.J3y), (J3z, P.J3z)]))

# display('Numerical Linearized Lateral Dynamics (step 4):')
# display(Math(vlatex(sp.simplify(A_num))))
# display(Math(vlatex(sp.simplify(B_num))))

# %%
