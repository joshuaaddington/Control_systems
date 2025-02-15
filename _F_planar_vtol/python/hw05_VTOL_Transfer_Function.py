#%%
from hw4_VTOL_linearization import *
import sympy as sp

# Longitudinal Dynamics
# Extract State Space Model
A_lon = A_lin_lon
B_lon = B_lin_lon

# Define Laplace Variable
s = sp.symbols('s')

# Define Output Matrices
C_lon = sp.Matrix([[1, 0]])
D_lon = sp.Matrix([[0]])

# Define Transfer Function
H_lon = C_lon @ (s* sp.eye(2) - A_lon).inv() @B_lon + D_lon
print("Longitudinal Transfer Function H(s)/F(s) = ")
display(Math(vlatex(H_lon)))
#%%
# Lateral Dynamics
# Extract State Space Model
A_lat = A_lin_lat
B_lat = B_lin_lat

# Define Output Matrices
C_lat = sp.Matrix([[1, 0, 0, 0],[0, 1, 0, 0]])
D_lat = sp.Matrix([[0],[0]])

# Define Transfer Function
H_lat = C_lat @ (s* sp.eye(4) - A_lat).inv() @B_lat + D_lat
print("Lateral Transfer Function H(s)/Tau(s) = ")
display(Math(vlatex(H_lat)))
ZoverTheta = sp.simplify(H_lat[0]/H_lat[1])
print("Transfer Function Z(s)/Theta(s) = ")
display(Math(vlatex(ZoverTheta)))
# %%
