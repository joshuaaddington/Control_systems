#%%
from hw04_massLinearization import *
import sympy as sp
from sympy import symbols

# Extract State Space Model
A = A_linear
B = B_linear

# Define Laplace Variable
s = symbols('s')

# Define Output Matrices
C = sp.Matrix([[1, 0]])
D = sp.Matrix([[0]])

# Define Transfer Function
H = C*(s* sp.eye(2) - A).inv()@B + D
# print("Transfer Function H(s)/ F(s) = ")
# display(Math(vlatex(H)))
# %%
