#%%
from hw4_blockbeam_linearization import *
import sympy as sp
from sympy import symbols, Matrix, laplace_transform
from sympy.physics.control.lti import TransferFunction

# Extract State Space Model
A = A_lin
B = B_lin

# Define Laplace Variable
s = symbols('s')

# Define Output Matrices
C = sp.Matrix([[1, 0, 0, 0],[0, 1, 0, 0]])
D = sp.Matrix([[0],[0]])

# Define Transfer Function
H = C*(s* sp.eye(4) - A).inv()@B + D
Zovertheta = H[0]/H[1]
print("Transfer Function H = ")
display(Math(vlatex(H)))
print("Transfer Function Z(s)/ Theta(s) = ")
display(Math(vlatex(Zovertheta)))
#%%