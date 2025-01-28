#%%
# Sympy solving for Kinetic Energy
import sympy as sp
from sympy.physics.vector import dynamicsymbols
from sympy.physics.vector.printing import vlatex
from IPython.display import Math, display

t, m1, m2 = sp.symbols('t, m')
z = dynamicsymbols('z')
theta = dynamicsymbols('theta')
q = sp.Matrix([[z]])

p = sp.Matrix([[z*sp.cos(theta)], [z*sp.sin(theta)], [0]])
v = p.diff(t)

R = sp.Matrix([[sp.cos(theta), -sp.sin(theta), 0], [sp.sin(theta), sp.cos(theta), 0], [0, 0, 1]])
Omeg_mat = sp.simplify(R.diff(t)*R.T)
omega = sp.Matrix([[Omeg_mat[2,1]], [Omeg_mat[0,2]], [Omeg_mat[1,0]]])

J = sp.diag(0, 0, m*z**2)

K = sp.simplify(0.5*m*v.T @ v + 0.5*omega.T @ R @ J @ R.T @ omega)
display(Math(vlatex(K)))

#%%

# # Animation for Homework 2
# import matplotlib.pyplot as plt
# import numpy as np
# import blockbeamParam as B
# from signalGenerator import signalGenerator
# from blockbeamAnimation import blockbeamAnimation
# from dataPlotter import dataPlotter

# animation = blockbeamAnimation()

# t = 0
# A = 5
# omega = 5
# Y = 3
# x = 0
# x_dot = 0

# while(t < 10):
#     x = A * np.sin(omega*t) + Y
#     x_dot = 0
#     state = np.array([[x], [x_dot]])
#     t = t + 0.01
#     plt.pause(0.001)

#     animation.update(state)

#%%
# Box on Beam Simulation
import matplotlib.pyplot as plt
import numpy as np
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter

# instantiate reference input classes
z_fakeValueGenerator = signalGenerator(amplitude=0.5, frequency=0.1)
theta_fakeValueGenerator = signalGenerator(amplitude=.25*np.pi, frequency=.5)
f_fakeValueGenerator = signalGenerator(amplitude=5, frequency=.5)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = blockbeamAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # set variables
    z = z_fakeValueGenerator.sin(t)
    theta = theta_fakeValueGenerator.sin(t)
    f = f_fakeValueGenerator.sawtooth(t)
    # update animation
    state = np.array([[z], [theta], [0.0], [0.0]])
    animation.update(state)
    dataPlot.update(t, state, f)
    # advance time by t_plot    
    t += P.t_plot  
    plt.pause(0.05)  # allow time for animation to draw

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()

# %%
