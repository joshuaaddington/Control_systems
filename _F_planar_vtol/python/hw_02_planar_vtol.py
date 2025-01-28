import matplotlib.pyplot as plt
import numpy as np
import VTOLParam as P
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter

# instantiate reference input classes
theta_fakeValueGenerator = signalGenerator(amplitude=2.0*np.pi, frequency=0.1)
z_fakeValueGenerator = signalGenerator(amplitude=0.5, frequency=0.1)
h_fakeValueGenerator = signalGenerator(amplitude=5, frequency=.5)
Tau_fakeValueGenerator = signalGenerator(amplitude=5, frequency=.5)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = VTOLAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # set variables
    theta = theta_fakeValueGenerator.sin(t)
    z = z_fakeValueGenerator.sin(t)
    h = h_fakeValueGenerator.sawtooth(t)
    Tau = np.array([[Tau_fakeValueGenerator.square(t)], [Tau_fakeValueGenerator.square(t)]])
    # update animation
    state = np.array([[theta], [z], [h]] )
    animation.update(state)
    dataPlot.update(t, state, Tau)
    # advance time by t_plot
    t += P.t_plot
    plt.pause(0.01)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()

#%%
import sympy as sp

# Define symbols
z, h, theta, t = sp.symbols('z h theta t')

z_dot = sp.diff(z, t)
h_dot = sp.diff(h, t)
theta_dot = sp.diff(theta, t)



# %%
