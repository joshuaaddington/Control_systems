import matplotlib.pyplot as plt
import numpy as np
import massParam as M
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter

#%%
# Simulation and animation of mass-spring-damper system
sine_plot = signalGenerator(amplitude=10.0, frequency=10, y_offset = 5)

# simulation plots and animation
dataplot = dataPlotter()
animation = massAnimation()

t = M.t_start  # time starts at t_start
while t < M.t_end:  # main simulation loop
    # set variables
    x = sine_plot.sin(t)
    F = sine_plot.sin(t)

    # update animation
    x_dot = 0.0
    state = np.array([[x], [x_dot]])  #state is made of x, and x_dot
    animation.update(state)
    dataplot.update(t, state, F)

    # advance time by t_plot
    t += M.t_plot  
    plt.pause(0.001)  # allow time for animation to draw

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()


#%%
# Animation for Homework 2
animation = massAnimation()

t = 0
A = 5
omega = 5
Y = 3
x = 0
x_dot = 0.0

while(t < 10):
    x = A * np.sin(omega*t) + Y
    x_dot = 0
    state = np.array([[x], [x_dot]])
    t = t + 0.01
    plt.pause(0.001)

    animation.update(state)