import matplotlib.pyplot as plt
import numpy as np
import massParam as P
from massDynamics import massDynamics
from ctrlDisturbanceObserver import ctrlDisturbanceObserver
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter
from dataPlotterObserver import dataPlotterObserver

# instantiate satellite, controller, and reference classes
mass = massDynamics(alpha=0.2)
controller = ctrlDisturbanceObserver()
reference = signalGenerator(amplitude=0.5, frequency=0.04)
disturbance = signalGenerator(amplitude=0.25)
noise = signalGenerator(amplitude=0.001)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = massAnimation()
dataPlotObserver = dataPlotterObserver()

t = P.t_start  # time starts at t_start
y = mass.h()  # output of system at start of simulation

while t < P.t_end:  # main simulation loop

    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    while t < t_next_plot:  # updates control and dynamics at faster simulation rate
        r = reference.square(t)  # reference input
        d = disturbance.step(t)  # input disturbance
        n = noise.random(t)  # simulate sensor noise
        u, x_hat, d_hat = controller.update(r, y + n)  # update controller
        y = mass.update(u + d)  # propagate system
        t += P.Ts  # advance time by Ts
        
    # update animation and data plots
    animation.update(mass.state)
    dataPlot.update(t, mass.state, u, r)
    dataPlotObserver.update(t, mass.state, x_hat, d, d_hat)
    plt.pause(0.001)  # the pause causes the figure to be displayed during the simulation

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
