import matplotlib.pyplot as plt
import numpy as np
import massParam as P
from massDynamics import massDynamics
from ctrlObserver import ctrlObserver
from signalGenerator import signalGenerator
from massAnimation import massAnimation
from dataPlotter import dataPlotter
from dataPlotterObserver import dataPlotterObserver

# instantiate satellite, controller, and reference classes
mass = massDynamics()
controller = ctrlObserver()
reference = signalGenerator(amplitude=0.5, frequency=0.1)
disturbance = signalGenerator(amplitude=0.25)
noise = signalGenerator(amplitude=0.01)

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
        d = 0.0 #disturbance.step(t)  # input disturbance
        n = 0.0 # noise.random(t)  # simulate sensor noise
        u, xhat = controller.update(r, y + n)  # update controller
        y = mass.update(u + d)  # propagate system
        t += P.Ts  # advance time by Ts

    # update animation and data plots
    animation.update(mass.state)
    dataPlot.update(t, mass.state, u, r)
    dataPlotObserver.update(t, mass.state, xhat)
    plt.pause(0.0001)  # the pause causes the figure to be displayed during the simulation

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
