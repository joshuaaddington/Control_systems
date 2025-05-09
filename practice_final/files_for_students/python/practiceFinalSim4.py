import matplotlib.pyplot as plt
import numpy as np
import rodMassParam as P
from signalGenerator import signalGenerator
from rodMassAnimation import rodMassAnimation
from dataPlotter import dataPlotter
from rodMassDynamics import rodMassDynamics
from dataPlotterObserver import dataPlotterObserver
from ctrlDisturbanceObserver import ctrlDisturbanceObserver

# instantiate system, controller, and reference classes
rodMass = rodMassDynamics()
controller = ctrlDisturbanceObserver()
reference = signalGenerator(amplitude=20*np.pi/180.0, frequency=0.1)
disturbance = signalGenerator(amplitude=0.5)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
dataPlotObserver = dataPlotterObserver()
animation = rodMassAnimation()

t = P.t_start
y = rodMass.h()
while t < P.t_end:
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        r = reference.square(t)
        d = disturbance.step(t)
        n = 0.0  #noise.random(t)
        u, xhat, dhat = controller.update(r, y + n)
        y = rodMass.update(u + d)
        t = t + P.Ts
    # update animation and data plots
    animation.update(rodMass.state)
    dataPlot.update(t, r, rodMass.state, u)
    dataPlotObserver.update(t, rodMass.state, xhat, d, dhat)
    plt.pause(0.05)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
