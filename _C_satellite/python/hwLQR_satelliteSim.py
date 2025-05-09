import matplotlib.pyplot as plt
import numpy as np
import satelliteParam as P
from signalGenerator import signalGenerator
from satelliteAnimation import satelliteAnimation
from dataPlotter import dataPlotter
from satelliteDynamics import satelliteDynamics
from ctrlLQR import ctrlLQR
from dataPlotterObserver import dataPlotterObserver

# instantiate satellite, controller, and reference classes
satellite = satelliteDynamics(alpha=0.0)
controller = ctrlLQR()
reference = signalGenerator(amplitude=15.0*np.pi/180.0,
                            frequency=0.03)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
dataPlotObserver = dataPlotterObserver()
animation = satelliteAnimation()

t = P.t_start  # time starts at t_start
y = satellite.h()  # output of system at start of simulation

while t < P.t_end:  # main simulation loop
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot

    # updates control and dynamics at faster simulation rate
    while t < t_next_plot:  
        r = reference.square(t)  # reference input
        d = 0.0  # input disturbance
        u, xhat = controller.update(r, y)  # update controller
        y = satellite.update(u + d)  # propagate system
        t += P.Ts  # advance time by Ts

    # update animation and data plots
    animation.update(satellite.state)
    dataPlot.update(t, satellite.state, u, r)
    dataPlotObserver.update(t, satellite.state, xhat, d, 0.0)

    # the pause causes the figure to display during simulation
    plt.pause(0.0001)  

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
