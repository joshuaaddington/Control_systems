import matplotlib.pyplot as plt
import numpy as np
import pendulumParam as P
from signalGenerator import signalGenerator
from pendulumAnimation import pendulumAnimation
from dataPlotter import dataPlotter
from pendulumDynamics import pendulumDynamics
from ctrlPD import ctrlPD

# instantiate pendulum, controller, and reference classes
pendulum = pendulumDynamics(alpha=0.1)
controller = ctrlPD()
reference = signalGenerator(amplitude=0.5, frequency=0.04)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = pendulumAnimation()

t = P.t_start  # time starts at t_start
y = pendulum.h()  # output of system at start of simulation
# add a random disturbance to the system
# disturbance = signalGenerator(amplitude=10, frequency=0.1)

# for part e), we can uncomment below
#pendulum.state[1,0] = 10.0*np.pi/180.0
#reference = signalGenerator(amplitude = 0.0, frequency=0.0)


while t < P.t_end:  # main simulation loop
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot

    while t < t_next_plot:
        r = reference.square(t)  # reference input
        # d = disturbance.sawtooth(t)
        x = pendulum.state  # use state instead of output
        u = controller.update(r, x)  # update controller
        y = pendulum.update(u)  # propagate system
        t += P.Ts  # advance time by Ts

    # update animation and data plots
    animation.update(pendulum.state)
    dataPlot.update(t, pendulum.state, u, r)
    plt.pause(0.0001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
