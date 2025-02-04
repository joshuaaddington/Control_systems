import matplotlib.pyplot as plt
import numpy as np
import VTOLParam as P
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter
from VTOLDynamics import Dynamics

# instantiate VTOL, controller, and reference classes
VTOL = Dynamics()
force = signalGenerator(amplitude=0.5, frequency=1.0)
torque = signalGenerator(amplitude=0.001, frequency=1.0, y_offset=-0.01)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = VTOLAnimation()

t = P.t_start # time starts at t_start
while t < P.t_end: # main simulation loop
# Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot

while t < t_next_plot: # updates control and dynamics at faster simulation rate
    f = (P.mc+2.0*P.mr)*P.g + force.sin(t)
    tau = torque.sin(t)

# converting force and torque to motor thrusts/forces
    motor_thrusts = P.mixing @ np.array([[f], [tau]])
    VTOL.update(motor_thrusts) # Propagate the dynamics
    t += P.Ts # advance time by Ts

# update animation and data plots
    animation.update(VTOL.state)
    dataPlot.update(t, VTOL.state, motor_thrusts)
    plt.pause(0.0001) # the pause causes the figure to be displayed during the simulation

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()