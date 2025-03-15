
import matplotlib.pyplot as plt
import VTOLParam as P
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter
from VTOLDynamics import VTOLDynamics
from _F_planar_vtol.python.ctrlStateFeedback3 import ctrlStateFeedback

# instantiate VTOL, controller, and reference classes
VTOL = VTOLDynamics()
controller = ctrlStateFeedback()
Zref = signalGenerator(amplitude=0.2, frequency=0.04, y_offset=.25)
Href = signalGenerator(amplitude=0.2, frequency=0.04, y_offset=.25)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = VTOLAnimation()

t = P.t_start  # time starts at t_start
y = VTOL.h()  # output of system at start of simulation

while t < P.t_end:  # main simulation loop
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot

    # updates control and dynamics at faster simulation rate
    while t < t_next_plot: 
        r = [Href.sin(t), Zref.square(t)]   # reference input
        x = VTOL.state
        u = controller.update(r, x)  # update controller
        y = VTOL.update(u)  # propagate system
        t += P.Ts  # advance time by Ts
    # update animation and data plots
    animation.update(VTOL.state)
    dataPlot.update(t, VTOL.state, u, r)
    plt.pause(0.0001)  

# Keeps the program from closing until user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
