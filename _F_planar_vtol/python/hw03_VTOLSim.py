import matplotlib.pyplot as plt
import VTOLParam as P
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter
from VTOLDynamics import VTOLDynamics

# instantiate VTOL, controller, and reference classes
VTOL = VTOLDynamics()
force = signalGenerator(amplitude=.5, frequency=1, y_offset = 11.5)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = VTOLAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    # updates control and dynamics at faster simulation rate
    while t < t_next_plot:  
        # Get referenced inputs from signal generators
        u = force.sin(t)        
        y = VTOL.update(u)  # Propagate the dynamics
        t += P.Ts  # advance time by Ts
    
    # update animation and data plots
    animation.update(VTOL.state)
    dataPlot.update(t, VTOL.state, u)

    # the pause causes the figure to be displayed during the simulation
    plt.pause(0.0001)  

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
