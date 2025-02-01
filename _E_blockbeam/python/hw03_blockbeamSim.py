import matplotlib.pyplot as plt
import blockbeamParam as P
from signalGenerator import signalGenerator
from blockbeamAnimation import blockbeamAnimation
from dataPlotter import dataPlotter
from blockbeamDynamics import blockbeamDynamics

# instantiate blockbeam, controller, and reference classes
blockbeam = blockbeamDynamics()
force = signalGenerator(amplitude=.5, frequency=1, y_offset = 11.5)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
animation = blockbeamAnimation()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    # updates control and dynamics at faster simulation rate
    while t < t_next_plot:  
        # Get referenced inputs from signal generators
        u = force.sin(t)        
        y = blockbeam.update(u)  # Propagate the dynamics
        t += P.Ts  # advance time by Ts
    
    # update animation and data plots
    animation.update(blockbeam.state)
    dataPlot.update(t, blockbeam.state, u)

    # the pause causes the figure to be displayed during the simulation
    plt.pause(0.0001)  

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
