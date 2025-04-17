import matplotlib.pyplot as plt
import numpy as np
import VTOLParam as P
from signalGenerator import signalGenerator
from VTOLAnimation import VTOLAnimation
from dataPlotter import dataPlotter
from VTOLDynamics import VTOLDynamics
from dataPlotterObserver import dataPlotterObserver
from ctrlLQRVTOL import ctrlLQRVTOL

# instantiate VTOL, controller, and reference classes
VTOL = VTOLDynamics(alpha=0.02)
controller = ctrlLQRVTOL()
z_reference = signalGenerator(amplitude=4.0, frequency=0.05, y_offset=5.0)
h_reference = signalGenerator(amplitude=3.0, frequency=0.06, y_offset=5.0)
F_disturbance = signalGenerator(amplitude=2.0)
tau_disturbance = signalGenerator(amplitude=0.5)
z_noise = signalGenerator(amplitude=0.1)
h_noise = signalGenerator(amplitude=0.1)
th_noise = signalGenerator(amplitude=0.001)

# instantiate the simulation plots and animation
dataPlot = dataPlotter()
dataPlotObserver = dataPlotterObserver()
animation = VTOLAnimation()

t = P.t_start  # time starts at t_start
y = VTOL.h()  # output of system at start of simulation
while t < P.t_end:  # main simulation loop
    # Propagate dynamics in between plot samples
    t_next_plot = t + P.t_plot
    while t < t_next_plot:  # updates control and dynamics at faster simulation rate
        h_ref = h_reference.square(t)
        z_ref = z_reference.square(t)
        r = np.array([[z_ref], [h_ref]])  # reference input
        n = np.array([[z_noise.random(t)],
                      [h_noise.random(t)],
                      [th_noise.random(t)]])
        d = np.array([[F_disturbance.step(t)],
                      [tau_disturbance.step(t)]])

        u, xhat_lat, xhat_lon, dhat_lat, dhat_lon = controller.update(r, y + n)  # update controller   
        y = VTOL.update(P.mixing @ (u + d))  # propagate system
        t += P.Ts  # advance time by Ts
        
    # update animation and data plots
    animation.update(VTOL.state, z_ref)
    dataPlot.update(t, VTOL.state, u, z_ref, h_ref)
    dataPlotObserver.update(t, VTOL.state, xhat_lat, xhat_lon, 
                            d[1,0], dhat_lat[0],
                            d[0,0], dhat_lon[0])
    plt.pause(0.0001)  # the pause causes the figure to be displayed during the simulation

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
