import matplotlib.pyplot as plt
import numpy as np
import hummingbirdParam as P
from signalGenerator import SignalGenerator
from hummingbirdAnimation import HummingbirdAnimation
from dataPlotter import DataPlotter
from hummingbirdDynamics import HummingbirdDynamics
from _hummingbird_sim.ctrlLonPD import *

# instantiate reference input classes
phi_ref = SignalGenerator(amplitude=1.5, frequency=0.05)
theta_ref = SignalGenerator(amplitude=0.5, frequency=0.05)
psi_ref = SignalGenerator(amplitude=0.5, frequency=.05)

# instantiate the simulation plots and animation
dataPlot = DataPlotter()
animation = HummingbirdAnimation()
hummingbird = HummingbirdDynamics()
equilbrium = ctrlEquilibrium()

t = P.t_start  # time starts at t_start
while t < P.t_end:  # main simulation loop
    # set variables
    phi = 0
    theta = 0
    psi = 0

    # define dummy state and reference values (since we aren't simulating yet)
    ref = np.array([[0], [0], [0]])
    
    # convert force and torque to pwm values
    pwm = equilbrium.update(hummingbird.state)
    y = hummingbird.update(pwm)

    # update animation and data plots
    animation.update(t, hummingbird.state)
    dataPlot.update(t, hummingbird.state, pwm, ref)

    t = t + P.t_plot  # advance time by t_plot
    plt.pause(0.05)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()
