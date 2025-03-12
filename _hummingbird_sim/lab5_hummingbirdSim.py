# import matplotlib.pyplot as plt
# import numpy as np
# import hummingbirdParam as P
# from signalGenerator import SignalGenerator
# from hummingbirdAnimation import HummingbirdAnimation
# from dataPlotter import DataPlotter
# from hummingbirdDynamics import HummingbirdDynamics
# from ctrlLonPD import *

# # instantiate reference input classes
# phi_ref = SignalGenerator(amplitude=1.5, frequency=0.05)
# theta_ref = SignalGenerator(amplitude=0.5, frequency=0.05)
# psi_ref = SignalGenerator(amplitude=0.5, frequency=.05)

# # instantiate the simulation plots and animation
# dataPlot = DataPlotter()
# animation = HummingbirdAnimation()
# hummingbird = HummingbirdDynamics()
# lonPD = ctrlLonPD()

# t = P.t_start  # time starts at t_start
# while t < P.t_end:  # main simulation loop
#     # set variables
#     phi = 0
#     psi = 0

#     # define dummy state and reference values (since we aren't simulating yet)
#     ref = np.array([[0, 0, 0]])
#     theta_d = 0
    
#     # convert force and torque to pwm values
#     pwm = lonPD.update(hummingbird.state, theta_d)
#     y = hummingbird.update(pwm)

#     # update animation and data plots
#     animation.update(t, hummingbird.state)
#     dataPlot.update(t, hummingbird.state, pwm, ref)

#     t = t + P.t_plot  # advance time by t_plot
#     plt.pause(0.05)

# # Keeps the program from closing until the user presses a button.
# print('Press key to close')
# plt.waitforbuttonpress()
# plt.close()

# \
import matplotlib.pyplot as plt
import numpy as np
import hummingbirdParam as P
from signalGenerator import SignalGenerator
from hummingbirdAnimation import HummingbirdAnimation
from dataPlotter import DataPlotter
from hummingbirdDynamicsOldSol import HummingbirdDynamics
from ctrlLonPD import ctrlLonPD

# instantiate pendulum, controller, and reference classes
hummingbird = HummingbirdDynamics()
controller = ctrlLonPD()
theta_ref = SignalGenerator(amplitude=15.*np.pi/180., frequency=0.05)

# instantiate the simulation plots and animation
dataPlot = DataPlotter()
animation = HummingbirdAnimation()

t = P.t_start  # time starts at t_start
y = hummingbird.h()
while t < P.t_end:  # main simulation loop

    # Propagate dynamics at rate Ts
    t_next_plot = t + P.t_plot
    while t < t_next_plot:
        r = np.array([[theta_ref.square(t)]])
        pwm, y_ref = controller.update(r, y)
        y = hummingbird.update(pwm)  # Propagate the dynamics
        t += P.Ts  # advance time by Ts

    # update animation and data plots at rate t_plot
    animation.update(t, hummingbird.state)
    dataPlot.update(t, hummingbird.state, pwm, y_ref)

    # the pause causes figure to be displayed during simulation
    plt.pause(0.001)

# Keeps the program from closing until the user presses a button.
print('Press key to close')
plt.waitforbuttonpress()
plt.close()