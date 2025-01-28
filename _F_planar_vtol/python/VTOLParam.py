# VTOL Parameter File
import numpy as np

# Physical parameters of the  VTOL known to the controller
mc =  4# kg
mr =   2# kg
Jc =  5 # kg m^2
d =   1# m
mu =  2 # kg/s
g = 9.8  # m/s^2
F_wind = 0 # wind disturbance force is zero in initial homeworks

# parameters for animation
length = 10.0

# Initial Conditions
z0 =  5 # initial lateral position
h0 =  5 # initial altitude
theta0 = 0 # initial roll angle
zdot0 = .5  # initial lateral velocity
hdot0 =  1 # initial climb rate
thetadot0 = 0  # initial roll rate
target0 = 0

# Simulation Parameters
t_start = 0 # Start time of simulation
t_end = 10  # End time of simulation
Ts =  .01 # sample time for simulation
t_plot = .01 # the plotting and animation is updated at this rate

# saturation limits
F_max = 100  # Max Force, N

# dirty derivative parameters
# sigma =   # cutoff freq for dirty derivative

# equilibrium force
# Fe =

# mixing matrix
unmixing = np.array([[1.0, 1.0], [-d, d]]) # converts fl and fr (LR) to force and torque (FT)
mixing = np.linalg.inv(unmixing) # converts force and torque (FT) to fl and fr (LR) 

