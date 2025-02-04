# Block on Beam Parameter File
import numpy as np
# import control as cnt

# Physical parameters of the  ballbeam known to the controller
m1 =  .35 # Mass of the block, kg
m2 =   2 # mass of beam, kg
ell =  .5 # length of beam, m
g = 9.8  # gravity at sea level, m/s^2

# parameters for animation
width = 0.05  # width of block
height = width*0.25  # height of block

# Initial Conditions
z0 = ell/2.0 # initial block position,m
theta0 = 0  # initial beam angle,rads
zdot0 =  0 # initial speed of block along beam, m/s
thetadot0 = 0 # initial angular speed of the beam,rads/s

# Simulation Parameters
t_start =  0 # Start time of simulation
t_end =  75 # End time of simulation
Ts =  .001 # sample time for simulation
t_plot =  .033 # the plotting and animation is updated at this rate

# saturation limits
F_max = 15  # Max Torque, N-m

# dirty derivative parameters
# sigma =   # cutoff freq for dirty derivative

# equilibrium force when block is in center of beam
ze = ell/2.0
Fe = m1*g*z0/ell + m2*g/2.0
