# Block on Beam Parameter File
import numpy as np
# import control as cnt

# Physical parameters of the  ballbeam known to the controller
m1 = 0.35  # Mass of the ball, kg
m2 = 2  # mass of beam, kg
length = 0.5  # length of beam, m
g = 9.81  # gravity at sea level, m/s^2

# equilibrium force when block is in center of beam
ze = length/2.0  # equilibrium position of ball, m
xe = np.array([[ze], [0], [0], [0]])  # equilibrium state
Fe = m1*g*ze/length + m2*g/2.0

# State Space Equations
A = np.array([[0, 0, 1, 0],
              [0, 0, 0, 1],
              [0, -1*g, 0, 0],
              [(-3*g*m1)/(length**2*m2 + 3*m1*ze**2), 0, 0, 0]])
B = np.array([[0], [0], [0], [(3*length)/(length**2*m2 + 3*m1*ze**2)]])
C = np.array([[1, 0, 0, 0]])  # output matrix

# parameters for animation
width = 0.05  # width of block
height = width*0.25  # height of block

# Initial Conditions
z0 = length/2.0  # initial ball position,m
theta0 = 0.0*np.pi/180  # initial beam angle,rads
zdot0 = 0.0  # initial speed of ball along beam, m/s
thetadot0 = 0.0  # initial angular speed of theh beam,rads/s

# Simulation Parameters
t_start = 0.0  # Start time of simulation
t_end = 75.0  # End time of simulation
Ts = 0.01  # sample time for simulation
t_plot = 0.033  # the plotting and animation is updated at this rate

# saturation limits
F_max = 150  # Max Force, N
