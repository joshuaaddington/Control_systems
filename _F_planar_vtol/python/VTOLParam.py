# VTOL Parameter File
import numpy as np

# Physical parameters of the  VTOL known to the controller
mc = 1.0  # kg
mr = 0.25  # kg
Jc = 0.0042  # kg m^2
d = 0.3  # m
mu = 0.1  # kg/s
g = 9.81  # m/s^2
F_wind = 0.0  # wind disturbance force is zero in initial homeworks

# parameters for animation
length = 10.0

# Initial Conditions
z0 = 0.0  # initial lateral position
h0 = 0.0  # initial altitude
theta0 = 0  # initial roll angle
zdot0 = 0  # initial lateral velocity
hdot0 = 0  # initial climb rate
thetadot0 = 0  # initial roll rate
target0 = 0

# Simulation Parameters
t_start = 0.0  # Start time of simulation
t_end = 50.0  # End time of simulation
Ts = 0.01  # sample time for simulation
t_plot = 0.033  # the plotting and animation is updated at this rate

# saturation limits
F_max = 10.0  # Max thrust produced by each motor, N
Tau_max = F_max * d  # Max torque produced by each motor, N*m

A_lon = np.array([
    [0, 1],
    [0, 0]
])

B_lon = np.array([
    [0],
    [1 / (mc + 2 * mr)]
])

C_lon = np.array([[1, 0]])  # output matrix

A_lat = np.array([
    [0, 0, 1, 0],
    [0, 0, 0, 1],
    [0, -g, -mu / (mc + 2*mr), 0],
    [0, 0, 0, 0]
])

B_lat = np.array([
    [0],
    [0],
    [0],
    [1 / (Jc + 2 * d**2 * mr)]
])

C_lat = np.array([[1, 0, 0, 0]])  # output matrix

# Print matrices to verify
print("A_lon:\n", A_lon)
print("B_lon:\n", B_lon)
print("A_lat:\n", A_lat)
print("B_lat:\n", B_lat)

# mixing matrix
unmixing = np.array([[1.0, 1.0], [d, -d]]) # converts fr and fl (RL) to force and torque (FT)
mixing = np.linalg.inv(unmixing) # converts force and torque (FT) to fr and fl (RL)

# equilibrium force
Fe = (mc + 2.0 * mr) * g
