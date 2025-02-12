import numpy as np
import gym
from gym import spaces
from pendulumDynamics import pendulumDynamics
import pendulumParam as P

class PendulumEnv(gym.Env):
    """Custom Environment for Inverted Pendulum on a Cart"""
    
    def __init__(self):
        super(PendulumEnv, self).__init__()

        # Initialize pendulum dynamics
        self.pendulum = pendulumDynamics()
        
        # Define action space (continuous force applied to the cart)
        self.action_space = spaces.Box(low=-P.F_max, high=P.F_max, shape=(1,), dtype=np.float32)

        # Define observation space: [z, theta, zdot, thetadot]
        high = np.array([np.inf, np.pi, np.inf, np.inf])  # Unbounded position and velocity
        self.observation_space = spaces.Box(low=-high, high=high, dtype=np.float32)

    def reset(self):
        """Resets the environment to the initial state"""
        self.pendulum.state = np.array([
            [P.z0],  
            [P.theta0],  
            [P.zdot0],  
            [P.thetadot0],  
        ])
        return self._get_obs()

    def step(self, action):
        """Takes an action and returns next state, reward, done flag"""
        u = np.clip(action, -P.F_max, P.F_max)  # Ensure action is within limits
        
        # Update pendulum dynamics
        self.pendulum.update(u)

        # Get new state
        obs = self._get_obs()
        
        # Define reward function (penalize deviation from upright position)
        theta = obs[1]
        reward = np.cos(theta)  # Encourages upright balance
        
        # Define termination condition
        done = abs(theta) > np.pi / 2  # Episode ends if pendulum falls

        return obs, reward, done, {}

    def render(self, mode='human'):
        """Placeholder for visualization (can integrate with PendulumAnimation)"""
        pass  # You can integrate this with PendulumAnimation if needed

    def _get_obs(self):
        """Helper function to extract state as a flat array"""
        state = self.pendulum.state.flatten()
        return np.array(state, dtype=np.float32)
