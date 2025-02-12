from pendulumEnv import *
from stable_baselines3 import PPO


env = PendulumEnv()
obs = env.reset()
print("Initial Observation:", obs)

action = np.array([0.0])  # No force applied
obs, reward, done, _ = env.step(action)
print("Next Observation:", obs, "Reward:", reward, "Done:", done)
