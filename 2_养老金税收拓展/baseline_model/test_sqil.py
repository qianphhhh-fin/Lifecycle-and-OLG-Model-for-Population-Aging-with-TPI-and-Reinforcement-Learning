import datasets
from imitation.data import huggingface_utils
from imitation.algorithms import sqil
from imitation.util.util import make_vec_env
import numpy as np
from stable_baselines3 import sac
import numpy as np
import gymnasium as gym
from stable_baselines3 import PPO
from stable_baselines3.common.evaluation import evaluate_policy
from stable_baselines3.ppo import MlpPolicy
from imitation.algorithms.adversarial.airl import AIRL
from imitation.data import rollout
from imitation.data.wrappers import RolloutInfoWrapper
from imitation.policies.serialize import load_policy
from imitation.rewards.reward_nets import BasicShapedRewardNet
from imitation.util.networks import RunningNorm
from imitation.util.util import make_vec_env

SEED = 42

venv = make_vec_env(
    "Pendulum-v1",
    rng=np.random.default_rng(seed=SEED),
)


# Download some expert trajectories from the HuggingFace Datasets Hub.
expert = load_policy("ppo", venv, path="models/ppo-Pendulum-v1.zip")

# reward, _ = evaluate_policy(expert, env, 10)
# print(reward)

rng = np.random.default_rng()
expert_trajectories = rollout.rollout(
    expert,
    venv,
    rollout.make_sample_until(min_timesteps=None, min_episodes=50),
    rng=rng,
)


from imitation.data import rollout

trajectory_stats = rollout.rollout_stats(expert_trajectories)

print(
    f"We have {trajectory_stats['n_traj']} trajectories. "
    f"The average length of each trajectory is {trajectory_stats['len_mean']}. "
    f"The average return of each trajectory is {trajectory_stats['return_mean']}."
)



sqil_trainer = sqil.SQIL(
    venv=venv,
    demonstrations=expert_trajectories,
    policy="MlpPolicy",
    rl_algo_class=sac.SAC,
    rl_kwargs=dict(seed=SEED,tensorboard_log="../tf-logs/" ),
)

reward_before_training, _ = evaluate_policy(sqil_trainer.policy, venv, 100)
print(f"Reward before training: {reward_before_training}")

sqil_trainer.train(
    total_timesteps=100_000,
)  # Note: set to 300_000 to obtain good results
reward_after_training, _ = evaluate_policy(sqil_trainer.policy, venv, 100)
print(f"Reward after training: {reward_after_training}")