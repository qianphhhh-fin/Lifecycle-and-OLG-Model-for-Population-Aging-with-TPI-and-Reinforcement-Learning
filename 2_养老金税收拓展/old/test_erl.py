import gymnasium as gym
from elegantrl.train.config import Config
from elegantrl.agents.AgentSAC import AgentSAC
from elegantrl.agents.AgentSAC import AgentModSAC
from utils.env_util import make_vec_env
from utils.vec_normalize import VecNormalize
from gymnasium.envs.registration import register
from v3_pensionfund_norm import PensionFundEnv
from elegantrl import get_gym_env_args
register(
    id='pensionfund-v3-norm',                                # call it whatever you want
    entry_point='v3_pensionfund_norm:PensionFundEnv', # module_name:class_name
)

if __name__ == '__main__':
        
    env = gym.make


    env_args = {
        "id": "pensionfund-v3-norm",
        "env_name": "pensionfund-v3-norm",
        "num_envs": 1,
        "max_step": 1000,
        "state_dim": 4,
        "action_dim": 4,
        "if_discrete": False,
    }
    # get_gym_env_args(env=PensionFundEnv(), if_print=True)
    args = Config(AgentSAC, env_class=env, env_args=env_args)

    args.max_step = 1000
    args.random_seed = 5
    # args.net_dims = [64, 64]
    # args.batch_size = int(256)
    args.buffer_init_size = int(10_000) # warmup
    args.reward_scale = 22.26  # RewardRange: -1800 < -200 < -50 < 0
    args.gamma = 0.95
    args.target_step = args.max_step
    args.eval_per_step = int(5000) 
    args.eval_times = 200
    args.gpu_id = 0 # if you have GPU
    args.break_step = 5_000_000 # break training if 'total_step > break_step'
    args.eval_seed = 3687851522
    args.num_workers = 8 # rollout workers number pre GPU (adjust it to get high GPU usage)
    args.learning_rate = 3e-04

    args.if_remove = True # remove the cwd folder? (True, False, None:ask me)


    from elegantrl.train.run import train_agent

    train_agent(args)