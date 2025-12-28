import gymnasium as gym
from elegantrl.train.config import Config
from elegantrl.agents.AgentSAC import AgentSAC
from elegantrl.train.config import Arguments,build_env
from elegantrl.agents.AgentSAC import AgentModSAC
from gymnasium.envs.registration import register
from elegantrl.train.run import train_agent
from elegantrl.agents.AgentPPO import AgentPPO
register(
    id='pensionfund-v3-norm',                                # call it whatever you want
    entry_point='v3_pensionfund_norm:PensionFundEnv', # module_name:class_name
)

if __name__ == '__main__':
        
    env_args = {
        "env_name": "pensionfund-v3-norm",
        "num_envs": 16,
        "max_step": 200,
        "state_dim": 4,
        "action_dim": 4,
        "if_discrete": False,
        'if_build_vec_env': True,
    }
    # get_gym_env_args(env=gym.make('pensionfund-v3-norm'), if_print=True)
    env_class = gym.make
    args = Config(AgentSAC, env_class=env_class, env_args=env_args)

    # Arguments(build_env('pensionfund-v3-norm'), AgentSAC())

    # args.net_dims = [64, 64]
    args.random_seed = 8
    args.batch_size = int(256)
    args.buffer_init_size = int(10_000) # warmup
    args.reward_scale = 2 ** -4 # RewardRange: -1800 < -200 < -50 < 0
    args.gamma = 0.95
    args.target_step = args.max_step
    args.eval_per_step = int(5000) 
    args.eval_times = 200
    args.gpu_id = 0 # if you have GPU
    args.break_step = 1_000_000 # break training if 'total_step > break_step'
    args.num_workers = 2 # rollout workers number pre GPU (adjust it to get high GPU usage)
    args.learning_rate = 1e-03

    args.if_remove = True # remove the cwd folder? (True, False, None:ask me)

    # eval env
    eval_env_args = {
        "env_name": "pensionfund-v3-norm",
        "num_envs": 1,
        "max_step": 200,
        "state_dim": 4,
        "action_dim": 4,
        "if_discrete": False,
        'if_build_vec_env': False,
    }
    args.eval_seed = 3687851522
    args.eval_env_class = gym.make  # eval_env = eval_env_class(*eval_env_args)
    args.eval_env_args = eval_env_args  # eval_env = eval_env_class(*eval_env_args)


    

    train_agent(args, if_single_process=False)