import gymnasium as gym
import numpy as np
from gymnasium.envs.registration import register
import os
import argparse
from sbx import SAC
# from stable_baselines3 import SAC
from utils.env_util import make_vec_env
from utils.AllocationEvalCallback import EvalCallback
# from stable_baselines3.common.vec_env import VecNormalize
from utils.vec_normalize import VecNormalize
from optax import cosine_decay_schedule, chain, adamw, GradientTransformation, constant_schedule, linear_schedule
import optax
from custom_sac_policy import CustomSACPolicy # <--- 新增导入

# 注册自定义环境
# 确保在代码开始时执行，以便 make_vec_env 可以找到它
# register(
#     id='cocco-v1',
#     entry_point='v1_cocco_env:CoccoEnv',
# )
register(
    id='cocco-v2',
    entry_point='v2_cocco_env:CoccoEnvV2',
)

def run(args):
    """
    主运行函数，执行包含奖励归一化的课程学习流程。
    """
    # --- 路径和模型名称设置 ---
    log_dir = './models/cocco_sac_curriculum/'
    base_model_name = 'cocco_sac_run'

    # 自动创建新的运行文件夹，方便实验管理
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    

    # 自动确定运行编号
    existing_runs = [d for d in os.listdir(log_dir) if d.startswith(base_model_name)]
    run_num = len(existing_runs)

        
    model_dir = os.path.join(log_dir, f"{base_model_name}_{run_num}")
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)
    print(f"模型和日志将保存到: {model_dir}")
    

    # --- 课程学习流程 ---

    # ==========================================================================
    # 阶段一: 预训练 - 环境更简单 (无随机死亡)
    # ==========================================================================        
    # 阶段一环境参数：关闭随机死亡
    env_kwargs_p1 = {'disable_random_death': True}
    
    # 1. 创建原始的 VecEnv 实例
    vec_env_p1_raw = make_vec_env("cocco-v2", n_envs=24, seed=42, env_kwargs=env_kwargs_p1)
    eval_env_p1_raw = make_vec_env("cocco-v2", n_envs=100, seed=123, env_kwargs=env_kwargs_p1)

    # 2. 使用 VecNormalize 封装环境以进行奖励归一化
    #    MODIFIED: norm_obs_index_mask=[1] 表示不归一化观测向量的第2个元素 (年龄)
    #    MODIFIED: norm_reward=False, 与 MATLAB 行为一致，不归一化奖励
    vec_env_p1 = VecNormalize(vec_env_p1_raw, norm_obs=True, norm_obs_index_mask=[2], norm_reward=True)
    eval_env_p1 = VecNormalize(eval_env_p1_raw, norm_obs=True, norm_obs_index_mask=[2], norm_reward=False)
    
    # 3. 预热 VecNormalize
    burn_in_steps = 2000
    print(f"\n开始对 VecNormalize 进行 {burn_in_steps} 步的预热...")
    vec_env_p1.reset()
    for _ in range(burn_in_steps):
        random_actions = np.array([vec_env_p1.action_space.sample() for _ in range(vec_env_p1.num_envs)])
        vec_env_p1.step(random_actions)
    print("VecNormalize 预热完成，统计数据已初步建立。")

    # 从环境中获取真实的折现因子
    model_beta = vec_env_p1.get_attr("beta")[0]

    # 3. 设置评估回调函数
    eval_callback_p1 = EvalCallback(
        eval_env_p1,
        best_model_save_path=os.path.join(model_dir, "best_model"),
        log_path=os.path.join(model_dir, "best_model"),
        eval_freq=1000, 
        n_eval_episodes=300,
        deterministic=True,
    )

        # 1. 定义学习率调度器
    initial_learning_rate = 1e-4
    final_learning_rate = 1e-6
    total_steps = args.timesteps_phase1  # 总训练步数

    # lr_schedule = constant_schedule(initial_learning_rate)
        # 余弦退火调度 (从 initial_learning_rate 衰减到 final_learning_rate)
    lr_schedule = cosine_decay_schedule(
        init_value=initial_learning_rate,
        decay_steps=total_steps,
        alpha=final_learning_rate / initial_learning_rate  # 最终学习率比例
    )

    # 4. 创建并配置 SAC 模型
    model = SAC(
        'MlpPolicy',
        # CustomSACPolicy, # <--- 核心修改: 使用我们的自定义多头策略类
        vec_env_p1, 
        verbose=1,
        gamma=model_beta,
        learning_rate=lr_schedule,
        # learning_rate=1e-04,
        buffer_size=2_000_000,
        batch_size=256,
        train_freq=(1, 'step'), 
        gradient_steps=1,
        learning_starts=10000,
                # target_entropy=0.01,
        # tau=0.01,
        # ent_coef=0.01,  # <--- 在这里设置一个固定的较小值
                # policy_kwargs={
                #     "net_arch": {"pi": [256], "qf": [512, 512]},
                #     "optimizer_class": optax.adamw,
                #     "optimizer_kwargs": {"weight_decay": 0.0001},
                #     "n_critics": 3,
                # },
        policy_kwargs={
        "net_arch": [256,256],
        # "head_net_arch": [128],          # <--- 新增: 每个分离头的网络结构
        # "optimizer_class": optax.adamw,
        # "optimizer_kwargs": {
        #         "weight_decay": 1e-5     # <-- 设置权重衰减系数
        #     },
        },
        tensorboard_log=f"./tf_logs/cocco_sac/run_{run_num}_phase1"
    )
    
    # 5. 开始训练
    model.learn(total_timesteps=args.timesteps_phase1, callback=eval_callback_p1, progress_bar=True)
    

    # 6. 保存最终产出
    final_model_path = os.path.join(model_dir, "final_model.zip")
    final_stats_path = os.path.join(model_dir, "vecnormalize_final.pkl")
    model.save(final_model_path)
    vec_env_p1.save(final_stats_path)
    print(f"最终模型已保存到 {final_model_path}")
    print(f"最终 VecNormalize 统计数据已保存到 {final_stats_path}")



def parse_args():
    """
    解析命令行参数。
    """
    parser = argparse.ArgumentParser("Cocco Model Training with Curriculum Learning and Reward Normalization")
    parser.add_argument("--skip-phase1", action="store_true", help="跳过阶段一，直接从已保存的phase1模型开始微调。")
    parser.add_argument("--timesteps-phase1", type=int, default=30_000_000, help="阶段一（预训练）的总时间步数。")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args)