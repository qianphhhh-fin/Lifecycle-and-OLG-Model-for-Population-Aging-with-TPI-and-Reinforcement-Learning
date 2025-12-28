import gymnasium as gym
import numpy as np
from gymnasium.envs.registration import register
import os
import argparse
from sbx import SAC
from utils.env_util import make_vec_env
from utils.AllocationEvalCallback import EvalCallback
from utils.vec_normalize import VecNormalize
from stable_baselines3.common.logger import configure

# 导入必要的库
import optax
from sbx.sac.sac import SACState # SAC 使用一个自定义的 TrainState

# 注册自定义环境
try:
    register(
        id='cocco-v2',
        entry_point='v2_cocco_env:CoccoEnvV2',
    )
except gym.error.Error:
    print("Environment 'cocco-v2' is already registered. Skipping.")


def set_new_learning_rate(model: SAC, new_lr: float) -> SAC:
    """
    为 SBX SAC 模型设置一个新的学习率。
    
    这个函数通过创建一个带有新学习率的全新优化器和训练状态 (TrainState)
    来实现这一点，同时保留所有已经训练好的网络权重。

    :param model: 要修改的 SBX SAC 模型实例。
    :param new_lr: 新的学习率。
    :return: 经过修改的模型实例。
    """
    print(f"\n[Updater] 正在为模型设置新的学习率...")
    print(f"  - 原始学习率: {model.learning_rate}")
    print(f"  - 新的学习率: {new_lr}")

    # 步骤 1: 创建一个带有新学习率的全新优化器
    new_optimizer = optax.adam(learning_rate=new_lr)

    # 步骤 2: 从旧的 state 中提取出需要保留的部分
    old_state = model.state
    
    # 步骤 3: 使用旧的参数和新的优化器创建一个全新的 SACState
    new_state = SACState.create(
        apply_fn=old_state.apply_fn,
        params=old_state.params,
        target_params=old_state.target_params,
        tx=new_optimizer,  # 注入新的优化器
    )
    
    # 步骤 4: 用新创建的 state 替换模型中旧的 state
    model.state = new_state
    
    # 步骤 5: 更新模型自身的 learning_rate 属性以保持一致性
    model.learning_rate = new_lr
    
    print("[Updater] 学习率更新完成。")
    return model


def run_resume(args):
    """
    主运行函数，用于加载现有模型并使用更小的学习率继续训练。
    """
    # --- 1. 路径和文件校验 ---
    if not os.path.isdir(args.run_dir):
        raise ValueError(f"指定的运行目录不存在: {args.run_dir}")

    model_dir = args.run_dir
    best_model_dir = os.path.join(model_dir, "best_model")
    save_best_model_dir = os.path.join(model_dir, "resumed_best_model","best_model")
    model_path = os.path.join(best_model_dir, "best_model.zip")
    stats_path = os.path.join(best_model_dir, "vecnormalize_best.pkl")

    if not os.path.exists(model_path):
        raise FileNotFoundError(f"模型文件未找到: {model_path}")
    if not os.path.exists(stats_path):
        raise FileNotFoundError(f"VecNormalize 统计文件未找到: {stats_path}")

    print(f"从以下目录恢复训练: {model_dir}")
    print(f"  - 加载模型自: {model_path}")
    print(f"  - 加载统计数据自: {stats_path}")

    # --- 2. 创建环境 ---
    env_kwargs_p1 = {'disable_random_death': True}
    
    vec_env_raw = make_vec_env("cocco-v2", n_envs=24, seed=42, env_kwargs=env_kwargs_p1)
    eval_env_raw = make_vec_env("cocco-v2", n_envs=100, seed=123, env_kwargs=env_kwargs_p1)

    # --- 3. 加载 VecNormalize 统计数据并封装环境 ---
    print("正在加载 VecNormalize 统计数据...")
    vec_env = VecNormalize.load(stats_path, vec_env_raw)
    vec_env.training = True

    eval_env = VecNormalize.load(stats_path, eval_env_raw)
    eval_env.training = False

    # --- 4. 加载 SAC 模型 ---
    print("正在加载 SAC 模型...")
    model = SAC.load(model_path, env=vec_env)
    
    # --- 5. 调用函数设置新的学习率 ---
    model = set_new_learning_rate(model, new_lr=args.new_learning_rate)
    
    # --- 6. 设置 Tensorboard 日志 ---
    run_name_base = os.path.basename(model_dir)
    tensorboard_log_dir = "./tf_logs/cocco_sac_resumed"
    model.tensorboard_log = tensorboard_log_dir
    print(f"Tensorboard 日志将保存到: {tensorboard_log_dir}/{run_name_base}_1")

    # --- 7. 设置评估回调函数 ---
    eval_callback = EvalCallback(
        eval_env,
        best_model_save_path=best_model_dir,
        log_path=best_model_dir,
        eval_freq=1000, 
        n_eval_episodes=300,
        deterministic=True,
    )

    # --- 8. 继续训练 ---
    print(f"\n开始额外 {args.additional_timesteps} 步的训练...")
    model.learn(
        total_timesteps=args.additional_timesteps,
        callback=eval_callback,
        progress_bar=True,
        reset_num_timesteps=False,
        tb_log_name=run_name_base
    )

    # --- 9. 保存续写训练后的最终产出 ---
    final_model_path = os.path.join(model_dir, "final_model_resumed.zip")
    final_stats_path = os.path.join(model_dir, "vecnormalize_final_resumed.pkl")
    model.save(final_model_path)
    vec_env.save(final_stats_path)
    
    print("\n训练已恢复并完成。")
    print(f"续写训练后的最终模型已保存到: {final_model_path}")
    print(f"续写训练后的最终 VecNormalize 统计数据已保存到: {final_stats_path}")


if __name__ == "__main__":
    # 使用一个简单的类来存储参数，方便直接修改
    class Args:
        # !!! 用户修改区域 !!!
        # 请将此路径修改为您要恢复的训练运行所在的目录
        run_dir = "models/cocco_sac_curriculum/cocco_sac_run_4"
        
        # 设置要额外训练的时间步数
        additional_timesteps = 5_000_000

    args = Args()
    run_resume(args)