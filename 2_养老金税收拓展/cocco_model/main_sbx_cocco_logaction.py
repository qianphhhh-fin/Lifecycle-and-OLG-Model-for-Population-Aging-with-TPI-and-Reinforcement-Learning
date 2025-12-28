# 文件名: main_sbx_cocco_logaction.py

import gymnasium as gym
from gymnasium.envs.registration import register
import os
import argparse
from sbx import SAC
from utils.env_util import make_vec_env
from utils.AllocationEvalCallback import EvalCallback
from stable_baselines3.common.vec_env import VecNormalize

# 注册新的 V3 环境
register(
    id='cocco-v3',
    entry_point='v3_cocco_env:CoccoEnvV3',
)

def run(args):
    """
    主运行函数，使用 V3 环境进行课程学习。
    V3 环境内置了对数动作空间、消费截断和仿射奖励变换。
    """
    # --- 路径和模型名称设置 ---
    log_dir = './models/cocco_sac_v3_logaction/'
    base_model_name = 'cocco_sac_run'

    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    
    if args.run_num is None:
        existing_runs = [d for d in os.listdir(log_dir) if d.startswith(base_model_name)]
        run_num = len(existing_runs)
    else:
        run_num = args.run_num
        
    model_dir = os.path.join(log_dir, f"{base_model_name}_{run_num}")
    if not os.path.exists(model_dir):
        os.makedirs(model_dir)
    print(f"模型和日志将保存到: {model_dir}")
    
    phase1_model_path = os.path.join(model_dir, "phase1_model.zip")
    vec_normalize_stats_path = os.path.join(model_dir, "vecnormalize_phase1.pkl")

    # --- 课程学习流程 ---

    # ==========================================================================
    # 阶段一: 预训练 (无随机死亡)
    # ==========================================================================
    if not args.skip_phase1:
        print("\n===== 开始阶段一: 预训练 (无随机死亡) - 环境 V3 =====")
        
        env_kwargs_p1 = {'disable_random_death': True}
        
        # 注意: VecNormalize 仍然是推荐的，它可以将经过我们仿射变换的奖励
        # (例如[-2, 1]) 进一步标准化到均值为0，标准差为1，这对于Actor-Critic算法的稳定性有额外好处。
        vec_env_p1_raw = make_vec_env("cocco-v3", n_envs=12, seed=42, env_kwargs=env_kwargs_p1)
        eval_env_p1_raw = make_vec_env("cocco-v3", n_envs=1, seed=123, env_kwargs=env_kwargs_p1)

        vec_env_p1 = VecNormalize(vec_env_p1_raw, norm_obs=False, norm_reward=True)
        eval_env_p1 = VecNormalize(eval_env_p1_raw, norm_obs=False, norm_reward=True, training=False)

        model_beta = vec_env_p1.get_attr("beta")[0]

        eval_callback_p1 = EvalCallback(
            eval_env_p1,
            best_model_save_path=os.path.join(model_dir, "best_model_phase1"),
            log_path=os.path.join(model_dir, "logs_phase1"),
            eval_freq=10000,
            n_eval_episodes=300,
            deterministic=True,
            discount_factor=model_beta
        )

        model = SAC(
            'MlpPolicy', vec_env_p1, verbose=1, gamma=model_beta,
            learning_rate=1e-5, buffer_size=1_000_000, batch_size=512,
            ent_coef='auto', train_freq=(1, 'step'), gradient_steps=1,
            learning_starts=1000, policy_kwargs=dict(net_arch=[128,128]),
            tensorboard_log=f"./tf_logs/cocco_sac_v3/run_{run_num}_phase1"
        )
        
        model.learn(total_timesteps=args.timesteps_phase1, callback=eval_callback_p1, progress_bar=True)
        
        model.save(phase1_model_path)
        vec_env_p1.save(vec_normalize_stats_path)
        print(f"阶段一模型已保存到: {phase1_model_path}")

    else:
        print("\n===== 跳过阶段一，直接进入阶段二 (环境 V3) =====")
        if not os.path.exists(phase1_model_path) or not os.path.exists(vec_normalize_stats_path):
            print(f"错误: 找不到阶段一模型或统计数据。请先完整运行阶段一。")
            return

    # ==========================================================================
    # 阶段二: 微调 (有随机死亡)
    # ==========================================================================
    print("\n===== 开始阶段二: 微调 (有随机死亡) - 环境 V3 =====")

    env_kwargs_p2 = {'disable_random_death': False}
    
    vec_env_p2_raw = make_vec_env("cocco-v3", n_envs=64, seed=42, env_kwargs=env_kwargs_p2)
    eval_env_p2_raw = make_vec_env("cocco-v3", n_envs=1, seed=123, env_kwargs=env_kwargs_p2)

    vec_env_p2 = VecNormalize.load(vec_normalize_stats_path, vec_env_p2_raw)
    eval_env_p2 = VecNormalize.load(vec_normalize_stats_path, eval_env_p2_raw)
    eval_env_p2.training = False

    model_beta = vec_env_p2.get_attr("beta")[0]
    
    eval_callback_p2 = EvalCallback(
        eval_env_p2,
        best_model_save_path=os.path.join(model_dir, "best_model_final"),
        log_path=os.path.join(model_dir, "logs_final"),
        eval_freq=15000, n_eval_episodes=500, deterministic=True,
        discount_factor=model_beta
    )

    model = SAC.load(phase1_model_path, env=vec_env_p2)
    
    print("已加载阶段一模型和统计数据，开始微调...")
    
    model.learn(
        total_timesteps=args.timesteps_phase2,
        callback=eval_callback_p2, progress_bar=True, reset_num_timesteps=False
    )
    
    final_model_path = os.path.join(model_dir, "final_model.zip")
    final_stats_path = os.path.join(model_dir, "vecnormalize_final.pkl")
    model.save(final_model_path)
    vec_env_p2.save(final_stats_path)
    print(f"最终模型已保存到 {final_model_path}")


def parse_args():
    parser = argparse.ArgumentParser("Cocco Model V3 (Log-Action) Training Script")
    parser.add_argument("--skip-phase1", action="store_true", help="跳过阶段一预训练。")
    parser.add_argument("--timesteps-phase1", type=int, default=1_000_000, help="阶段一的总时间步数。")
    parser.add_argument("--timesteps-phase2", type=int, default=500_000, help="阶段二的总时间步数。")
    parser.add_argument("--run-num", type=int, default=None, help="手动指定运行编号。")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    run(args)