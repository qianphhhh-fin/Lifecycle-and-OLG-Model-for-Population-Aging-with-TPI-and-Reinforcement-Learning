#
# main_sbx_cocco_resume_train.py
#
import gymnasium as gym
import numpy as np
from gymnasium.envs.registration import register
import os
import jax
from sbx import SAC
from sbx.sac.sac import ConstantEntropyCoef
from flax.training.train_state import TrainState
import optax

from utils.env_util import make_vec_env
from utils.AllocationEvalCallback import EvalCallback
from utils.vec_normalize import VecNormalize
from custom_sac_policy import CustomSACPolicy

# 注册自定义环境
try:
    register(id='cocco-v2', entry_point='v2_cocco_env:CoccoEnvV2')
except gym.error.Error:
    print("Environment 'cocco-v2' is already registered. Skipping.")

class Args:
    """
    用于微调训练的参数配置。
    """
    def __init__(self):
        # --- 输入路径 (需要加载的模型) ---
        model_file = "cocco_sac_run_1"
        self.rl_model = "models/cocco_sac_curriculum/"+  model_file + "/best_model/best_model.zip"
        self.rl_model_rb = "models/cocco_sac_curriculum/"+  model_file + "/best_model/best_model_rb.pkl"
        self.stats = "models/cocco_sac_curriculum/"+  model_file + "/best_model/vecnormalize_best.pkl"
        
        # --- 微调超参数 ---
        self.finetune_timesteps = 2_000_000  # 微调的总步数
        self.finetune_lr = 1e-5             # 微调时使用的新学习率
        # self.finetune_ent_coef = 1e-6       # 微调时使用的接近于零的熵系数

def run_finetuning(args):
    """
    主运行函数，执行模型加载和微调流程。
    """
    print("--- 开始执行 SAC 模型微调流程 ---")
    
    # --- 1. 路径设置 ---
    # 基于原始模型路径，创建一个新的微调文件夹
    original_model_dir = os.path.dirname(os.path.dirname(args.rl_model)) # e.g., .../cocco_sac_run_6
    finetune_model_dir = f"{original_model_dir}_finetuned"
    
    if not os.path.exists(finetune_model_dir):
        os.makedirs(finetune_model_dir)
    print(f"微调后的模型和日志将保存到: {finetune_model_dir}")

    # --- 2. 创建并加载环境 ---
    # 必须使用与原始训练完全相同的环境配置
    env_kwargs = {'disable_random_death': True}
    
    print("创建用于加载的 Dummy VecEnv...")
    # a. 创建原始的 VecEnv 实例
    vec_env_raw = make_vec_env("cocco-v2", n_envs=24, seed=42, env_kwargs=env_kwargs)
    eval_env_raw = make_vec_env("cocco-v2", n_envs=100, seed=123, env_kwargs=env_kwargs)

    # b. 加载归一化统计数据
    print(f"正在从 {args.stats} 加载 VecNormalize 统计数据...")
    vec_env = VecNormalize.load(args.stats, vec_env_raw)
    eval_env = VecNormalize.load(args.stats, eval_env_raw)
    
    # 确保环境处于评估模式（不更新归一化统计数据）
    vec_env.training = False
    eval_env.training = False
    print("环境和统计数据加载完成。")

    # --- 3. 加载预训练模型 ---
    print(f"正在从 {args.rl_model} 加载 SAC 模型...")
    # 加载模型时，我们暂时不设置 tensorboard_log，因为之后会手动设置
    model = SAC.load(
        args.rl_model, 
        env=vec_env, 
        custom_objects={"policy_class": CustomSACPolicy}
    )
    print("模型加载完成。")
    
    # [CRITICAL MODIFICATION] 加载单独保存的 Replay Buffer

    print(f"正在从 {args.rl_model_rb} 加载 Replay Buffer...")
    model.load_replay_buffer(args.rl_model_rb)
    print(f"Replay Buffer 加载完成，当前大小: {model.replay_buffer.pos} / {model.replay_buffer.buffer_size}")


    # [MODIFIED] 手动设置新的 TensorBoard 日志路径
    model.tensorboard_log = model.tensorboard_log + "_resume"
    print(f"微调的 TensorBoard 日志将保存到: {model.tensorboard_log}")


    # --- 4. 修改模型超参数以进行微调 ---
    print("\n--- 正在修改模型超参数以进行微调 ---")
    
    # a. 修改学习率
    print(f"  将 Actor 和 Critic 的学习率设置为: {args.finetune_lr}")
    model.lr_schedule = lambda _: args.finetune_lr  # 更新学习率调度器

    # 手动重建优化器状态以正确注入超参数
    old_actor_opt_state_tuple = model.policy.actor_state.opt_state
    old_qf_opt_state_tuple = model.policy.qf_state.opt_state

    optimizer_class = model.policy.optimizer_class
    optimizer_kwargs = model.policy.optimizer_kwargs or {}
    
    new_actor_optimizer = optax.inject_hyperparams(optimizer_class)(learning_rate=args.finetune_lr, **optimizer_kwargs)
    new_qf_optimizer = optax.inject_hyperparams(optimizer_class)(learning_rate=args.finetune_lr, **optimizer_kwargs)

    new_actor_opt_state_tuple = new_actor_optimizer.init(model.policy.actor_state.params)
    new_qf_opt_state_tuple = new_qf_optimizer.init(model.policy.qf_state.params)
    
    # 从旧状态元组的第一个元素中提取 count, mu, nu，并赋给新状态元组的第一个元素
    new_actor_opt_state_tuple[0].count = old_actor_opt_state_tuple[0].count
    new_actor_opt_state_tuple[0].mu = old_actor_opt_state_tuple[0].mu
    new_actor_opt_state_tuple[0].nu = old_actor_opt_state_tuple[0].nu
    
    new_qf_opt_state_tuple[0].count = old_qf_opt_state_tuple[0].count
    new_qf_opt_state_tuple[0].mu = old_qf_opt_state_tuple[0].mu
    new_qf_opt_state_tuple[0].nu = old_qf_opt_state_tuple[0].nu

    # 用重建后的新状态和优化器替换掉模型中的旧状态
    model.policy.actor_state = model.policy.actor_state.replace(
        tx=new_actor_optimizer,
        opt_state=new_actor_opt_state_tuple
    )
    model.policy.qf_state = model.policy.qf_state.replace(
        tx=new_qf_optimizer,
        opt_state=new_qf_opt_state_tuple
    )
    
    # b. 修改熵系数 (ent_coef)
    # print(f"  将 ent_coef 固定为: {args.finetune_ent_coef}")
    # new_ent_coef_module = ConstantEntropyCoef(ent_coef_init=args.finetune_ent_coef)
    # key, ent_key = jax.random.split(model.key)
    # model.key = key
    # new_ent_coef_state = TrainState.create(
    #     apply_fn=new_ent_coef_module.apply,
    #     params=new_ent_coef_module.init(ent_key)["params"],
    #     tx=optax.adam(learning_rate=0.0),
    # )
    # model.ent_coef = new_ent_coef_module
    # model.ent_coef_state = new_ent_coef_state
    
    print("超参数修改完成。")

        # --- 5. 设置新的评估回调 ---
    print("\n--- 正在设置新的评估回调 ---")
    eval_callback_finetune = EvalCallback(
        eval_env,
        best_model_save_path=os.path.join(finetune_model_dir, "best_model_finetuned"),
        log_path=os.path.join(finetune_model_dir, "best_model_finetuned"),
        eval_freq=1000, 
        n_eval_episodes=300,
        deterministic=True,
    )



    # --- 6. 开始微调训练 ---
    print(f"\n--- 开始微调训练 {args.finetune_timesteps} 步... ---")
    # 如果加载了Buffer，learning_starts可以设为0，因为Buffer已经满了
    # 如果没有加载Buffer，您可能需要保留一个较大的learning_starts
    if os.path.exists(args.rl_model_rb):
        model.learning_starts = 0
    else:
        model.learning_starts = 500_000 
        
    model.learn(
        total_timesteps=args.finetune_timesteps, 
        callback=eval_callback_finetune, 
        progress_bar=True,
        reset_num_timesteps=False, # 关键：不要重置时间步计数器
        tb_log_name="SAC" # 传递一个基础名称给 learn 方法
    )

    # --- 7. 保存最终产出 ---
    final_model_path = os.path.join(finetune_model_dir, "final_model_finetuned.zip")
    model.save(final_model_path)
    print(f"\n微调完成。最终模型已保存到 {final_model_path}")


if __name__ == "__main__":
    args = Args()
    run_finetuning(args)