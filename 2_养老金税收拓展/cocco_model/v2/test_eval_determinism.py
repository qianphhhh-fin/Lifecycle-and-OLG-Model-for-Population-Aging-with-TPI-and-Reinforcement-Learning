import time # 导入 time 模块
import gymnasium as gym
import numpy as np
import pandas as pd
from gymnasium.envs.registration import register

from sbx import SAC
from utils.vec_normalize import VecNormalize
from stable_baselines3.common.vec_env import DummyVecEnv
from utils.AllocationEvalCallback import evaluate_policy

# --- [用户配置区] ---
# 请将以下路径和参数修改为你自己的配置
MODEL_PATH = "models/cocco_sac_curriculum/cocco_sac_run_13/best_model/best_model.zip"
STATS_PATH = "models/cocco_sac_curriculum/cocco_sac_run_13/best_model/vecnormalize_best.pkl"
N_EVAL_ENVS = 100         # 使用并行的环境数量 (可以设置为 1 来对比串行性能)
N_EVAL_EPISODES = 300    # 评估的总轮次 (选择一个非 N_EVAL_ENVS 整数倍的数来测试边缘情况)
EVAL_SEED = 123          # 用于评估的固定随机种子
NUM_RUNS = 3             # 重复运行测试的次数
# --- [配置区结束] ---


def run_determinism_test():
    """
    执行 evaluate_policy 函数的确定性测试，并增加计时功能。
    """
    print("=" * 60)
    print("开始执行 evaluate_policy 的确定性与性能测试...")
    print(f"  - 模型路径: {MODEL_PATH}")
    print(f"  - 并行环境数 (n_envs): {N_EVAL_ENVS}")
    print(f"  - 评估 Episodes: {N_EVAL_EPISODES}")
    print(f"  - 随机种子: {EVAL_SEED}")
    print(f"  - 测试运行次数: {NUM_RUNS}")
    print("=" * 60)

    # 1. 注册环境 (使用 try-except 避免重复注册错误)
    try:
        register(id='cocco-v2', entry_point='v2_cocco_env:CoccoEnvV2')
    except gym.error.Error:
        print("环境 'cocco-v2' 已注册，跳过。")

    # 2. 加载模型和环境
    print("\n[步骤 1/3] 加载模型和创建评估环境...")
    model = SAC.load(MODEL_PATH, device='cpu')
    
    # 创建并行的、未封装的评估环境
    eval_env_raw = DummyVecEnv([lambda: gym.make("cocco-v2")] * N_EVAL_ENVS)
    
    # 加载 VecNormalize 统计数据并封装环境
    eval_env = VecNormalize.load(STATS_PATH, eval_env_raw)
    eval_env.training = False 
    eval_env.norm_reward = False 
    print("模型和环境加载成功。")

    # 3. 多次运行 evaluate_policy
    all_results = []
    run_times = [] # --- [MODIFIED] 新增列表用于存储每次运行的时间
    print("\n[步骤 2/3] 多次运行 evaluate_policy 函数...")
    for i in range(NUM_RUNS):
        print(f"  > 开始第 {i + 1}/{NUM_RUNS} 次运行...")
        
        # --- [MODIFIED] 增加计时逻辑 ---
        start_time = time.time()
        results = evaluate_policy(
            model=model,
            env=eval_env,
            n_eval_episodes=N_EVAL_EPISODES,
            deterministic=True,
            seed=EVAL_SEED
        )
        end_time = time.time()
        elapsed_time = end_time - start_time
        run_times.append(elapsed_time)
        # --- [END MODIFICATION] ---
        
        all_results.append(results)
        print(f"  > 第 {i + 1} 次运行完成。耗时: {elapsed_time:.4f} 秒. Mean Reward: {results['mean_reward']:.8f}")

    # 4. 严格比较结果
    print("\n[步骤 3/3] 严格比较所有运行结果...")
    baseline_results = all_results[0]
    is_deterministic = True

    for i in range(1, NUM_RUNS):
        current_results = all_results[i]
        print(f"\n--- 比较 Run 1 和 Run {i + 1} ---")

        # 比较 mean_reward
        try:
            np.testing.assert_equal(baseline_results['mean_reward'], current_results['mean_reward'])
            print(f"  [PASS] mean_reward 匹配: {baseline_results['mean_reward']:.8f}")
        except AssertionError:
            print(f"  [FAIL] mean_reward 不匹配! "
                  f"Run 1: {baseline_results['mean_reward']:.8f}, Run {i+1}: {current_results['mean_reward']:.8f}")
            is_deterministic = False

        # 比较 std_reward
        try:
            np.testing.assert_equal(baseline_results['std_reward'], current_results['std_reward'])
            print(f"  [PASS] std_reward 匹配: {baseline_results['std_reward']:.8f}")
        except AssertionError:
            print(f"  [FAIL] std_reward 不匹配! "
                  f"Run 1: {baseline_results['std_reward']:.8f}, Run {i+1}: {current_results['std_reward']:.8f}")
            is_deterministic = False

        # 比较 mean_profiles DataFrame
        try:
            pd.testing.assert_frame_equal(baseline_results['mean_profiles'], current_results['mean_profiles'])
            print(f"  [PASS] mean_profiles DataFrame 完全匹配。")
        except AssertionError as e:
            print(f"  [FAIL] mean_profiles DataFrame 不匹配!")
            print(e)
            is_deterministic = False

    # --- [MODIFIED] 增加最终的性能总结 ---
    print("=" * 60)
    print("测试摘要:")
    if is_deterministic:
        print("  ✅ 确定性测试: 成功! 所有运行结果完全一致。")
    else:
        print("  ❌ 确定性测试: 失败! 运行结果存在差异。")
    
    mean_time = np.mean(run_times)
    std_time = np.std(run_times)
    print(f"  ⏱️  性能测试: 平均耗时 {mean_time:.4f} 秒 (标准差: {std_time:.4f} 秒)。")
    print("=" * 60)
    # --- [END MODIFICATION] ---


if __name__ == "__main__":
    run_determinism_test()