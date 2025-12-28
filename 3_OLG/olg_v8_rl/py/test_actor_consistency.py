import os
# 解决OpenMP库冲突问题
os.environ['KMP_DUPLICATE_LIB_OK'] = 'TRUE'

import numpy as np
import json
from pathlib import Path

# 切换到 SB3
try:
    from stable_baselines3 import SAC
    SB3_AVAILABLE = True
except ImportError:
    SB3_AVAILABLE = False
    print("❌ 错误: Stable Baselines 3 (SB3) 未安装。请运行: pip install stable-baselines3")
    exit()

# 导入MATLAB Engine
try:
    import matlab.engine
    MATLAB_AVAILABLE = True
except ImportError:
    MATLAB_AVAILABLE = False
    print("❌ 错误: MATLAB Engine 未安装。请运行: pip install matlabengine")
    exit()

def test_actor_consistency():
    """
    测试并比较原始Python SAC模型与MATLAB中重建的模型的输出是否一致。
    """
    print("=" * 80)
    print("🔬 SAC Actor 输出一致性测试 (Python vs MATLAB)")
    print("=" * 80)

    # --- 1. 定义模型和参数路径 (更新为 SB3 路径) ---
    best_model_path = Path('./py/best_model_sb3/best_model.zip')
    params_path = Path('./py/best_model_sb3/best_model_params.json')
    matlab_function_dir = Path.cwd() # 假设 .m 文件在根目录

    # 检查文件是否存在
    if not best_model_path.exists():
        print(f"❌ 错误: 找不到SAC模型文件: {best_model_path}")
        return
    if not params_path.exists():
        print(f"❌ 错误: 找不到参数文件: {params_path}")
        return

    # --- 2. 构建一个测试用的观测(Observation)向量 ---
    # 这个向量的结构必须与 reconstruct_rl.m 和 OLG 环境中的完全一致
    # [k, kpps, age_idx, eps_idx, r-1, w, TR, b_payg, pps_tax, is_retired]
    test_obs = np.array([
        0.5,      # k:      资产
        0.2,      # kpps:   个人养老金资产
        25,       # age_idx: 年龄 (索引从1开始)
        4,        # eps_idx: 生产力冲击 (索引从1开始)
        1.03 - 1, # r-1:    无风险利率
        2.0,      # w:      工资率
        0.1,      # TR:     转移支付
        0.4,      # b_payg: PAYG养老金
        0.1,      # pps_tax: PPS税率
        0         # is_retired: 是否退休 (0=否, 1=是)
    ], dtype=np.float32)
    print("🧪 使用测试观测向量:")
    print(test_obs)
    print("-" * 40)

    # --- 3. 从Python SAC模型获取预测 ---
    print("🐍 步骤 1: 从原始Python (SB3)模型获取预测...")
    try:
        model = SAC.load(best_model_path)
        # 使用确定性预测 (deterministic=True)
        py_action, _ = model.predict(test_obs, deterministic=True)
        print(f"✅ Python (SB3) 模型预测成功!")
        print(f"   -> 输出动作: {py_action}")
    except Exception as e:
        print(f"❌ 错误: 加载或预测Python模型时出错: {e}")
        return
    print("-" * 40)

    # --- 4. 从MATLAB重建的网络获取预测 ---
    print("🇲 步骤 2: 从MATLAB重建的网络获取预测...")
    mat_action = None
    try:
        print("   - 正在启动 MATLAB Engine...")
        eng = matlab.engine.start_matlab()
        print("   - MATLAB Engine 启动成功!")
        # 将工作目录添加到MATLAB路径
        eng.addpath(str(matlab_function_dir), nargout=0)
        print(f"   - 已将 '{matlab_function_dir}' 添加到MATLAB路径")
        
        # 将numpy数组转换为MATLAB可以接收的格式
        matlab_obs = matlab.double(test_obs.tolist())
        
        # 调用我们创建的MATLAB函数
        # nargout=1 确保函数返回一个值
        print("   - 正在调用 'predict_from_reconstructed_net.m'...")
        mat_action_raw = eng.predict_from_reconstructed_net(matlab_obs, nargout=1)
        
        # 将MATLAB返回的结果转换为numpy数组
        mat_action = np.array(mat_action_raw[0])
        print(f"✅ MATLAB 模型预测成功!")
        print(f"   -> 输出动作: {mat_action}")
        
        # 关闭MATLAB引擎
        eng.quit()
        print("   - MATLAB Engine 已关闭。")

    except Exception as e:
        print(f"❌ 错误: 调用MATLAB时出错: {e}")
        if 'eng' in locals():
            eng.quit()
        return
    print("-" * 40)
    
    # --- 5. 比较结果 ---
    print("📊 步骤 3: 比较结果...")
    if py_action is not None and mat_action is not None:
        difference = np.abs(py_action - mat_action)
        total_diff = np.sum(difference)
        
        print(f"   - Python (SB3)    输出: {py_action}")
        print(f"   - MATLAB (重建) 输出: {mat_action}")
        print(f"   - 绝对差值:            {difference}")
        print(f"   - 总绝对差值:          {total_diff:.10f}")
        
        # 设定一个非常小的容忍度来判断是否一致
        tolerance = 1e-6
        if total_diff < tolerance:
            print("\n✅ 结论: 两个模型的输出高度一致！模型重建成功。")
        else:
            print("\n⚠️ 结论: 两个模型的输出存在显著差异。请检查模型导出/重建过程。")
    else:
        print("❌ 未能完成比较，因为其中一个步骤失败了。")
        
    print("=" * 80)

if __name__ == "__main__":
    test_actor_consistency() 