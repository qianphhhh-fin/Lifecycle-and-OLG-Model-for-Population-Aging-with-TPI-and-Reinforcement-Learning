import sys
import os
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from main_olg_v9_utils_simplified import OLG_V9_Utils
import numpy as np

# 初始化参数
cS = OLG_V9_Utils.ParameterValues_HuggettStyle()
leLogGridV, leTrProbM, leProb1V = OLG_V9_Utils.EarningProcess_olgm(cS)

print('效率状态数:', cS.nw)
print('初始概率分布 leProb1V:', leProb1V)
print('转移概率矩阵 leTrProbM:')
print(leTrProbM)

# 测试一些随机转移
np.random.seed(42)
current_eps_idx = 1  # 从状态3开始（1-based）
print(f'\n从状态 {current_eps_idx+1} 开始转移测试:')

for i in range(10):
    trans_probs = leTrProbM[current_eps_idx, :]
    new_eps_idx = np.random.choice(len(trans_probs), p=trans_probs)
    print(f'步骤 {i+1}: {current_eps_idx+1} -> {new_eps_idx+1}')
    current_eps_idx = new_eps_idx 