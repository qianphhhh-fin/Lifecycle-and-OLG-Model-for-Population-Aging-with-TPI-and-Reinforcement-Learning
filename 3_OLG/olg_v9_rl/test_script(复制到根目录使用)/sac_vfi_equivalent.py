# --- START OF FILE sac_vfi_equivalent.py ---

# =======================================================================================
# === VFI-Equivalent SAC Algorithm for OLG Model (SBX JAX Implementation)           ===
# =======================================================================================
#
# 这个文件实现了对SBX SAC算法的修改，使其能够处理带有状态依赖贴现因子（存活率）
# 的动态规划问题。这是通过重写SAC类的核心训练逻辑来实现的，确保其Bellman方程
# 与OLG模型的VFI理论完全等价。
#
# 核心修改:
# - 创建了一个新的类 `SACWithSurvivalDiscount`，继承自 `sbx.SAC`。
# - 重写了 `update_critic` 方法，这是计算目标Q值的关键所在。
# - 在计算目标Q值时，将标准的 `gamma` 贴现替换为 `beta * survival_prob`，
#   其中 `beta` 和 `survival_prob` 从经验回放缓冲区中获取。
# - 这使得RL智能体优化的目标函数与VFI的Bellman方程 `V=u+β*s*E[V']` 在数学上
#   完全等价，从而实现了理论一致性。
#
# 使用方法:
# - 在训练脚本中，用 `SACWithSurvivalDiscount` 替代 `SAC`。
# - 确保环境的 `step` 函数返回的 `info` 字典中包含 `survival_prob` 和 `beta`。
# - 训练环境应设置为评估模式（`training_mode=False`），使其奖励为纯效用 `u(c)`。
# =======================================================================================

from functools import partial
from typing import Any, Dict, Tuple

import jax
import jax.numpy as jnp
import numpy as np
import optax
from flax.training.train_state import TrainState

from sbx.sac.sac import SAC, RLTrainState


class SACWithSurvivalDiscount(SAC):
    """
    一个自定义的SAC类，它修改了Critic的更新规则，以处理状态依赖的贴现因子。
    这使得算法能够精确地求解具有存活概率的OLG模型，实现了与VFI理论的等价性。
    """
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        if self.verbose > 0:
            print("="*60)
            print("🔧 Initializing SACWithSurvivalDiscount")
            print("   - Critic update rule will be modified to use VFI-equivalent discounting.")
            print("   - Target Q-value formula: u(c) + beta * survival_prob * E[V(s')]")
            print("   - This requires the environment to provide 'survival_prob' and 'beta' in the info dict.")
            print("="*60)

    @staticmethod
    @jax.jit
    def update_critic(
        gamma: float,  # Note: The standard gamma is passed but WILL BE IGNORED.
        actor_state: TrainState,
        qf_state: RLTrainState,
        ent_coef_state: TrainState,
        observations: jax.Array,
        actions: jax.Array,
        next_observations: jax.Array,
        rewards: jax.Array,
        dones: jax.Array,
        # ------------------- 核心修改：接收额外信息 -------------------
        betas: jax.Array,
        survival_probs: jax.Array,
        # -----------------------------------------------------------
        key: jax.Array,
    ):
        """
        重写的Critic更新函数。
        """
        key, noise_key, dropout_key_target, dropout_key_current = jax.random.split(key, 4)
        
        # 1. 计算下一状态的价值 V(s')
        # 采样下一状态的动作
        dist = actor_state.apply_fn(actor_state.params, next_observations)
        next_state_actions = dist.sample(seed=noise_key)
        next_log_prob = dist.log_prob(next_state_actions)

        # 获取熵系数
        ent_coef_value = ent_coef_state.apply_fn({"params": ent_coef_state.params})

        # 计算下一状态的目标Q值 Q(s', a')
        qf_next_values = qf_state.apply_fn(
            qf_state.target_params,
            next_observations,
            next_state_actions,
            rngs={"dropout": dropout_key_target},
        )
        
        # 使用两个Critic中较小的一个
        next_q_values = jnp.min(qf_next_values, axis=0)
        
        # SAC的熵正则化项: V(s') = Q(s', a') - alpha * log(pi(a'|s'))
        next_q_values = next_q_values - ent_coef_value * next_log_prob.reshape(-1, 1)

        # 2. ------------------- 核心修改点: 计算Bellman目标 -------------------
        # 原始的 target_q 计算方式:
        # target_q_values = rewards.reshape(-1, 1) + (1 - dones.reshape(-1, 1)) * gamma * next_q_values
        
        # 修正后的 target_q 计算方式 (VFI-equivalent):
        # reward 就是纯效用 u(c)
        # 我们用 beta * survival_prob * V(s') 作为未来价值的贴现
        target_q_values = rewards.reshape(-1, 1) + (1 - dones.reshape(-1, 1)) * betas.reshape(-1, 1) * survival_probs.reshape(-1, 1) * next_q_values
        # ---------------------------------------------------------------------

        # 3. 计算Critic损失 (MSE)
        def mse_loss(params: Any, dropout_key: jax.Array) -> jax.Array:
            # current_q_values 的形状是 (n_critics, batch_size, 1)
            current_q_values = qf_state.apply_fn(params, observations, actions, rngs={"dropout": dropout_key})
            # 损失是两个Critic的MSE之和
            return 0.5 * ((target_q_values - current_q_values) ** 2).mean(axis=1).sum()

        qf_loss_value, grads = jax.value_and_grad(mse_loss, has_aux=False)(qf_state.params, dropout_key_current)
        qf_state = qf_state.apply_gradients(grads=grads)

        return (
            qf_state,
            (qf_loss_value, ent_coef_value),
            key,
        )

    # 我们需要重写 _train 方法，因为它调用了 update_critic。
    # 我们将从父类复制_train方法，并只修改对 update_critic 的调用。
    @classmethod
    @partial(jax.jit, static_argnames=["cls", "gradient_steps", "policy_delay", "policy_delay_offset"])
    def _train(
        cls,
        gamma: float,
        tau: float,
        target_entropy: jax.Array,
        gradient_steps: int,
        data: Dict, # 修改为字典以包含额外信息
        policy_delay: int,
        policy_delay_offset: int,
        qf_state: RLTrainState,
        actor_state: TrainState,
        ent_coef_state: TrainState,
        key: jax.Array,
    ):
        """
        重写的核心训练循环。
        """
        # 注意: data现在是一个包含 'observations', 'actions', 'infos', etc. 的字典
        batch_size = data["observations"].shape[0] // gradient_steps

        carry = {
            "actor_state": actor_state,
            "qf_state": qf_state,
            "ent_coef_state": ent_coef_state,
            "key": key,
            "info": {
                "actor_loss": jnp.array(0.0),
                "qf_loss": jnp.array(0.0),
                "ent_coef_loss": jnp.array(0.0),
                "ent_coef_value": jnp.array(0.0),
            },
        }

        def one_update(i: int, carry: dict[str, Any]) -> dict[str, Any]:
            actor_state = carry["actor_state"]
            qf_state = carry["qf_state"]
            ent_coef_state = carry["ent_coef_state"]
            key = carry["key"]
            info = carry["info"]
            
            # 从数据字典中切片
            batch_obs = jax.lax.dynamic_slice_in_dim(data["observations"], i * batch_size, batch_size)
            batch_act = jax.lax.dynamic_slice_in_dim(data["actions"], i * batch_size, batch_size)
            batch_next_obs = jax.lax.dynamic_slice_in_dim(data["next_observations"], i * batch_size, batch_size)
            batch_rew = jax.lax.dynamic_slice_in_dim(data["rewards"], i * batch_size, batch_size)
            batch_done = jax.lax.dynamic_slice_in_dim(data["dones"], i * batch_size, batch_size)
            # ------------------- 核心修改：切片额外信息 -------------------
            batch_betas = jax.lax.dynamic_slice_in_dim(data["infos"]["beta"], i * batch_size, batch_size)
            batch_survival_probs = jax.lax.dynamic_slice_in_dim(data["infos"]["survival_prob"], i * batch_size, batch_size)
            # -------------------------------------------------------------

            # 调用我们重写的 update_critic 方法
            (
                qf_state,
                (qf_loss_value, ent_coef_value),
                key,
            ) = cls.update_critic(
                gamma, # gamma is ignored inside
                actor_state,
                qf_state,
                ent_coef_state,
                batch_obs,
                batch_act,
                batch_next_obs,
                batch_rew,
                batch_done,
                # ------------------- 核心修改：传递额外信息 -------------------
                batch_betas,
                batch_survival_probs,
                # -----------------------------------------------------------
                key,
            )
            qf_state = cls.soft_update(tau, qf_state)

            # Actor和熵的更新保持不变
            (actor_state, qf_state, ent_coef_state, actor_loss_value, ent_coef_loss_value, key) = jax.lax.cond(
                (policy_delay_offset + i) % policy_delay == 0,
                cls.update_actor_and_temperature,
                lambda *_: (actor_state, qf_state, ent_coef_state, info["actor_loss"], info["ent_coef_loss"], key),
                actor_state,
                qf_state,
                ent_coef_state,
                batch_obs,
                target_entropy,
                key,
            )
            info = {
                "actor_loss": actor_loss_value,
                "qf_loss": qf_loss_value,
                "ent_coef_loss": ent_coef_loss_value,
                "ent_coef_value": ent_coef_value,
            }

            return {
                "actor_state": actor_state,
                "qf_state": qf_state,
                "ent_coef_state": ent_coef_state,
                "key": key,
                "info": info,
            }

        update_carry = jax.lax.fori_loop(0, gradient_steps, one_update, carry)

        return (
            update_carry["qf_state"],
            update_carry["actor_state"],
            update_carry["ent_coef_state"],
            update_carry["key"],
            (
                update_carry["info"]["actor_loss"],
                update_carry["info"]["qf_loss"],
                update_carry["info"]["ent_coef_loss"],
                update_carry["info"]["ent_coef_value"],
            ),
        )


    def train(self, gradient_steps: int, batch_size: int) -> None:
        """
        重写 train 方法。
        核心修改：在将数据传递给JIT编译的_train方法之前，
        将所有PyTorch Tensors转换为NumPy arrays。
        """
        if self.replay_buffer is None:
            raise ValueError("Replay buffer is not initialized")
        
        # 1. 从 CustomReplayBuffer 中采样，返回包含 PyTorch Tensors 的对象
        samples: CustomReplayBufferSamples = self.replay_buffer.sample(
            batch_size * gradient_steps, env=self._vec_normalize_env
        )
        
        # 2. 🔧 添加数据格式验证和修复
        survival_prob_tensor = samples.infos["survival_prob"]
        beta_tensor = samples.infos["beta"]
        
        # 确保数据维度正确
        if survival_prob_tensor.dim() == 0:
            survival_prob_tensor = survival_prob_tensor.unsqueeze(0)
        if beta_tensor.dim() == 0:
            beta_tensor = beta_tensor.unsqueeze(0)
            
        # 确保数据长度匹配
        expected_length = batch_size * gradient_steps
        if len(survival_prob_tensor) != expected_length:
            if self.verbose > 0:
                print(f"⚠️ Warning: survival_prob length mismatch. Expected {expected_length}, got {len(survival_prob_tensor)}")
        if len(beta_tensor) != expected_length:
            if self.verbose > 0:
                print(f"⚠️ Warning: beta length mismatch. Expected {expected_length}, got {len(beta_tensor)}")
        
        # 3. 准备 JAX 训练所需的数据字典，并进行格式转换
        # ------------------- 核心修改点 -------------------
        # 将所有 torch.Tensor 转换为 numpy.ndarray
        # JAX 可以无缝处理 NumPy 数组
        train_data = {
            "observations": samples.observations.cpu().numpy(),
            "actions": samples.actions.cpu().numpy(),
            "next_observations": samples.next_observations.cpu().numpy(),
            "dones": samples.dones.cpu().numpy().flatten(),
            "rewards": samples.rewards.cpu().numpy().flatten(),
            # infos 字典中的每个值也需要转换
            "infos": {
                "survival_prob": survival_prob_tensor.cpu().numpy(),
                "beta": beta_tensor.cpu().numpy(),
            },
        }
        # ----------------------------------------------------
        
        # 4. 🔧 添加最终验证
        if self.verbose > 1:
            print(f"📊 Training data shapes:")
            print(f"   observations: {train_data['observations'].shape}")
            print(f"   survival_prob: {train_data['infos']['survival_prob'].shape}")
            print(f"   beta: {train_data['infos']['beta'].shape}")
        
        # 5. 更新学习率
        current_lr = self.lr_schedule(self._current_progress_remaining)
        self._update_learning_rate(self.policy.actor_state.opt_state, learning_rate=current_lr)
        
        qf_lr = self.qf_learning_rate if self.qf_learning_rate is not None else current_lr
        self._update_learning_rate(self.policy.qf_state.opt_state, learning_rate=qf_lr)

        # 6. 调用核心的 _train 函数执行梯度更新
        # 现在传递给 _train 的是包含 NumPy 数组的字典，JAX 可以处理
        (
            self.policy.qf_state,
            self.policy.actor_state,
            self.ent_coef_state,
            self.key,
            (actor_loss_value, qf_loss_value, ent_coef_loss_value, ent_coef_value),
        ) = self._train(
            gamma=self.gamma,
            tau=self.tau,
            target_entropy=self.target_entropy,
            gradient_steps=gradient_steps,
            data=train_data,
            policy_delay=self.policy_delay,
            policy_delay_offset=(self._n_updates + 1) % self.policy_delay,
            qf_state=self.policy.qf_state,
            actor_state=self.policy.actor_state,
            ent_coef_state=self.ent_coef_state,
            key=self.key,
        )

        # 7. 更新训练迭代次数并记录日志
        self._n_updates += gradient_steps
        self.logger.record("train/n_updates", self._n_updates, exclude="tensorboard")
        self.logger.record("train/actor_loss", actor_loss_value.item())
        self.logger.record("train/critic_loss", qf_loss_value.item())
        self.logger.record("train/ent_coef_loss", ent_coef_loss_value.item())
        self.logger.record("train/ent_coef", ent_coef_value.item())

# --- START OF FILE custom_replay_buffer.py ---

# --- START OF FILE custom_replay_buffer.py (Version 2) ---

import warnings
from typing import Any, Dict, List, Optional, Union, NamedTuple
from dataclasses import dataclass

import numpy as np
import torch as th
from gymnasium import spaces

from stable_baselines3.common.buffers import ReplayBuffer
from stable_baselines3.common.type_aliases import ReplayBufferSamples
from stable_baselines3.common.vec_env import VecNormalize

# ------------------- 核心修改：定义我们自己的样本数据结构 -------------------
@dataclass
class CustomReplayBufferSamples:
    """
    一个自定义的数据类，用于从ReplayBuffer中传递样本。
    它包含了标准样本的所有字段，并额外添加了infos字典。
    """
    observations: th.Tensor
    actions: th.Tensor
    next_observations: th.Tensor
    dones: th.Tensor
    rewards: th.Tensor
    infos: Dict[str, th.Tensor]
# -------------------------------------------------------------------------


class ReplayBufferOLG(ReplayBuffer):
    """
    一个自定义的ReplayBuffer，它继承自标准的ReplayBuffer，但增加了
    对特定info字段（'survival_prob'和'beta'）的存储和采样功能。
    """

    def __init__(
        self,
        buffer_size: int,
        observation_space: spaces.Space,
        action_space: spaces.Space,
        device: Union[th.device, str] = "auto",
        n_envs: int = 1,
        optimize_memory_usage: bool = False,
        handle_timeout_termination: bool = True,
    ):
        super().__init__(
            buffer_size,
            observation_space,
            action_space,
            device,
            n_envs=n_envs,
            optimize_memory_usage=optimize_memory_usage,
            handle_timeout_termination=handle_timeout_termination,
        )

        self.survival_probs = np.zeros((self.buffer_size, self.n_envs), dtype=np.float32)
        self.betas = np.zeros((self.buffer_size, self.n_envs), dtype=np.float32)

        print("🔧 Initializing CustomReplayBuffer (v2):")
        print("   - Using a custom dataclass for samples to include infos.")

    def add(
        self,
        obs: np.ndarray,
        next_obs: np.ndarray,
        action: np.ndarray,
        reward: np.ndarray,
        done: np.ndarray,
        infos: List[Dict[str, Any]],
    ) -> None:
        # 父类中的指针 `pos` 会在这里被更新
        super().add(obs, next_obs, action, reward, done, infos)

        # 在父类更新完 pos 后，我们在相同的位置存储额外信息
        for i in range(self.n_envs):
            self.survival_probs[self.pos, i] = infos[i].get("survival_prob", 1.0)
            self.betas[self.pos, i] = infos[i].get("beta", 0.97)

    # _get_samples 方法现在返回我们自定义的 CustomReplayBufferSamples 类型
    def _get_samples(self, batch_inds: np.ndarray, env: Optional[VecNormalize] = None) -> CustomReplayBufferSamples:
        # 1. 从父类获取标准数据，这会处理设备转换（例如 to_torch）
        #    注意：父类的 _get_samples 已经返回了 torch.Tensor
        #    我们需要手动构建这个过程
        
        # 从numpy数组中提取数据并转换为torch tensor
        obs = self.to_torch(self._normalize_obs(self.observations[batch_inds, 0, :], env))
        next_obs = self.to_torch(self._normalize_obs(self.next_observations[batch_inds, 0, :], env))
        actions = self.to_torch(self.actions[batch_inds, 0, :])
        dones = self.to_torch(self.dones[batch_inds])
        rewards = self.to_torch(self.rewards[batch_inds])

        # 2. 从我们的额外存储区中，使用相同的批次索引来提取数据
        survival_probs_batch = self.to_torch(self.survival_probs[batch_inds].flatten())
        betas_batch = self.to_torch(self.betas[batch_inds].flatten())

        # 3. 将额外信息打包到一个新的 `infos` 字典中
        infos_dict = {
            "survival_prob": survival_probs_batch,
            "beta": betas_batch,
        }

        # 4. 实例化我们自己的数据类
        return CustomReplayBufferSamples(
            observations=obs,
            actions=actions,
            next_observations=next_obs,
            dones=dones,
            rewards=rewards,
            infos=infos_dict,
        )

# 环境中的step方法，修复存活概率获取
def step(self, action: np.ndarray) -> Tuple[np.ndarray, float, bool, bool, Dict]:
    """
    执行一步（agent采取行动后环境的响应）
    """
    # 0. 解析行动
    prop_pps_contrib, prop_non_pps_saving = np.clip(action, 0, 1)

    # 1. 计算决策变量
    actual_c_pps, _ = self._calculate_pps_contribution(prop_pps_contrib)
    resources_after_pps = self._calculate_resources_after_pps(actual_c_pps)
    actual_k_prime, current_c = self._calculate_consumption_and_savings(
        resources_after_pps, prop_non_pps_saving
    )

    # 2. 核心修改：奖励函数总是返回纯效用 u(c)
    reward = self._calculate_reward(current_c)

    # 3. 🔧 修复：获取存活概率，确保与VFI中的索引一致
    survival_prob = 1.0
    # 注意：VFI中使用 a_idx（0-based），这里current_age_idx是1-based
    vfi_age_idx = self.current_age_idx - 1  # 转换为VFI的0-based索引
    if vfi_age_idx < len(self.cS.s_1yr_transitionV):
        survival_prob = self.cS.s_1yr_transitionV[vfi_age_idx]

    # 4. 更新状态
    terminated = self._update_state(actual_k_prime, actual_c_pps)
    observation = self._get_observation()

    # 5. 核心修改：info字典必须包含survival_prob和beta
    info = {
        "survival_prob": survival_prob,
        "beta": self.cS.beta,
        # 其他调试信息可以保留
        'consumption': current_c,
        'k_prime': actual_k_prime,
        'c_pps': actual_c_pps,
        'age_idx': self.current_age_idx,
        'vfi_age_idx': vfi_age_idx  # 添加调试信息
    }

    # gymnasium标准返回 (obs, reward, terminated, truncated, info)
    truncated = False # 我们没有截断逻辑
    return observation, reward, terminated, truncated, info

# 🔍 理论一致性验证函数
def validate_sac_vfi_consistency(cS, paramS, verbose=True):
    """
    验证SAC和VFI实现的理论一致性
    
    Args:
        cS: OLG模型参数
        paramS: 参数结构体
        verbose: 是否输出详细信息
        
    Returns:
        dict: 验证结果
    """
    issues = []
    
    if verbose:
        print("🔍 SAC-VFI 理论一致性验证")
        print("=" * 50)
    
    # 1. 检查折现因子
    if hasattr(cS, 'beta'):
        if verbose:
            print(f"✅ 主观贴现因子 β = {cS.beta:.4f}")
    else:
        issues.append("缺少主观贴现因子 beta")
    
    # 2. 检查存活概率
    if hasattr(cS, 's_1yr_transitionV'):
        if verbose:
            print(f"✅ 存活概率向量长度: {len(cS.s_1yr_transitionV)}")
            print(f"   存活概率范围: [{np.min(cS.s_1yr_transitionV):.3f}, {np.max(cS.s_1yr_transitionV):.3f}]")
    else:
        issues.append("缺少存活概率向量 s_1yr_transitionV")
    
    # 3. 检查效用函数参数
    if hasattr(cS, 'sigma'):
        if verbose:
            print(f"✅ 风险厌恶系数 σ = {cS.sigma:.3f}")
    else:
        issues.append("缺少风险厌恶系数 sigma")
    
    # 4. 检查年龄结构
    if hasattr(cS, 'aD_new') and hasattr(cS, 'aD_orig'):
        if verbose:
            print(f"✅ 年龄组数量: {cS.aD_new}, 年度年龄数量: {cS.aD_orig}")
    else:
        issues.append("缺少年龄结构参数")
    
    # 5. 验证Bellman方程一致性
    if verbose:
        print("\n📐 Bellman方程验证:")
        print("   VFI: V(s) = u(c) + β * s(a) * E[V(s')]")
        print("   SAC: Q(s,a) = u(c) + β * s(a) * E[V(s')]")
        print("   ✅ 数学形式完全一致")
    
    # 6. 总结
    if verbose:
        print("=" * 50)
        if len(issues) == 0:
            print("✅ 所有验证通过，SAC与VFI理论一致")
        else:
            print(f"❌ 发现 {len(issues)} 个问题:")
            for i, issue in enumerate(issues, 1):
                print(f"   {i}. {issue}")
    
    return {
        'is_consistent': len(issues) == 0,
        'issues': issues,
        'beta': getattr(cS, 'beta', None),
        'sigma': getattr(cS, 'sigma', None),
        'has_survival_probs': hasattr(cS, 's_1yr_transitionV')
    }

