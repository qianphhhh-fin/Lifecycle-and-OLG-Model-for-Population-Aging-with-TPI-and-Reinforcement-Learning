from collections.abc import Sequence
from typing import Any, Callable, Optional, Union

import flax.linen as nn
import jax
import jax.numpy as jnp
import numpy as np
import optax
import tensorflow_probability.substrates.jax as tfp
from gymnasium import spaces
from stable_baselines3.common.type_aliases import Schedule
from flax.training.train_state import TrainState

# --- 新增导入 ---
from sbx.common.policies import VectorCritic
# ----------------

from sbx.sac.policies import SACPolicy
from sbx.common.distributions import TanhTransformedDistribution
from sbx.common.policies import Flatten
from sbx.common.type_aliases import RLTrainState


tfd = tfp.distributions


class MultiHeadSquashedGaussianActor(nn.Module):
    """
    一个具有分离头部的自定义Actor网络，用于处理具有不同梯度信号的动作。

    - 共享主干网络处理输入状态。
    - 两个独立的头部网络分别输出 c_prop 和 alpha 的均值与对数标准差。
    """
    net_arch: Sequence[int]
    head_net_arch: Sequence[int]  # 为每个头部定义的网络结构
    action_dim: int
    log_std_min: float = -20
    log_std_max: float = 2
    activation_fn: Callable[[jnp.ndarray], jnp.ndarray] = nn.relu

    def get_std(self):
        return jnp.array(0.0)

    @nn.compact
    def __call__(self, x: jnp.ndarray) -> tfd.Distribution:
        # 1. 共享主干网络
        shared_features = Flatten()(x)
        for n_units in self.net_arch:
            shared_features = nn.Dense(n_units)(shared_features)
            shared_features = self.activation_fn(shared_features)

        # 2. 分离的头部
        # 确保 action_dim 为 2
        if self.action_dim != 2:
            raise ValueError(f"MultiHeadSquashedGaussianActor is designed for action_dim=2, but got {self.action_dim}")

        # --- 头部 1: c_prop (动作维度 0) ---
        c_prop_head = shared_features
        for n_units in self.head_net_arch:
            c_prop_head = nn.Dense(n_units)(c_prop_head)
            c_prop_head = self.activation_fn(c_prop_head)
        mean_c_prop = nn.Dense(1)(c_prop_head)
        log_std_c_prop = nn.Dense(1)(c_prop_head)

        # --- 头部 2: alpha (动作维度 1) ---
        alpha_head = shared_features
        for n_units in self.head_net_arch:
            alpha_head = nn.Dense(n_units)(alpha_head)
            alpha_head = self.activation_fn(alpha_head)
        mean_alpha = nn.Dense(1)(alpha_head)
        log_std_alpha = nn.Dense(1)(alpha_head)

        # 3. 合并头部的输出
        mean = jnp.concatenate([mean_c_prop, mean_alpha], axis=-1)
        log_std = jnp.concatenate([log_std_c_prop, log_std_alpha], axis=-1)
        
        # 4. 创建分布 (与原始 sbx Actor 相同)
        log_std = jnp.clip(log_std, self.log_std_min, self.log_std_max)
        dist = TanhTransformedDistribution(
            tfd.MultivariateNormalDiag(loc=mean, scale_diag=jnp.exp(log_std)),
        )
        return dist


class CustomSACPolicy(SACPolicy):
    """
    一个使用我们自定义多头Actor的SAC策略类。
    """
    def __init__(
        self,
        observation_space: spaces.Space,
        action_space: spaces.Box,
        lr_schedule: Schedule,
        head_net_arch: Optional[Sequence[int]] = None,
        **kwargs,
    ):
        if head_net_arch is None:
            head_net_arch = [64]
        self.head_net_arch = head_net_arch

        if "actor_class" in kwargs:
            del kwargs["actor_class"]
        # 也移除 vector_critic_class 以防万一
        if "vector_critic_class" in kwargs:
            del kwargs["vector_critic_class"]
        
        super().__init__(
            observation_space, 
            action_space, 
            lr_schedule, 
            **kwargs
        )
        
        # --- 核心修改 ---
        # 手动覆盖 actor_class 和 vector_critic_class 属性
        self.actor_class = MultiHeadSquashedGaussianActor
        self.vector_critic_class = VectorCritic
        # -----------------

    def build(self, key: jax.Array, lr_schedule: Schedule, qf_learning_rate: float) -> jax.Array:
        # 在创建Actor时，需要将 head_net_arch 参数传递进去
        # 覆盖原始的 build 方法以注入新参数
        
        key, actor_key, qf_key, dropout_key = jax.random.split(key, 4)
        key, self.key = jax.random.split(key, 2)
        self.reset_noise()

        if isinstance(self.observation_space, spaces.Dict):
            obs = jnp.array([spaces.flatten(self.observation_space, self.observation_space.sample())])
        else:
            obs = jnp.array([self.observation_space.sample()])
        action = jnp.array([self.action_space.sample()])

        self.actor = self.actor_class(
            action_dim=int(np.prod(self.action_space.shape)),
            net_arch=self.net_arch_pi,
            head_net_arch=self.head_net_arch,  # <--- 注入新参数
            activation_fn=self.activation_fn,
        )

        self.actor.reset_noise = self.reset_noise
        
        optimizer_class = optax.inject_hyperparams(self.optimizer_class)(learning_rate=lr_schedule(1), **self.optimizer_kwargs)

        self.actor_state = TrainState.create(
            apply_fn=self.actor.apply,
            params=self.actor.init(actor_key, obs),
            tx=optimizer_class,
        )

        self.qf = self.vector_critic_class(
            dropout_rate=self.dropout_rate,
            use_layer_norm=self.layer_norm,
            net_arch=self.net_arch_qf,
            n_critics=self.n_critics,
            activation_fn=self.activation_fn,
        )
        
        optimizer_class_qf = optax.inject_hyperparams(self.optimizer_class)(
            learning_rate=qf_learning_rate, **self.optimizer_kwargs
        )

        self.qf_state = RLTrainState.create(
            apply_fn=self.qf.apply,
            params=self.qf.init(
                {"params": qf_key, "dropout": dropout_key},
                obs,
                action,
            ),
            target_params=self.qf.init(
                {"params": qf_key, "dropout": dropout_key},
                obs,
                action,
            ),
            tx=optimizer_class_qf,
        )

        self.actor.apply = jax.jit(self.actor.apply)
        self.qf.apply = jax.jit(
            self.qf.apply,
            static_argnames=("dropout_rate", "use_layer_norm"),
        )

        return key