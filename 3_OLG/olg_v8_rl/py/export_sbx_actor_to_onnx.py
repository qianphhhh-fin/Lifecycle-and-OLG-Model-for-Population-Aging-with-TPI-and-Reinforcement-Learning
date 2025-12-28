"""
export_sbx_actor_to_onnx.py - 将SBX训练的Actor网络导出为ONNX格式

此脚本从SBX模型中提取Actor网络参数，并将其保存为ONNX格式
参考: https://stable-baselines3.readthedocs.io/en/master/guide/export.html
"""

import numpy as np
import os
import pickle
import json
from typing import Dict, Any, List, Tuple
import sys

try:
    from sbx import SAC
    import jax
    import jax.numpy as jnp
    import torch as th
    import onnx
    HAVE_SBX = True
except ImportError:
    print("警告: 未找到SBX或必要的包，将使用有限功能模式")
    HAVE_SBX = False

class OnnxableSACPolicy(th.nn.Module):
    """用于ONNX导出的包装类"""
    def __init__(self, model):
        super().__init__()
        self.model = model
        
    def forward(self, observation: th.Tensor) -> th.Tensor:
        """只导出actor网络的前向传播"""
        # 使用确定性策略
        with th.no_grad():
            # 这里我们只需要动作，不需要其他返回值
            action = self.model.predict(observation.numpy(), deterministic=True)[0]
            return th.as_tensor(action)

def export_model_to_onnx(model_path: str, output_dir: str) -> Tuple[str, Dict[str, Any]]:
    """
    将SBX模型导出为ONNX格式，并提取模型参数
    
    Args:
        model_path: SBX模型路径
        output_dir: 输出目录
        
    Returns:
        Tuple[str, Dict]: ONNX文件路径和模型参数字典
    """
    if not HAVE_SBX:
        raise ImportError("需要安装SBX、PyTorch和ONNX包才能导出模型")
    
    print(f"加载模型: {model_path}")
    model = SAC.load(model_path)
    
    # 提取模型参数
    model_params = {}
    
    # 获取观察空间和动作空间信息
    observation_space = model.observation_space
    action_space = model.action_space
    
    # 提取维度信息
    if hasattr(observation_space, 'shape'):
        model_params['input_dim'] = observation_space.shape[0]
    else:
        model_params['input_dim'] = 1
        
    if hasattr(action_space, 'shape'):
        model_params['output_dim'] = action_space.shape[0]
    else:
        model_params['output_dim'] = 1
    
    # 尝试提取网络架构信息
    try:
        if hasattr(model, 'actor') and hasattr(model.actor, 'latent_pi_net'):
            # 提取actor网络架构
            actor_net = model.actor.latent_pi_net
            net_arch = []
            for module in actor_net.modules():
                if isinstance(module, th.nn.Linear):
                    net_arch.append(module.out_features)
            model_params['net_arch'] = net_arch[:-1]  # 排除输出层
        else:
            # 默认网络架构
            model_params['net_arch'] = [256, 256]
    except Exception as e:
        print(f"提取网络架构失败: {e}")
        model_params['net_arch'] = [256, 256]  # 默认值
    
    # 创建ONNX导出包装器
    onnxable_model = OnnxableSACPolicy(model)
    
    # 创建示例输入
    observation_size = model.observation_space.shape
    dummy_input = th.zeros(1, *observation_size)
    
    # ONNX输出路径
    onnx_path = os.path.join(output_dir, 'sbx_actor.onnx')
    
    # 导出到ONNX
    print(f"导出模型到ONNX: {onnx_path}")
    th.onnx.export(
        onnxable_model,
        dummy_input,
        onnx_path,
        opset_version=14,
        input_names=["observation"],
        output_names=["action"],
        dynamic_axes={
            "observation": {0: "batch_size"},
            "action": {0: "batch_size"}
        }
    )
    
    # 验证ONNX模型
    onnx_model = onnx.load(onnx_path)
    onnx.checker.check_model(onnx_model)
    print("ONNX模型验证通过")
    
    # 从ONNX模型中提取额外信息
    for input_info in onnx_model.graph.input:
        if input_info.name == "observation":
            # 获取输入维度
            input_shape = []
            for dim in input_info.type.tensor_type.shape.dim:
                if dim.dim_param:  # 动态维度
                    input_shape.append(-1)
                else:
                    input_shape.append(dim.dim_value)
            input_shape = input_shape[1:]  # 去掉批次维度
            model_params['input_shape'] = input_shape
    
    for output_info in onnx_model.graph.output:
        if output_info.name == "action":
            # 获取输出维度
            output_shape = []
            for dim in output_info.type.tensor_type.shape.dim:
                if dim.dim_param:  # 动态维度
                    output_shape.append(-1)
                else:
                    output_shape.append(dim.dim_value)
            output_shape = output_shape[1:]  # 去掉批次维度
            model_params['output_shape'] = output_shape

    if hasattr(model.action_space, 'low') and hasattr(model.action_space, 'high'):
        model_params['action_space_low'] = model.action_space.low.tolist()
        model_params['action_space_high'] = model.action_space.high.tolist()

    # 保存模型参数到JSON文件
    params_path = os.path.join(output_dir, 'sbx_actor_params.json')
    with open(params_path, 'w') as f:
        json.dump(model_params, f, indent=4)
    print(f"模型参数已保存到: {params_path}")
    
    return onnx_path, model_params

def main():
    """主函数"""
    # 模型路径
    model_path = './py/best_model_sbx/best_model'  # 修改为您的模型路径
    
    # 输出路径
    output_dir = './py/matlab_export'
    os.makedirs(output_dir, exist_ok=True)
    
    try:
        # 导出为ONNX并提取参数
        onnx_path, model_params = export_model_to_onnx(model_path, output_dir)
        
        # 打印模型参数
        print("\n模型参数摘要:")
        print(f"输入维度: {model_params['input_dim']}")
        print(f"输出维度: {model_params['output_dim']}")
        print(f"网络架构: {model_params['net_arch']}")
        
        print("\n导出完成！")
        
    except ImportError as e:
        print(f"错误: {e}")
        print("请确保已安装SBX、PyTorch和ONNX包")
    except Exception as e:
        print(f"导出过程中出错: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    main() 