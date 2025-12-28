classdef NormalizedCoccoEnv < rl.env.MATLABEnvironment
    % 描述:
    % 一个环境封装器，用于对观测和奖励进行移动平均归一化。
    % 模仿 Python stable-baselines3 中 VecNormalize 的核心功能。
    
    properties
        env             % 底层的 CoccoEnv 实例
        
        obs_rms         % RunningMeanStd 对象用于观测
        ret_rms         % RunningMeanStd 对象用于回报
        
        returns (1,1) double = 0.0 % 当前回合的累计折现回报
        gamma   (1,1) double = 0.95 % 折现因子
        
        clip_obs (1,1) double = 999999
        clip_reward (1,1) double = 999999
        epsilon (1,1) double = 1e-8
        
        norm_obs_idx    % 需要归一化的观测维度索引
        is_training (1,1) logical = true % 是否更新统计数据
    end
    
    methods
        function this = NormalizedCoccoEnv(env, varargin)
            % 构造函数
            p = inputParser;
            addRequired(p, 'env');
            addParameter(p, 'gamma', 0.95);
            addParameter(p, 'norm_obs_idx', []);
            parse(p, env, varargin{:});

            % 继承观测和动作空间定义
            this@rl.env.MATLABEnvironment(env.getObservationInfo(), env.getActionInfo());
            
            this.env = env;
            this.gamma = p.Results.gamma;
            this.norm_obs_idx = p.Results.norm_obs_idx;
            
            % 初始化 RunningMeanStd 对象
            obs_dim = env.getObservationInfo().Dimension(1);
            this.obs_rms = RunningMeanStd([obs_dim, 1]);
            this.ret_rms = RunningMeanStd([1, 1]);
        end
        
        function [InitialObservation, LoggedSignal] = reset(this)
            [rawObs, LoggedSignal] = this.env.reset();
            
            if this.is_training
                this.obs_rms.update(rawObs);
            end
            this.returns = 0;
            
            InitialObservation = this.normalize_obs(rawObs);
        end
        
        function [NextObs, Reward, IsDone, LoggedSignal] = step(this, Action)
            [rawNextObs, rawReward, IsDone, LoggedSignal] = this.env.step(Action);

            if this.is_training
                % 更新回报统计
                this.returns = rawReward + this.gamma * this.returns;
                this.ret_rms.update(this.returns);
                if IsDone
                    this.returns = 0;
                end
                
                % 更新观测统计
                this.obs_rms.update(rawNextObs);
            end

            NextObs = this.normalize_obs(rawNextObs);
            Reward = this.normalize_reward(rawReward);
        end
        
        function obs_norm = normalize_obs(this, obs)
            obs_norm = (obs - this.obs_rms.Mean) ./ sqrt(this.obs_rms.Variance + this.epsilon);
            % 对不需要归一化的维度，恢复原值
            if ~isempty(this.norm_obs_idx)
                all_idx = 1:length(obs);
                unnorm_idx = setdiff(all_idx, this.norm_obs_idx);
                obs_norm(unnorm_idx) = obs(unnorm_idx);
            end
            obs_norm = max(-this.clip_obs, min(this.clip_obs, obs_norm));
        end

        function obs = unnormalize_obs(this, obs_norm)
             obs = (obs_norm .* sqrt(this.obs_rms.Variance + this.epsilon)) + this.obs_rms.Mean;
             % 对不需要归一化的维度，恢复原值
             if ~isempty(this.norm_obs_idx)
                all_idx = 1:length(obs_norm);
                unnorm_idx = setdiff(all_idx, this.norm_obs_idx);
                obs(unnorm_idx) = obs_norm(unnorm_idx);
            end
        end
        
        function reward_norm = normalize_reward(this, reward)
            reward_norm = reward ./ sqrt(this.ret_rms.Variance + this.epsilon);
            reward_norm = max(-this.clip_reward, min(this.clip_reward, reward_norm));
        end
        
        function saveStats(this, filepath)
            % obs_rms_data = this.obs_rms;
            % ret_rms_data = this.ret_rms;
            eval_env = this;
            save(filepath, 'eval_env');
        end
        
        function loadStats(this, filepath)
            data = load(filepath);
            this.obs_rms = data.obs_rms_data;
            this.ret_rms = data.ret_rms_data;
        end
    end
end